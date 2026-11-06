#!/usr/bin/env python3
"""
headlamp_road_illuminance.py
----------------------------
Road illuminance from vehicle headlamp photometry (IES LM-63 file,
Type A or Type B, H-V angles), including lamp mounting position, aim,
toe-in/out, vehicle pitch and roll.

Only numpy is needed; matplotlib is optional (--plot).

-------------------------------------------------------------------------------
COORDINATES (ISO 8855 vehicle frame, vehicle at rest on a flat road)
    x : forward along the road            (0 = lamp reference plane)
    y : to the LEFT of vehicle centreline (right side is negative y)
    z : up                                (road surface z = 0)

PHOTOMETRIC ANGLES (automotive convention, as seen from the driver)
    H : horizontal angle, positive to the RIGHT
    V : vertical angle,   positive UP
    H = V = 0 is the lamp's photometric reference axis (HV point).

    Two H-V systems exist. Choose the one your goniometer / export uses:

    --angle-system azel  (default)
        H measured in the horizontal plane first, then V as elevation.
        Lines of constant H are vertical planes.
        => direction = (cosV cosH, -cosV sinH, sinV)

    --angle-system elaz
        V is the tilt of a plane about the lamp's lateral axis, then H is
        measured inside that tilted plane.
        Lines of constant V are planes through the lateral axis.
        => direction = (cosH cosV, -sinH, cosH sinV)

    Near the HV point both agree; at large |H| they differ (e.g. at H = 40 deg,
    V = -2 deg the "elaz" plane tilt is -2.6 deg). This matters near the cut-off,
    so check against a known goniometer reading if unsure.

MOUNTING AND ORIENTATION
    --height          lamp reference-centre height above road [m]
    --separation      lateral distance between left and right lamp centres [m]
    --aim             vertical aim of the lamp axis in PERCENT (ECE style):
                      -1.0 = 1 % down (10 cm per 10 m). Negative = down.
    --toe             horizontal aim [deg], positive = OUTWARD
                      (left lamp turns left, right lamp turns right)
    --pitch           vehicle pitch [deg], positive = nose UP (e.g. rear load)
    --roll            vehicle roll  [deg], positive = left side UP
    --pivot-x         x of the pitch pivot (e.g. -2.7 for a rear-axle pivot
                      2.7 m behind the lamps). Pitch then also changes lamp height.

    Lamps:
      single IES file  -> used for both lamps
      --ies-right FILE -> separate file for the right lamp
      --mirror-right   -> right lamp uses the left file mirrored in H
      --only left|right -> evaluate one lamp only
      --mirror-all     -> mirror H for both lamps (RHT <-> LHT conversion)

ILLUMINANCE METHOD (--method)
    hspr (default) : HSPR / GTB reference method (TU Darmstadt HSPR software).
                     Illuminance on a vertical plane facing the vehicle at --eval-height,
                     E = I cos(alpha) / d^2. Aim and toe shift the photometric data by an
                     angle; H, V are read on the 0.05 deg grid. Thresholds count E > limit.
    road           : horizontal illuminance on the road surface, E = I h / d^3.

EVALUATION
    Grid of illuminance on the road, x in [x-min, x-max], y in [y-min, y-max].
    Reports:
      * E at a list of road points (--points "x,y;x,y;...")
      * Range: farthest distance on given lateral lines where E >= threshold
      * Width: lateral extent of the E >= threshold area at given distances
      * Geometric cut-off distance of the lamp axis (h / tan(aim))

EXAMPLE
    python headlamp_road_illuminance.py lowbeam.ies --height 0.65 \
        --separation 1.40 --aim -1.0 --pitch 0.5 --pivot-x -2.7 \
        --thresholds 1 3 --range-y 0 -3.5 3.5 --plot --csv road.csv
"""

import argparse
import math
import sys

import numpy as np


# =============================================================================
# IES reader for Type A / Type B (H-V) photometry
# =============================================================================
class HeadlampIES:
    def __init__(self, path, mirror=False):
        self.path = path
        self._parse(path)
        self._prepare(mirror)

    def _parse(self, path):
        with open(path, "r", errors="replace") as fh:
            lines = fh.read().splitlines()
        idx = 0
        while idx < len(lines) and not lines[idx].strip().upper().startswith("TILT"):
            idx += 1
        if idx == len(lines):
            raise ValueError(f"{path}: no TILT= line, not a valid IES file.")
        tilt = lines[idx].split("=", 1)[1].strip().upper()
        nums = [float(t) for t in " ".join(lines[idx + 1:]).replace(",", " ").split()]
        pos = 0
        if tilt == "INCLUDE":
            pos += 1
            n = int(nums[pos]); pos += 1 + 2 * n
        (self.n_lamps, self.lm_per_lamp, self.multiplier, n_v, n_h,
         self.ptype, _units, _w, _l, _h) = nums[pos:pos + 10]
        pos += 10
        self.ballast, _f, self.watts = nums[pos:pos + 3]
        pos += 3
        n_v, n_h = int(n_v), int(n_h)
        self.V = np.array(nums[pos:pos + n_v]); pos += n_v
        self.H = np.array(nums[pos:pos + n_h]); pos += n_h
        cd = np.array(nums[pos:pos + n_v * n_h])
        if cd.size != n_v * n_h:
            raise ValueError(f"{path}: candela table incomplete.")
        # rows = horizontal angles, columns = vertical angles
        self.cd = cd.reshape(n_h, n_v) * self.multiplier * self.ballast

        if int(self.ptype) == 1:
            raise ValueError(
                f"{path}: Type C photometry. Headlamp mode needs H-V data "
                "(Type A or B). Re-export the file with H-V angles.")

    def _prepare(self, mirror):
        H, V, cd = self.H, self.V, self.cd
        # Laterally symmetric file (H from 0 only): mirror to negative H
        if H[0] >= 0 and H[-1] > 0 and len(H) > 1:
            if H[0] == 0:
                H = np.concatenate([-H[:0:-1], H]); cd = np.vstack([cd[:0:-1], cd])
            else:
                H = np.concatenate([-H[::-1], H]); cd = np.vstack([cd[::-1], cd])
        if mirror:
            H = -H[::-1]; cd = cd[::-1]
        # make sure both axes ascend
        oh, ov = np.argsort(H), np.argsort(V)
        self.Hs, self.Vs, self.cds = H[oh], V[ov], cd[oh][:, ov]

    def intensity(self, H, V):
        """Bilinear interpolation, 0 outside the measured range."""
        H = np.asarray(H, float); V = np.asarray(V, float)
        hs, vs, cd = self.Hs, self.Vs, self.cds
        Hc = np.clip(H, hs[0], hs[-1]); Vc = np.clip(V, vs[0], vs[-1])
        i = np.clip(np.searchsorted(hs, Hc, side="right") - 1, 0, len(hs) - 2)
        j = np.clip(np.searchsorted(vs, Vc, side="right") - 1, 0, len(vs) - 2)
        th = (Hc - hs[i]) / (hs[i + 1] - hs[i])
        tv = (Vc - vs[j]) / (vs[j + 1] - vs[j])
        I = ((1 - th) * (1 - tv) * cd[i, j] + th * (1 - tv) * cd[i + 1, j]
             + (1 - th) * tv * cd[i, j + 1] + th * tv * cd[i + 1, j + 1])
        outside = (H < hs[0]) | (H > hs[-1]) | (V < vs[0]) | (V > vs[-1])
        return np.where(outside, 0.0, I)

    def flux(self, system, dh=0.25, dv=0.25):
        """Luminous flux within the measured H-V range [lm]."""
        H = np.arange(self.Hs[0], self.Hs[-1], dh) + dh / 2
        V = np.arange(self.Vs[0], self.Vs[-1], dv) + dv / 2
        HH, VV = np.meshgrid(H, V)
        I = self.intensity(HH, VV)
        # solid-angle element: cos(V) dH dV (azel) or cos(H) dH dV (elaz)
        w = np.cos(np.radians(VV if system == "azel" else HH))
        return float(np.sum(I * w) * math.radians(dh) * math.radians(dv))


# =============================================================================
# Geometry
# =============================================================================
def rot_x(a):   # roll, positive = +y (left) side up
    c, s = math.cos(a), math.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def rot_nose_up(a):  # pitch, positive = +x (nose) up
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])


def lamp_frame(aim_pct, yaw_deg):
    """Lamp axes in the vehicle frame: columns = forward, left, up."""
    el = math.atan(aim_pct / 100.0)          # aim in percent -> radians
    ps = math.radians(yaw_deg)               # positive = turned to the left
    f = np.array([math.cos(el) * math.cos(ps), math.cos(el) * math.sin(ps), math.sin(el)])
    l = np.array([-math.sin(ps), math.cos(ps), 0.0])
    u = np.cross(f, l)
    return np.column_stack([f, l, u])


def build_lamps(args):
    """Return list of dicts: name, position (3,), axes (3x3), ies."""
    left_ies = HeadlampIES(args.ies, mirror=args.mirror_all)
    if args.ies_right:
        right_ies = HeadlampIES(args.ies_right, mirror=args.mirror_all)
    elif args.mirror_right:
        right_ies = HeadlampIES(args.ies, mirror=not args.mirror_all)
    else:
        right_ies = left_ies

    R_veh = rot_x(math.radians(args.roll)) @ rot_nose_up(math.radians(args.pitch))
    pivot = np.array([args.pivot_x, 0.0, 0.0])

    lamps = []
    half = args.separation / 2
    spec = [("left", +half, +args.toe, left_ies), ("right", -half, -args.toe, right_ies)]
    for name, y, yaw, ies in spec:
        if args.only and args.only != name:
            continue
        p0 = np.array([0.0, y, args.height])
        pos = pivot + R_veh @ (p0 - pivot)
        hspr = args.method == "hspr"
        aim_pct = args.aim if args.aim_unit == "pct" else math.tan(math.radians(args.aim)) * 100
        aim_deg = args.aim if args.aim_unit == "deg" else math.degrees(math.atan(args.aim / 100))
        # HSPR shifts the photometric data by the aim/toe angles; road mode tilts the lamp axis
        axes = R_veh @ (lamp_frame(0.0, 0.0) if hspr else lamp_frame(aim_pct, yaw))
        lamps.append(dict(name=name, pos=pos, axes=axes, ies=ies, method=args.method,
                          zp=args.eval_height if hspr else 0.0,
                          aimV=aim_deg if hspr else 0.0, aimH=-yaw if hspr else 0.0))
    return lamps


def angles_HV(vec_local, system):
    """vec_local[..., 0:3] = (forward, left, up) components -> H, V in degrees."""
    f, l, u = vec_local[..., 0], vec_local[..., 1], vec_local[..., 2]
    if system == "azel":
        H = np.degrees(np.arctan2(-l, f))
        V = np.degrees(np.arctan2(u, np.hypot(f, l)))
    else:  # elaz
        V = np.degrees(np.arctan2(u, f))
        H = np.degrees(np.arctan2(-l, np.hypot(f, u)))
    return H, V


HSPR_STEP = 0.05


def road_illuminance(lamps, X, Y, system, scale):
    """Illuminance at road points (X, Y).

    method "hspr": HSPR / GTB reference method. Illuminance on a vertical plane facing
        the vehicle at the evaluation height, E = I cos(alpha) / d^2 with cos(alpha) =
        forward distance / d. H and V are snapped to the 0.05 deg grid and clamped to
        H +-45 deg, V +-15 deg; aim and toe shift the photometric data by an angle.
    method "road": horizontal illuminance on the road surface, E = I h / d^3.
    """
    E = np.zeros(np.broadcast(X, Y).shape)
    for lp in lamps:
        px, py, pz = lp["pos"]
        hspr = lp["method"] == "hspr"
        v = np.stack(np.broadcast_arrays(X - px, Y - py, np.full_like(E, lp["zp"] - pz)), axis=-1)
        D = np.linalg.norm(v, axis=-1)
        loc = v @ lp["axes"]                 # components along lamp f, l, u
        H, V = angles_HV(loc, system)
        if hspr:
            H = np.clip(np.round(H / HSPR_STEP) * HSPR_STEP, -45, 45)
            V = np.clip(np.round(V / HSPR_STEP) * HSPR_STEP, -15, 15)
            I = lp["ies"].intensity(H - lp["aimH"], V - lp["aimV"]) * scale
        else:
            I = lp["ies"].intensity(H, V) * scale
        # behind the lamp -> no light
        I = np.where(loc[..., 0] > 0, I, 0.0)
        E += I * (v[..., 0] if hspr else pz) / D ** 3
    return E


def meets(lamps, E, thr):
    """HSPR counts E > threshold; road mode counts E >= threshold."""
    return E > thr if lamps and lamps[0]["method"] == "hspr" else E >= thr


# =============================================================================
# Evaluation helpers
# =============================================================================
def farthest_range(lamps, y0, thr, x_max, system, scale, dx=0.05):
    xs = np.arange(dx, x_max + dx, dx)
    E = road_illuminance(lamps, xs, np.full_like(xs, y0), system, scale)
    hit = np.nonzero(meets(lamps, E, thr))[0]
    return None if hit.size == 0 else xs[hit[-1]]


def width_at(lamps, x0, thr, y_lim, system, scale, dy=0.02):
    ys = np.arange(-y_lim, y_lim + dy, dy)
    E = road_illuminance(lamps, np.full_like(ys, x0), ys, system, scale)
    hit = np.nonzero(meets(lamps, E, thr))[0]
    if hit.size == 0:
        return None
    return ys[hit[0]], ys[hit[-1]]           # (rightmost, leftmost) in y


def side(y):
    return "centre" if abs(y) < 1e-9 else (f"{abs(y):.2f} m {'left' if y > 0 else 'right'}")


# =============================================================================
# CLI
# =============================================================================
def main():
    p = argparse.ArgumentParser(
        description="Road illuminance from headlamp IES photometry.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("ies", help="IES file (left lamp, or both lamps)")
    p.add_argument("--ies-right", help="separate IES file for the right lamp")
    p.add_argument("--mirror-right", action="store_true",
                   help="right lamp = left file mirrored in H")
    p.add_argument("--mirror-all", action="store_true",
                   help="mirror H for all lamps (RHT <-> LHT)")
    p.add_argument("--only", choices=["left", "right"], help="evaluate one lamp only")
    p.add_argument("--angle-system", choices=["azel", "elaz"], default="azel")
    p.add_argument("--method", choices=["hspr", "road"], default="hspr",
                   help="hspr = GTB/HSPR vertical-plane illuminance (default); "
                        "road = horizontal illuminance on the road surface")
    p.add_argument("--eval-height", type=float, default=0.0,
                   help="HSPR evaluation height above the road [m] (zones A, B, D: 0; C, E: 0.25)")

    g = p.add_argument_group("mounting")
    g.add_argument("--height", type=float, required=True, help="lamp centre height [m]")
    g.add_argument("--separation", type=float, default=1.4, help="lamp centre distance [m]")
    g.add_argument("--aim", type=float, default=0.0, help="vertical aim, negative = down (unit: --aim-unit)")
    g.add_argument("--aim-unit", choices=["pct", "deg"], default="pct",
                   help="unit of --aim: pct (ECE, -1.0 = 1 %% down) or deg (HSPR)")
    g.add_argument("--toe", type=float, default=0.0, help="horizontal aim [deg], + = outward")
    g.add_argument("--pitch", type=float, default=0.0, help="vehicle pitch [deg], + = nose up")
    g.add_argument("--roll", type=float, default=0.0, help="vehicle roll [deg], + = left side up")
    g.add_argument("--pivot-x", type=float, default=0.0, help="pitch pivot x [m] (rear axle < 0)")
    g.add_argument("--scale", type=float, default=1.0,
                   help="intensity factor (voltage, ageing, tolerance)")
    g.add_argument("--lumens", type=float, help="rescale each lamp to this flux [lm]")

    e = p.add_argument_group("evaluation")
    e.add_argument("--x-min", type=float, default=0.0)
    e.add_argument("--x-max", type=float, default=150.0)
    e.add_argument("--y-min", type=float, default=-15.0, help="right limit (negative)")
    e.add_argument("--y-max", type=float, default=15.0, help="left limit (positive)")
    e.add_argument("--step-x", type=float, default=0.5)
    e.add_argument("--step-y", type=float, default=0.25)
    e.add_argument("--thresholds", type=float, nargs="+", default=[1.0, 3.0], help="lx")
    e.add_argument("--range-y", type=float, nargs="+", default=[0.0, -3.5, 3.5],
                   help="lateral lines for range evaluation [m]")
    e.add_argument("--width-x", type=float, nargs="+", default=[10, 20, 30, 50],
                   help="distances for width evaluation [m]")
    e.add_argument("--points", default=None, help='road points "x,y;x,y;..." [m]')
    e.add_argument("--road-edges", type=float, nargs=2, metavar=("RIGHT", "LEFT"),
                   help="road edge y positions for the plot, e.g. -1.75 5.25")
    e.add_argument("--csv", help="save grid to CSV")
    e.add_argument("--plot", action="store_true", help="isolux plot (matplotlib)")
    args = p.parse_args()

    lamps = build_lamps(args)
    sysname = args.angle_system

    # optional flux rescale (per lamp file)
    scale = args.scale
    flux = lamps[0]["ies"].flux(sysname)
    if args.lumens:
        scale *= args.lumens / flux

    # ---------------- header
    print(f"\nIES           : {args.ies}" + (f" | right: {args.ies_right}" if args.ies_right else ""))
    print(f"Angle system  : {sysname}  (H + right, V + up)")
    print("Method        : " + (f"HSPR, vertical-plane illuminance at {args.eval_height} m height"
                               if args.method == "hspr" else "road surface, horizontal illuminance"))
    print(f"Flux in file  : {flux:.0f} lm (within measured H-V range)"
          + (f" -> scaled to {args.lumens:.0f} lm" if args.lumens else ""))
    print(f"Mounting      : h = {args.height} m, separation = {args.separation} m, "
          f"aim = {args.aim:+.2f} {'%' if args.aim_unit == 'pct' else 'deg'}, toe = {args.toe:+.2f} deg")
    print(f"Vehicle       : pitch = {args.pitch:+.2f} deg, roll = {args.roll:+.2f} deg, "
          f"pivot x = {args.pivot_x} m")
    for lp in lamps:
        x, y, z = lp["pos"]
        f = lp["axes"][:, 0]
        el = math.degrees(math.asin(f[2])) + lp["aimV"]
        az = math.degrees(math.atan2(f[1], f[0]))
        hit = (f"axis meets road at {z / math.tan(-math.radians(el)):.1f} m"
               if el < 0 else "axis does not meet road (aimed level/up)")
        print(f"  {lp['name']:5s} lamp : pos ({x:+.3f}, {y:+.3f}, {z:.3f}) m, "
              f"axis elev {el:+.3f} deg, azim {az:+.3f} deg (+left) -> {hit}")

    # ---------------- points
    if args.points:
        print("\nPoint illuminance:")
        for pt in args.points.split(";"):
            x, y = map(float, pt.split(","))
            Ep = road_illuminance(lamps, np.array([x]), np.array([y]), sysname, scale)[0]
            print(f"  x = {x:6.1f} m, y = {y:+6.2f} m ({side(y)}): E = {Ep:8.2f} lx")

    # ---------------- range and width
    print("\nRange (farthest distance where E reaches the threshold):")
    for thr in args.thresholds:
        for y0 in args.range_y:
            r = farthest_range(lamps, y0, thr, args.x_max, sysname, scale)
            rs = f"{r:6.1f} m" if r is not None else "   not reached"
            print(f"  {thr:5.1f} lx on line {side(y0):>14s}: {rs}")

    print("\nWidth of illuminated area (E reaches the threshold):")
    ylim = max(abs(args.y_min), abs(args.y_max))
    for thr in args.thresholds:
        for x0 in args.width_x:
            w = width_at(lamps, x0, thr, ylim, sysname, scale)
            if w is None:
                print(f"  {thr:5.1f} lx at {x0:5.1f} m: none")
            else:
                yr, yl = w
                print(f"  {thr:5.1f} lx at {x0:5.1f} m: from {side(yr)} to {side(yl)} "
                      f"(width {yl - yr:.2f} m)")

    # ---------------- grid
    xs = np.arange(args.x_min, args.x_max + 1e-9, args.step_x)
    ys = np.arange(args.y_min, args.y_max + 1e-9, args.step_y)
    X, Y = np.meshgrid(xs, ys)
    E = road_illuminance(lamps, X, Y, sysname, scale)
    k = np.unravel_index(np.argmax(E), E.shape)
    print(f"\nGrid {len(xs)} x {len(ys)}: E_max = {E.max():.1f} lx at "
          f"x = {X[k]:.1f} m, y = {Y[k]:+.2f} m")

    if args.csv:
        with open(args.csv, "w") as fh:
            fh.write("y\\x," + ",".join(f"{x:.2f}" for x in xs) + "\n")
            for j in range(len(ys)):
                fh.write(f"{ys[j]:.2f}," + ",".join(f"{v:.3f}" for v in E[j]) + "\n")
        print(f"Grid saved to {args.csv}")

    if args.plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(12, 5))
        levels = [0.5, 1, 3, 5, 10, 20, 50, 100]
        levels = [l for l in levels if l < E.max()]
        cf = ax.contourf(X, Y, np.clip(E, 0.1, None), levels=[0.1] + levels + [E.max()],
                         cmap="magma")
        cs = ax.contour(X, Y, E, levels=levels, colors="w", linewidths=0.6)
        ax.clabel(cs, fmt="%g lx", fontsize=8)
        fig.colorbar(cf, label="E [lx]")
        for lp in lamps:
            ax.plot(lp["pos"][0], lp["pos"][1], "c>", ms=8)
        if args.road_edges:
            for yy in args.road_edges:
                ax.axhline(yy, color="0.7", ls="--", lw=1)
        ax.axhline(0, color="0.5", ls=":", lw=0.8)
        ax.set_xlabel("x, distance ahead [m]")
        ax.set_ylabel("y [m]  (+ left / - right)")
        ax.set_title(f"Headlamp road isolux | h {args.height} m, aim {args.aim:+.1f} %, "
                     f"pitch {args.pitch:+.1f} deg")
        plt.tight_layout(); plt.show()


if __name__ == "__main__":
    main()
