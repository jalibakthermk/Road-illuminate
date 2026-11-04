# Headlamp Road Illuminance

`headlamp_road_illuminance.py` calculates the horizontal illuminance on a flat road surface produced by a vehicle's headlamps. It starts from the headlamp's intensity distribution, given as an IES (LM-63) file in H-V angles.

It takes into account the real installation of the lamps on the vehicle:

- mounting height and the lateral distance between the lamps
- vertical aim and horizontal toe
- vehicle pitch and roll, for example from load
- separate or mirrored files for the left and right lamps

**Inputs:** one or two IES files plus the mounting and evaluation options on the command line.

**Outputs:** the lamp axis geometry, illuminance at chosen road points, range to an illuminance threshold, width of the illuminated area, the maximum illuminance, and an optional CSV grid and isolux plot.

Typical uses are comparing beam patterns, checking how aim or load changes the road light, and quick design iterations without a full simulation tool.

## Contents

1. [Requirements](#1-requirements)
2. [Quick start](#2-quick-start)
3. [Coordinate system and sign conventions](#3-coordinate-system-and-sign-conventions)
4. [Photometric angle systems](#4-photometric-angle-systems)
5. [Mounting and vehicle attitude](#5-mounting-and-vehicle-attitude)
6. [Lamp file options](#6-lamp-file-options)
7. [Evaluation options](#7-evaluation-options)
8. [Understanding the output](#8-understanding-the-output)
9. [Worked examples](#9-worked-examples)
10. [Calculation method](#10-calculation-method)
11. [Validation, limitations and troubleshooting](#11-validation-limitations-and-troubleshooting)
12. [Command-line reference](#12-command-line-reference)

## 1. Requirements

- Python 3.8 or newer
- numpy (required)
- matplotlib (optional, only for `--plot`)

```
pip install -r requirements.txt
```

No further installation is needed. Run the script from a terminal in its folder.

The IES file must use H-V photometry (Type A or Type B). Type C files are rejected with a message asking you to re-export them, because Type C angles describe a downward-pointing luminaire, not a forward beam.

## 2. Quick start

The minimum command needs only the IES file and the lamp height:

```
python headlamp_road_illuminance.py lowbeam.ies --height 0.65
```

A more complete run for a loaded vehicle with a plot and a CSV export:

```
python headlamp_road_illuminance.py lowbeam.ies --height 0.65 \
    --separation 1.40 --aim -1.0 --pitch 0.5 --pivot-x -2.7 \
    --thresholds 1 3 --range-y 0 -3.5 3.5 --plot --csv road.csv
```

This places two lamps 0.65 m high and 1.40 m apart, aimed 1 % down. The vehicle nose is raised 0.5° about a rear axle 2.7 m behind the lamps. The script reports the 1 lx and 3 lx range on three lines and saves the grid. Run `python headlamp_road_illuminance.py -h` to see every option.

## 3. Coordinate system and sign conventions

The road uses the ISO 8855 vehicle frame with the vehicle at rest on a flat road.

| Axis | Direction | Zero point |
| --- | --- | --- |
| x | Forward along the road | Lamp reference plane (front of the vehicle) |
| y | To the left of the vehicle centreline; right is negative | Vehicle centreline |
| z | Up | Road surface |

Photometric angles follow the automotive convention, as seen from the driver's seat:

| Angle | Positive direction | Zero |
| --- | --- | --- |
| H | Right | Lamp reference axis |
| V | Up | Lamp reference axis |

The two systems run in opposite lateral directions: +y is left on the road, but +H is right in the photometry. The script handles the conversion. Enter road positions in y (for example `-3.5` for 3.5 m right) and read IES data in H as measured.

## 4. Photometric angle systems

An H-V file can be measured in two ways. Tell the script which one your goniometer or optical software uses, with `--angle-system`.

| Option | How the direction is built | Lines of constant angle | Direction vector (forward, left, up) |
| --- | --- | --- | --- |
| `azel` (default) | H is turned in the horizontal plane first, then V is raised as elevation | Constant H = vertical planes | (cosV·cosH, −cosV·sinH, sinV) |
| `elaz` | V tilts a plane about the lamp's lateral axis, then H is measured inside that plane | Constant V = planes through the lateral axis | (cosH·cosV, −sinH, cosH·sinV) |

Both systems agree on the H = 0 and V = 0 lines. They separate as |H| grows. At H = 40°, V = −2° (azel), the elaz plane tilt is −2.6°. Near a sharp cut-off this shifts the light on the road noticeably.

**How to choose:**

1. Check the goniometer manual or the export settings of your simulation tool.
2. If unsure, take one measured value at a large H just below the cut-off. Run both systems and keep the one whose road result matches the screen measurement.

The mapping of the IES labels "Type A" and "Type B" to these two systems is not consistent between tools, so the script does not guess it from the file.

## 5. Mounting and vehicle attitude

| Option | Unit | Default | Meaning | Typical values |
| --- | --- | --- | --- | --- |
| `--height` | m | required | Height of the lamp reference centre above the road | 0.50–1.20 (cars), up to 1.5 (trucks) |
| `--separation` | m | 1.40 | Lateral distance between the left and right lamp centres | 1.2–1.6 |
| `--aim` | % | 0.0 | Vertical aim of the lamp axis; negative = down. −1.0 % = 10 cm drop per 10 m | −0.5 to −2.5 for low beam |
| `--toe` | deg | 0.0 | Horizontal aim; positive = outward (left lamp turns left, right lamp turns right) | −2 to +2 |
| `--pitch` | deg | 0.0 | Vehicle pitch; positive = nose up, as with rear load | −1 to +2 |
| `--roll` | deg | 0.0 | Vehicle roll; positive = left side up | −3 to +3 |
| `--pivot-x` | m | 0.0 | x position of the pitch pivot. 0 = pitch rotates the lamps in place; −2.7 = pivot at a rear axle 2.7 m back, so nose-up pitch also raises the lamps | −2.5 to −3.0 for the rear axle |
| `--scale` | – | 1.0 | Multiplies all intensities, for example for supply voltage, ageing or tolerance | 0.8–1.1 |
| `--lumens` | lm | none | Rescales each lamp file so its flux equals this value | – |

**Order of rotations.** Each lamp is first aimed (aim, then toe) in the vehicle frame. The whole vehicle is then pitched about `--pivot-x` and rolled about the x axis. Both the lamp positions and the lamp axes move with the vehicle.

**Aim versus pitch.** Aim is the lamp's adjustment relative to the vehicle. Pitch is the vehicle's change relative to the road. They add on the road: an aim of −1 % (−0.57°) with a nose-up pitch of +0.5° leaves the axis only 0.07° below horizontal.

## 6. Lamp file options

| Case | Options | What happens |
| --- | --- | --- |
| One file for both lamps | `lamp.ies` only | The same pattern is placed at the left and right positions |
| Different left and right files | `left.ies --ies-right right.ies` | Each lamp uses its own measurement |
| Right lamp is a mirror image of the left | `left.ies --mirror-right` | The right lamp uses the left file with H reversed |
| Evaluate one lamp | `--only left` or `--only right` | The other lamp is ignored |
| Convert between traffic sides | `--mirror-all` | H is reversed for all lamps; a right-hand-traffic beam becomes a left-hand-traffic beam |

If a file covers only H ≥ 0 (laterally symmetric beam, for example a high beam), the script mirrors it to negative H automatically.

## 7. Evaluation options

| Option | Default | Meaning |
| --- | --- | --- |
| `--x-min`, `--x-max` | 0, 150 m | Distance range of the road grid ahead of the lamps |
| `--y-min`, `--y-max` | −15, +15 m | Lateral range of the grid; negative = right |
| `--step-x`, `--step-y` | 0.5, 0.25 m | Grid spacing |
| `--thresholds` | 1 3 lx | Illuminance levels used for range and width |
| `--range-y` | 0 −3.5 3.5 m | Lateral lines along which the range is found |
| `--width-x` | 10 20 30 50 m | Distances at which the lit width is measured |
| `--points` | none | Individual road points, written as `"x,y;x,y"` |
| `--road-edges` | none | Two y values (right, left) drawn as dashed road edges in the plot |
| `--csv` | none | File name for the grid export |
| `--plot` | off | Shows the isolux plot |

**Range** is the farthest distance along a line y = constant where E is still at or above the threshold. It is found on a 5 cm step, so it is accurate to about ±0.05 m.

**Width** is measured across the road at one distance, on a 2 cm step. It gives the outermost right and left points that reach the threshold. If the lit area has a gap in the middle, the width still spans the outermost points.

## 8. Understanding the output

The report is printed in this order.

1. **Header:** the files, the angle system, the flux found in the file, and the mounting values.
2. **Lamp lines:** one line per lamp. Each gives the lamp position after pitch and roll, the axis elevation and azimuth, and where the axis meets the road.
3. **Point illuminance:** E at each `--points` location.
4. **Range:** one line per threshold and lateral line. "not reached" means E never reaches the threshold on that line.
5. **Width:** the right and left limits and the total width at each `--width-x` distance.
6. **Grid summary:** the maximum illuminance and where it occurs, plus the CSV file name if requested.

Example lamp line:

```
left  lamp : pos (-0.006, +0.700, 0.674) m, axis elev -0.073 deg, azim +0.000 deg (+left) -> axis meets road at 529.1 m
```

Here the 0.5° nose-up pitch has almost cancelled the −1 % aim, so the axis meets the road very far ahead. A result like this is a warning sign of glare under load.

**Flux in file** is integrated over the measured H-V range only. Light outside that range, for example far below the lamp, is not included.

**CSV layout:** the first row holds the x values, the first column holds the y values, and each cell is E in lx.

**Plot:** a top view with distance ahead on the horizontal axis and y on the vertical axis (left is up). Labelled isolux lines are drawn at 0.5, 1, 3, 5, 10, 20, 50 and 100 lx. Arrows mark the lamps, and the dotted line is the vehicle centreline.

## 9. Worked examples

**A. Low beam, unloaded car, nominal aim**

```
python headlamp_road_illuminance.py lowbeam.ies --height 0.65 --separation 1.40 --aim -1.0 --plot
```

The axis meets the road at 65 m (0.65 m ÷ 0.01). Check the 1 lx range on the right edge line and the width at 20–30 m.

**B. Same car with a full rear load**

```
python headlamp_road_illuminance.py lowbeam.ies --height 0.65 --aim -1.0 --pitch 0.8 --pivot-x -2.7
```

Compare the range with example A. A large increase on the left line (oncoming lane) points to glare risk. That is the case for automatic or manual levelling.

**C. Different left and right lamps, one lamp at a time**

```
python headlamp_road_illuminance.py L.ies --ies-right R.ies --height 0.70 --only left --csv left.csv
python headlamp_road_illuminance.py L.ies --ies-right R.ies --height 0.70 --only right --csv right.csv
```

Use this to see each lamp's share of the combined pattern.

**D. Convert a right-hand-traffic beam for left-hand traffic**

```
python headlamp_road_illuminance.py lowbeam_RHT.ies --height 0.65 --aim -1.0 --mirror-all
```

**E. Specific road points and a custom lane layout**

```
python headlamp_road_illuminance.py lowbeam.ies --height 0.65 --aim -1.0 \
    --points "25,0;50,-1.75;75,-3.5" --range-y -1.75 1.75 5.25 --road-edges -1.75 5.25 --plot
```

This is a 3.5 m own lane with the vehicle centred in it (right edge at −1.75 m) and a 3.5 m oncoming lane on the left.

## 10. Calculation method

For every road point and every lamp, the script carries out five steps.

1. **Vector to the point.** v = P − L, from the lamp position L to the road point P = (x, y, 0).
2. **Lamp frame.** v is projected on the lamp's forward, left and up axes, which already include aim, toe, pitch and roll.
3. **Angles.** H and V are found from the local components with the chosen angle system (section 4). Points behind the lamp get zero intensity.
4. **Intensity.** I(H, V) is interpolated bilinearly from the candela table and multiplied by `--scale` and the flux factor. Directions outside the measured range give zero.
5. **Illuminance.** The horizontal illuminance from one lamp is:

```math
E = \frac{I(H,V)\,\cos\theta}{D^{2}} = \frac{I(H,V)\,h_L}{D^{3}}
```

Here D is the distance from lamp to point, h_L is the lamp height after pitch and roll, and θ is the angle of incidence measured from the road normal. The contributions of all lamps are added.

The lamp axes are built as:

```math
\mathbf{f} = (\cos\varepsilon\cos\psi,\; \cos\varepsilon\sin\psi,\; \sin\varepsilon), \quad \varepsilon = \arctan(\mathrm{aim}/100)
```

In this formula ψ is the toe angle (positive to the left in vehicle coordinates). The left axis is l = (−sin ψ, cos ψ, 0) and the up axis is u = f × l. The vehicle rotation R = R_roll · R_pitch is then applied to the axes, and to the lamp positions about the pivot.

The flux in the file is found by integrating I over the measured range. The solid-angle element is cos V dH dV for azel and cos H dH dV for elaz.

## 11. Validation, limitations and troubleshooting

**Validation done.** The script was tested with a synthetic low-beam file (Gaussian beam peaking 5° right and 1° down).

| Check | Result |
| --- | --- |
| Integrated flux vs. fine reference integral | 478.19 vs 478.19 lm (< 0.01 %) |
| Point illuminance vs. hand calculation (20 m ahead, 3 m right) | 1.015 vs 1.023 lx; the 0.8 % gap comes from interpolating the 1° × 0.5° table |
| Aim −1 % at h = 0.65 m | Axis meets road at 65.0 m |
| Right-peaked beam | 0.90 lx on the right vs 0.30 lx on the left at 30 m; reversed with `--mirror-all` |
| Pitch +1° about pivot −2.7 m | Lamp height 0.697 m (expected 0.697 m) |
| Toe +2°, roll +2° | Left lamp turned left and raised, right lamp turned right and lowered |

**Limitations**

- The road is flat and level. There is no slope, camber or crest.
- Only horizontal illuminance on the road is calculated. Road luminance (with an r-table), vertical illuminance on objects and glare at a driver's eye are not.
- Bilinear interpolation follows the resolution of the file. Use a fine V step (0.1–0.25°) near the cut-off for accurate range values.
- Light outside the measured H-V range is treated as zero.
- TILT=INCLUDE data is skipped, not applied.

**Troubleshooting**

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| "Type C photometry" error | File exported as a luminaire (Type C) | Re-export with H-V angles |
| Light appears on the wrong side | File uses a different H sign | Add `--mirror-all` |
| Cut-off on the road is in the wrong place at large H | Wrong angle system | Try `--angle-system elaz` |
| Range is "not reached" everywhere | Aim too low, threshold too high, or units of the file wrong | Check `--aim`, the multiplier in the file and `--scale` |
| Very long range on the left line | Too little down aim, or pitch cancelling aim | Check the lamp line in the header |
| `--plot` fails | matplotlib not installed | `pip install matplotlib` |

## 12. Command-line reference

| Option | Value | Default | Purpose |
| --- | --- | --- | --- |
| `ies` | file | required | IES file for the left lamp, or both lamps |
| `--ies-right` | file | – | Separate file for the right lamp |
| `--mirror-right` | flag | off | Right lamp = left file mirrored in H |
| `--mirror-all` | flag | off | Mirror H for all lamps (traffic side) |
| `--only` | left / right | both | Evaluate one lamp |
| `--angle-system` | azel / elaz | azel | H-V convention of the file |
| `--height` | m | required | Lamp centre height |
| `--separation` | m | 1.40 | Lamp centre distance |
| `--aim` | % | 0.0 | Vertical aim, negative = down |
| `--toe` | deg | 0.0 | Horizontal aim, positive = outward |
| `--pitch` | deg | 0.0 | Vehicle pitch, positive = nose up |
| `--roll` | deg | 0.0 | Vehicle roll, positive = left up |
| `--pivot-x` | m | 0.0 | Pitch pivot position |
| `--scale` | factor | 1.0 | Intensity multiplier |
| `--lumens` | lm | – | Rescale each lamp to this flux |
| `--x-min`, `--x-max` | m | 0, 150 | Grid distance range |
| `--y-min`, `--y-max` | m | −15, 15 | Grid lateral range |
| `--step-x`, `--step-y` | m | 0.5, 0.25 | Grid spacing |
| `--thresholds` | lx … | 1 3 | Levels for range and width |
| `--range-y` | m … | 0 −3.5 3.5 | Lines for range |
| `--width-x` | m … | 10 20 30 50 | Distances for width |
| `--points` | "x,y;…" | – | Road points to report |
| `--road-edges` | m m | – | Road edges drawn in the plot |
| `--csv` | file | – | Save the grid |
| `--plot` | flag | off | Show the isolux plot |

## License

MIT – see [LICENSE](LICENSE).
