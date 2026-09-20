# Four exterior zipper families

This figure set displays **DD (1,0), DD (1,1), DO (1,0), and DO (1,1)**. It shows only certified stable parameter cells. Empty space has no classification. The reflected DD (0,1) presentation is omitted from the display; its independently certified cells contribute to DD (1,0) through the proved symmetry.

The three shared probes recover the exact parameters from the earlier DD (1,1) poster strip:

| Probe | Endpoint parameter |
|---|---|
| A | `0.5 + 0.2i` |
| B | `0.5 + 0.4i` |
| C | `0.34090875 + 0.43484625i` |

The new figure compares these same three probes across the four selected presentations. The probe markers identify parameters; they do not themselves assert stable membership. A and B have fresh positive box certificates in all four families. C has positive certificates in DD (1,1) and DO (1,1), and an exact secondary-contact certificate in DO (1,0). No membership claim is attached to the DD (1,0) curve at C.

## Outputs

| File | Intended use |
|---|---|
| `four_family_projection.pdf` | 15.6 × 10.2 cm paper overview: two DD/DO views |
| `four_family_DD.pdf` | 15.6 × 14.5 cm exterior slices and six DD specimens |
| `four_family_DO.pdf` | 15.6 × 14.5 cm exterior slices and six DO specimens |
| `four_family_gallery_paper.pdf` | Combined two-page gallery |
| `four_family_poster.pdf` | 60 × 42 cm poster plate, four rows with 3D placement, exterior slice, and three curves |
| `DD11_island_exterior.pdf` | Enlarged exterior image of the small certified island neighborhood |

PNG companions are previews. Labels remain vector text in the PDFs. Dense scientific meshes and polylines are rasterized at 600 dpi in the paper PDFs and 300 dpi in the poster plate.

## Geometry and certificates

The full four-dimensional atlas is described by the canonical contraction coefficients `(a,b)`, or exterior coefficients `(c1,c2)=(1/a,1/b)` with both moduli greater than one. The display is

```
theta1 = arg(a) = -arg(c1)
theta2 = arg(b) = -arg(c2)
rho = 2 / (abs(a) + abs(b))
w = abs(a) / (abs(a) + abs(b))
X = (4 + rho*cos(theta1))*cos(theta2)
Y = (4 + rho*cos(theta1))*sin(theta2)
Z = rho*sin(theta1)
```

At a fixed weight `w`, this is a one-to-one torus-coordinate embedding for `0 < rho < 4`. When weights are displayed together, `w` is the omitted coordinate and is encoded by color. Projected coincidences need not be four-dimensional intersections. The hole of the toroidal display is a coordinate feature.

The faint reference grid has `rho=2`. Only its half with `theta1` between `pi/2` and `3*pi/2` is drawn, because all four selected canonical charts have `a=-q` with positive real part of `q` in the shown stable region. This reference grid is not a stability classification.

`projection_geometry.py` includes the exact phase correction needed to normalize the DO endpoint maps to canonical translations. In particular, the unsigned exterior slice coordinate `c=1/q` is different from the canonical first coefficient `c1=-1/q`; the canonical opposite coefficient also has a phase correction.

Original positive masks are read from the sibling `zippers/` directory and reconstructed from the outward interval CSV records, including the reflected DD companion. Each 3D mesh face corresponds to a 2 × 2 block for which all four original closed cells have positive certificates. Both the original and conjugate cells are shown. The additional DD (1,1) island square and its conjugate are included in the mesh even though they are subpixel at the overview scale. Its exterior inset displays all 4,096 positive cells separately.

The 2D exterior plots use all positive original cells. They show the finite window `0.96 <= Re(c) <= 4.65`, `-2.70 <= Im(c) <= 0.13`; the complete slice is unbounded as `q` approaches zero. The illustrated upper half of the q-plane maps to the lower half of the c-plane. Reflection supplies the conjugate half.

Planar quadrilateral faces and pcolormesh corner chords are **display approximations** of the nonlinear images of certified q-cells. Certification belongs to the original closed q-cells. Neither coarsening nor projection is a new stability test. No background layer labels unproved cells, necessary-bound disks, or omitted presentations.

## Curve approximation

A and B use depth-20 polygonal approximants. The generic uniform bound is

```
r^20 * abs(q - 0.5) / (1-r),  r=max(abs(q),abs(1-q)).
```

It is less than `1.823282e-6` at A and `1.492703e-4` at B. C uses family-specific adaptive cylinder subdivision, with exact product stopping and a sharpened global chord bound derived from an exact depth-14 reference polygon. Every C curve has uniform geometric approximation error less than `5e-4`; exact bounds and source-partition audits are in `specimen_checks/adaptive_curve_metadata.json`.

The adaptive NPZ files retain the exact source traversal's `first_level_join_index`. That index, rather than proximity to the spatial join, determines red/blue first-piece coloring. No geometric simplification is used.

## Rebuild

Run from any directory:

```bash
python build_four_family_figures.py
```

Requires Python, NumPy, Matplotlib, Pillow, and pypdf. The sibling original masks/CSV records and `specimen_checks` adaptive NPZ files must remain in place. `--only projection`, `--only paper`, `--only poster`, and `--only island` rebuild individual groups. `--preview` uses lower-depth A/B curves and must not be used for the final release. All published outputs in this revision were rebuilt without that flag.
