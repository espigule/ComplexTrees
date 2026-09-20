# Reproducing the twelve marked zipper panels

The mathematical statements and exact return identities are in
`../manuscript/math/stable_zipper_slices.tex`. The OO01/OO10 proof is in `../manuscript/math/oo01_oo10_exact_disk.tex`.
The original seven-return verification records are retained.

## Files and meaning

- `zipper_stability_atlas.pdf`: 60 × 32 cm poster plate with twelve parameter
  panels and twelve linked canonical-curve examples.
- `zipper_stability_atlas_paper.pdf`: 15.6 × 12.8 cm manuscript parameter plate.
- `zipper_curve_gallery.pdf`: 15.6 × 13.3 cm manuscript curve plate.
- `zipper_stability_metadata.json`: plot meaning, curve construction and
  analytic error bound.
- `zipper_stability_run_manifest.json`: exact commands, compiler, assumptions,
  interval records, masks, byte counts and SHA-256 digests.

The horizontal parameter axis is Re(q), and the vertical parameter axis is
Im(q). Each panel shows the same upper-half window. Complex conjugation gives
the lower half of the same marked family. DD, DO and OO specify the two plane
orientations; the columns specify the endpoint-order bits e=(e1,e2). These are
twelve marked presentations, not a claim of twelve inequivalent unmarked
attractor classes.

Blue is a proved subset of the stable set. For five presentations each blue
grid cell represents an entire closed dyadic rectangle certified by outward
interval bounds. For OO00, OO01 and OO10 the full open diameter disks are
established exactly by invariant-triangle theorems. Every real parameter 0<q<1 is stable in all
twelve presentations.

Gray is unresolved, including the interiors of the four presentations not
covered by the seven-row return table. The boundary between blue and gray is
a boundary of this certificate cover and must not be called the true stable
set boundary. Pink is the part of the strict contraction lens outside or on
the common necessary diameter disk; it is proved nonembedded. White is
outside the strict contraction lens and is not classified as an attractor.
All valid zipper parameters give connected attractors, including pink ones.

The round parameter dot marks q*=3/8+i/4. Each curve is constructed from that
same exact parameter and the row/column's own maps, using sixteen zipper
refinements (65,537 vertices). The red and blue pieces join at q*. The exact
uniform polygonal approximation error is bounded by

    (sqrt(29)/8)^16 (sqrt(5)/8) / (1 - sqrt(29)/8) < 0.00152.

The geometric display is a polygonal approximation. Its stability labels
come from five independent point certificate runs and three exact OO
disk theorems, not from whether the plotted polyline appears simple.

## Running the verifier

The implementation uses IEEE 754 binary64 operations with one `nextafter`
step of outward widening for each elementary interval bound. It assumes
round-to-nearest arithmetic, a correctly rounded square root and gradual
underflow. Do not enable fast-math or floating-point contraction.

```sh
g++ -O2 -std=c++17 -fno-fast-math -ffp-contract=off \
  certify_zipper_cells.cpp -o certify_zipper_cells
./certify_zipper_cells DD11 DD11_cells.csv 512 12000 40
python build_zipper_stability_figures.py
```

The manifest supplies the seven full-grid commands and the point/rectangle
commands. The main grid is 512 by 256 on [0,1] × [0,1/2], with cell side
1/512. The limits are 12,000 pair visits and maximum word depth 40. Hitting
either limit is unresolved. Only cells that finish with every required
pair strictly separated are assigned status 1.

The companion pair DD01/DD10 is completed using the proved symmetry
q -> 1-conjugate(q). Its displayed masks unite direct and reflected companion
certificates. The earlier OO01/OO10 masks are retained with their reflection
completion, but these two displayed panels now use their exact full-disk
classification instead of the numerical inner covers.

CSV `minimum_squared_gap` is the smallest outward lower separation margin
among discarded leaf disk pairs. It is **not** a positive distance between
the original first-level attractor pieces, whose common endpoint is q.

## The wiggle-island point

The star L in DD11 is the exact rational parameter

    L = 0.34090875 + 0.43484625 i.

Calegari's *Wiggle Island*, arXiv:2205.11442v2, establishes that its stable
component is separate from the real-interval component. This topological
claim is cited to that theorem; it does not follow merely from our common
grid or its unresolved region.

`DD11_island_point.csv` is our independently replayed outward-interval
certificate for a closed dyadic rectangle containing L. Its run used 35,291
pair checks, maximum word depth 36, and a minimum leaf squared gap larger
than 2.84 × 10^-13. The exact rectangle endpoints are recorded as decimals
and hexadecimal binary64 values in the run manifest.

`DD11_island_zoom.csv` additionally certifies all 4,096 cells of a small
64-by-64 dyadic neighborhood around L, with cell side 2^-26. This entire
square lies inside the stable set. It was retained as certificate data,
not drawn as if it showed the boundary of the whole island.

## Caption suitable for the paper

Certified stable subsets of the twelve marked binary zipper families.
Rows give the plane orientations and columns the endpoint-order signatures.
Horizontal and vertical coordinates are Re(q) and Im(q), respectively; the
lower-half diagrams follow by complex conjugation. Blue cells certify every
parameter in a dyadic rectangle; the three OO disks and the real interval follow
from exact theorems. Gray parameters are unresolved, so blue-gray boundaries
are not asserted to be stable-set boundaries. Pink parameters belong to the
strict contraction lens but fail the necessary embeddedness disk bound;
white parameters lie outside that lens. The dot marks q*=3/8+i/4, used in
the curve gallery. The star L marks the Koch–wiggle island parameter cited
to Calegari; an independent interval rectangle certificate containing L is
included with the source.
