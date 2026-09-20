# Finite separation certificates and parameter neighborhoods

`tools/verify_rational.py` is a small, separately runnable verifier for explicitly specified rational complex multipliers. It uses only Python's standard library. The browser renderer and floating-point search are independent from it.

## Scope

Input the two complex multipliers as exact decimal or fraction **strings** and each branch's orientation as 0 (direct) or 1 (reversing). The translations are -1 and +1. A successful version-2 result certifies the rational center and an explicit closed product of complex disks about its two multipliers. It does not certify connectedness, tile or boundary membership. Decimal values rounded from trigonometry describe a different rational parameter; a displayed slider point is covered only after its distance to the exact center has itself been bounded rigorously.

Example `quarter_turn_rational.json`:

```json
{"multipliers": [["0", "3/5"], ["0", "3/5"]], "orientation": [0, 0]}
```

From `4D/`:

```sh
python3 tools/verify_rational.py tools/quarter_turn_rational.json --output tools/quarter_turn_certificate.json
python3 tools/check_certificate.py tools/quarter_turn_certificate.json
python3 tests/rational_verifier_test.py
```

## Why rejection proves disconnection

For each multiplier, the verifier encloses its Euclidean modulus above by a rational number `q_i < 1`. The upper bound is obtained by an exact integer square root, followed by a rational square comparison and one upward grid step when necessary. With `q = max q_i`, the closed disk of radius `R = 1/(1-q)` centred at zero is invariant under both maps, since their translations have modulus one.

A word `u` maps that disk into a disk with exactly computed rational centre and radius bounded by `R product(q_i : i in u)`. The code composes maps outermost first and carries conjugation parity. At each common depth it retains pairs whose outermost symbols differ. A pair is discarded only if the exact squared distance between centres is strictly greater than the exact square of the sum of enclosing radii. Discarding a parent also discards all its descendants because their attractor cylinders are subsets of the parent enclosure.

If every pair is discarded, the first-level attractor pieces are disjoint and their union is disconnected. If any pair survives, or the finite work budget cannot cover the next frontier, the result is unresolved. All integer and rational operations in these comparisons are exact. This is a one-sided finite disconnection test; it has no positive connectedness branch.

## A complete finite cover

The successful JSON stores the discarded word pairs in `separation_cover.leaves`, with words composed outermost first. They form a complete prefix antichain in the four-branch tree of paired suffixes, starting with `-` and `+`. Leaves may occur at different depths. This finite cover proves separation for **every** infinite cross-prefix address pair.

`check_certificate.py` reads that JSON without importing the search or repeating its exploration. It verifies the prefix cover is complete, evaluates word centers using independent rational real 2-by-2 matrices, checks the invariant enclosing disk, and verifies every separation inequality. Removing a leaf, duplicating a leaf, changing coefficients without their checksum, shrinking the norm bounds, or inflating the gap is rejected. A saved hash alone is never sufficient evidence.

## An explicit disconnected neighborhood

For a separated leaf with center distance `d` and sum `s` of the enclosing radii, the search uses the positive rational lower bound

`delta_leaf = (d²-s²)/(sqrt_upper(d²)+s) <= d-s`.

The minimum `delta` over the complete leaf cover bounds the distance between the first-level attractor pieces below. On a fixed orientation chart with both systems' contraction ratios at most `qbar < 1`, their contact gaps differ by at most `2 eta/(1-qbar)²`, where `eta = max(|lambda_minus-mu_minus|, |lambda_plus-mu_plus|)`.

Set `qbar=(1+q)/2`, where `q` is the rational outward contraction bound at the center, and choose

`eta = min((1-q)/2, delta(1-qbar)²/4, min_s max(|Re lambda_s|,|Im lambda_s|)/2)`.

Every parameter in that closed product of complex disks is a nonzero strict contraction pair and has contact gap at least `delta/2`. The checker verifies these perturbation inequalities directly. The neighborhood has four real dimensions and fixed orientations; it is **not** an unverified angular grid cell.

## Reproduction and limitations

The output records canonical rational coefficients, orientation bits, their SHA-256, verifier version, norm bounds, containing radius, depth, work limit and pair counts. A successful result additionally records the complete separation cover, rational gap, minimum squared margin and parameter-neighborhood bounds. The separate checker validates the saved proof object. Rerunning the search on the `input` object also reproduces it.

The 80-bit rational norm grid can fail to resolve a contraction arbitrarily close to one. Such an input is rejected as out of the supported computational range. Work and depth limits produce unresolved results. This verifier has targeted mathematical review and regression tests; it is not a formal proof-assistant development. Version changes require rerunning the suite and reviewing the soundness argument.
