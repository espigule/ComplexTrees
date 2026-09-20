# Exact finite exclusion for a normalized planar pair

This directory contains a small research computation and an independent verifier for

\[
F(z)=-1+aJ_p(z),\qquad G(z)=1+bJ_q(z),\qquad
J_0(z)=z,\quad J_1(z)=\overline z,
\]

where the real and imaginary parts of `a` and `b` are rational and both multipliers have modulus less than one. The three orientation classes are represented by `(p,q)=(0,0),(0,1),(1,1)`. The implementation also accepts `(1,0)`.

The two possible outputs are **DISCONNECTED** and **UNRESOLVED**. The latter means that a finite depth or pair budget was exhausted. It never means connected. These benchmark records exercise the finite exclusion criterion. The ambient atlas, exact families and zipper analysis are developed in the accompanying article.

## Reproduce the supplied calculations

Python 3.10 or later is sufficient. There are no external dependencies. From this directory:

```bash
python3 -m unittest -v test_exact_exclusion.py
python3 exact_exclusion.py benchmarks --output-dir certificates --max-depth 8 --max-pairs 20000
python3 verify_certificate.py certificates/quarter_pp.json certificates/imaginary_two_thirds_pp.json certificates/imaginary_two_thirds_pm.json certificates/imaginary_two_thirds_mm.json certificates/interval_half_pp.json certificates/interval_half_pm.json certificates/interval_half_mm.json
```

To inspect one named example or give different rational inputs:

```bash
python3 exact_exclusion.py list
python3 exact_exclusion.py benchmark imaginary_two_thirds_pm --output example.json
python3 exact_exclusion.py run --a-real 1/3 --a-imag 2/5 --b-real=-2/7 --b-imag 1/4 --p 0 --q 1 --output custom.json
python3 verify_certificate.py custom.json
```

The verifier is self-contained and imports no code from the producer. It can be copied and used with a JSON certificate on its own. The producer composes affine coefficients forward through words; the verifier reconstructs each cylinder centre independently by applying the individual letters from right to left to zero.

## Why exhaustion is a proof

The program obtains a rational `q_bound` with

\[
0\leqq_{\mathrm{bound}}<1,\qquad |a|^2\leqq_{\mathrm{bound}}^2,\qquad |b|^2\leqq_{\mathrm{bound}}^2.
\]

It uses an exact rational square root when one exists. Otherwise an integer square-root calculation followed by an exact ceiling test gives a dyadic upper bound. No floating-point square roots or numerical tolerances are used.

Let `R=1/(1-q_bound)`. Both maps preserve the closed disk centred at zero with radius `R`, because `1+q_bound*R=R`. Thus the attractor lies in that disk. A word `u` of length `n` has the form

\[
T_u(z)=B_u+A_uJ_{r_u}(z),\qquad |A_u|\leqq_{\mathrm{bound}}^n,
\]

so its cylinder is contained in the closed disk of centre `B_u` and radius `q_bound^n*R`.

Start with the pair `(F,G)`. A pair of equally long words `(u,v)` is discarded only if the exact rational inequality

\[
|B_u-B_v|^2>\bigl(2q_{\mathrm{bound}}^nR\bigr)^2
\]

holds. Otherwise replace it with all four pairs `(uF,vF)`, `(uF,vG)`, `(uG,vF)`, `(uG,vG)`, unless a stated limit has been reached. The children cover the parent cylinder product. If a finite tree has only discarded leaves, every possible intersection of `FK` and `GK` has been excluded. Those two nonempty compact pieces are disjoint, so their union `K` is disconnected.

For composition, if `T(z)=B+A J_r(z)` and `S(z)=D+C J_s(z)`, the producer uses

\[
(T\circ S)(z)=B+A J_r(D)+A J_r(C)J_{r\mathbin\oplus s}(z).
\]

The verifier checks the contraction and invariant-disk bounds, every excluded leaf inequality, the recorded rational centres and gaps, all four children of each internal node, absence of missing or duplicate leaves, frontier reasons, computation limits and tree statistics. A claimed disconnected result with a surviving frontier is rejected.

## Computed benchmark results

These values were actually produced and independently verified with depth limit 8 and pair-test budget 20,000.

| Multipliers | Orientation `(p,q)` | Result from this algorithm | Pair tests | Excluded leaves | Frontier leaves | Maximum depth |
|---|---:|---|---:|---:|---:|---:|
| `a=b=1/4` | `(0,0)` | DISCONNECTED | 1 | 1 | 0 | 1 |
| `a=b=2i/3` | `(0,0)` | DISCONNECTED | 69 | 52 | 0 | 6 |
| `a=b=2i/3` | `(0,1)` | DISCONNECTED | 69 | 52 | 0 | 6 |
| `a=b=2i/3` | `(1,1)` | DISCONNECTED | 69 | 52 | 0 | 6 |
| `a=b=1/2` | `(0,0)` | UNRESOLVED | 29 | 21 | 1 | 8 |
| `a=b=1/2` | `(0,1)` | UNRESOLVED | 29 | 21 | 1 | 8 |
| `a=b=1/2` | `(1,1)` | UNRESOLVED | 29 | 21 | 1 | 8 |

For `a=b=2i/3`, the verified common contraction bound is `q_bound=2/3` and `R=3`. This bound is distinct from the established atlas coordinate `rho=2/(|a|+|b|)`, which equals `3/2` for this example. All three runs split 17 pairs and exclude 52 terminal pairs. Their equal statistics are specific to these symmetric examples; they do not assert that the orientation classes agree in general.

For `a=b=1/2`, the attractor is independently known to be `[-2,2]`: conjugation fixes the real line, and the two images of that interval are `[-2,0]` and `[0,2]`. The algorithm deliberately leaves the shared-endpoint pair unresolved. This external exact argument supplies connectedness; finite survival supplies none.

The six test methods include 248 rational word-composition checks across all four parity combinations, each evaluated at three rational complex points, as well as independently reconstructed centres. They also test exact contraction bounds, known interval survival, zero multipliers, depth and budget exhaustion, and rejection of deliberately corrupted certificates.

The JSON files are deterministic and contain no timestamps. Rerunning the producer with the same inputs and limits reproduces them byte for byte.
