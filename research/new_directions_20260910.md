# New construction directions examined on 2026-09-10

These are inputs to the non-6-colorability search, not a claimed solution.

## The June 2026 independence-ratio paper

[Ducz and Varga, arXiv:2606.28157v1](https://arxiv.org/html/2606.28157) prove existence of a finite unit-distance graph with independence ratio below 1/4. The explicit seed is a 29-point configuration whose **geometric fractional** chromatic number exceeds 4.0007. Subsequent amenability and additive blow-ups establish the independence-ratio result; the paper does not enumerate the enormous resulting graph. Its two added points are ordinary degree-one vertices, so they do not directly raise the integer chromatic number of the seed.

Write `w=(1+i*sqrt(3))/2`, `z=(5+i*sqrt(11))/6`, and `r=sqrt((415+79*sqrt(33))/8)`, choosing positive real roots. The new points are

```
p = 3 + 17*w/8 - 7*z/8 + 2*w*z
    + sqrt(5)*(-1/4 + w/8 - z/8 + w*z/4),
q = 11/4 + 13*w/8 - z/8 + 2*w*z + i*r*(-w+z+w*z)/8.
```

The transferable search idea is to intersect circles centered at existing vertices, using radii already occurring as distances in the configuration, and screen augmentations against constrained colorings. [Sections 4 and 6](https://arxiv.org/html/2606.28157).

## Acquired exact artifacts

The authors' [supplement](https://users.renyi.hu/~akos/ep1070/) contains symbolic coordinates, congruences, a rational dual witness, and a verification script. The linked archive was downloaded to `sources/ducz-varga-2606.28157/snail.zip`:

```
bytes: 45637428
sha256: ff6c9fc9df606ee15be9b4e3d14e7144e3d829fbd7fbeca3bc21bdd50e43667c
```

Selected files were extracted with individual hashes in `manifest.json`. The NumPy symbolic coordinate payload was inspected with `pickletools` and decoded with a restricted allowlist of NumPy and SymPy constructors. A text-only copy is `vertices-exact-expressions.json`. Its first two points are `p,q`; its remaining 27 points are the seed. This differs from the paper's vertex numbering.

Fresh symbolic verification reduced the squared lengths of all **51 supplied edges** to exactly 1. The two added points both connect only to file index 4 in the supplied adjacency. Results are in `geometry-inspection.json`. No induced-graph completeness claim is made, and the geometric fractional dual certificate was **not** rerun. The script's published rational lower bound is `4000716307/1000000018`.

## Two deductions that constrain how to use this paper

First, geometric/fractional lower bounds cannot themselves reach this task's target. The paper recalls the Croft construction, giving `chi_f(R^2) < 4.36`. Therefore any lower bound using only `n/alpha(G)` or fractional chromatic number remains below 4.36 and can certify integer chromatic number at most 5. A non-6 certificate must use stronger coloring constraints. This conclusion follows from the bound and is not a conjecture about integer chromatic number. [Section 2.1](https://arxiv.org/html/2606.28157).

Second, the complicated new radical in `q` still permits the previously derived 5-coloring obstruction. Under the coordinate ring reduction modulo 11,

```
sqrt(3) -> +/-5, sqrt(5) -> +/-4, sqrt(11) -> 0,
r^2 = (415+79*sqrt(33))/8 -> 1,
r -> +/-1.
```

All displayed denominators are invertible modulo 11. Consequently this extended ring still maps to `F_11`, and all exact unit-distance edges map to the 5-colorable graph on `F_11^2`. For the entire number field, simple-root lifting gives the same residue-field extension of the local argument in [Madore, Corollary 3.4](https://arxiv.org/html/1509.07023). Thus merely using `p,q` to enlarge a coordinate field cannot reach six colors. This calculation is derived here and was not stated in the June paper.

## An exact circle-intersection operation that escapes that obstruction

Let `a,b` be distinct existing points, `d=b-a`, and let their exact squared separation be the rational number `s`, with `0<s<4`. Their two common unit neighbors are

```
c_plus  = (a+b)/2 + (i*d/2)*sqrt((4-s)/s),
c_minus = (a+b)/2 - (i*d/2)*sqrt((4-s)/s).
```

Indeed the two summands in `c_plus-a = d/2 + i*d*sqrt((4-s)/s)/2` are perpendicular and have squared lengths `s/4` and `(4-s)/4`. Their sum is 1; the same identity holds for the other three center-neighbor pairs. An implementation should nevertheless check all four norm identities independently for each accepted pair.

The new square root fails to exist in `F_11` when `((4-s)/s) mod 11` is nonsquare. For nonzero `s mod 11`, the excluded residues are `s in {5,6,9,10}`. This is a field-choice filter, not evidence of high chromatic number.

Two particularly cheap choices introduce only the already supported generator `sqrt(2)`:

| Exact separation squared | Radical | Common-neighbor offset |
| --- | --- | --- |
| `4/3` | `sqrt(2)` | `+/- i*(b-a)*sqrt(2)/2` |
| `8/3` | `sqrt(2)/2` | `+/- i*(b-a)*sqrt(2)/4` |

Both choices have nonsquare reduction modulo 11. Adding **both** intersection points couples the old field to the new field through explicit unit edges to both centers. New points of degree only two cannot by themselves obstruct six colors, so the useful next step is to close further or select intersection points with additional neighbors and compare the resulting 6-core. This avoids treating a collection of pendant additions as chromatic progress.

## Verified primary sources on dense constructions

[Alon, Bloom, Gowers, Litt, Sawin, Shankar, Tsimerman, Wang, and Wood, arXiv:2605.20695](https://arxiv.org/html/2605.20695) give a mathematical proof of the 2026 disproof of Erdos's `n^(1+o(1))` unit-distance bound. Their construction uses CM fields of increasing degree, bounded root discriminant, controlled small prime splitting, many norm-one elements, and a bounded window in the full Minkowski embedding. Fixing one complex embedding turns the finite window into planar points. The increase in field degree is essential to the asymptotic method; a small fixed-degree picture is not its asymptotic counterexample.

[Sawin, arXiv:2605.20579](https://arxiv.org/html/2605.20579) makes an exponent explicit: at least `n^1.014114/C` ordered unit-distance pairs for arbitrarily large `n`. His Lemma 2 gives a usable finite template. If an injectively projected lattice `Lambda` has `M` short vectors projecting to unit length, a translated radius-`R` window has average directed unit degree at least `(1-1/R)^(2d)*M` for a suitable translation. Later lemmas choose ideals and relative norms to generate those vectors.

Neither paper proves any corresponding chromatic-number lower bound. High average degree alone does not force high chromatic number; even bipartite graphs can have arbitrarily large average degree. For this project, dense arithmetic windows are a possible candidate source only after exact field construction and a modular coloring screen. The directly implementable step with the current exact arithmetic is the rational-separation intersection closure above, not a claim to have instantiated a class-field tower.

## Next computation

Enumerate pairs in the exact V1939/Heule seed whose squared separation is `4/3` or `8/3`, add their two exact common unit neighbors, then extract exact unit edges and the 6-core. Preserve the pair provenance and radical used for every new point. SAT-check six colors with a bounded timer, retaining full models or UNKNOWN; test five only as a diagnostic. Repeat a closure round only if the first round produces substantial new edges surviving the 6-core. The geometry and SAT certificates remain separate.

## First execution of the intersection search

Implemented in `hn/mq_intersections.py`. The first complete rounds, with no vertex-cap truncation, produced:

| Seed | Accepted center pairs | New points | Total points | Exact edges | New points surviving 6-core |
| --- | ---: | ---: | ---: | ---: | ---: |
| V1939 | 6,225 | 12,450 | 14,389 | 116,742 | 12,354 |
| Heule union | 12,450 | 24,900 | 28,777 | 233,550 | 24,708 |

All accepted center pairs had separation squared `4/3`; none had `8/3`. Independent geometry checks passed with `tools/verify_mq.py`. Both six-color and five-color Glucose4 runs reached their 60-second caps for both graphs; all four outcomes are **UNKNOWN_TIMEOUT**, not lower-bound evidence. Runs, exact coordinates, center-pair provenance, CNFs, and verification reports are in `runs/mq-intersections-v1939-r1` and `runs/mq-intersections-heule-r1`.

Exact point-set comparison and embedded edge comparison showed that the Heule closure consists of the V1939 closure and its `theta4` image, sharing one point, with 66 cross edges. This structural check is recorded in `layer-analysis.json`.

A separate exact distance screen for rotations introducing `sqrt(2)` found 42 Heule-seed pairs at squared separation `3/4`, and no V1939 pairs there. For a pair `(a,b)` at this separation, rotating a copy around `a` by `rho=1/3+2*i*sqrt(2)/3` guarantees a unit edge from `b` to `a+rho*(b-a)`. None of the tested seed pairs had squared separation `9/4`, `19/4`, `33/4`, or `51/4`. This suggests using actual endpoints as rotation centers instead of the arbitrary origin. The screen is saved in `runs/mq-rotation-distance-screen.json`.
