# Literature and exact search constraints

Research date: 2026-09-10. This document is a research input, not a certificate of non-6-colorability. The requested success condition remains an exact finite Euclidean realization and an independently checked refutation of its 6-coloring formula.

## Current published bounds and artifacts

The primary sources inspected support **5 <= chi(R^2) <= 7**. A 6-chromatic finite unit-distance graph was not found in this literature search; the requested non-6-colorable graph would establish the stronger lower bound 7.

| Primary source | Relevant result / artifact |
| --- | --- |
| [de Grey, arXiv:1804.02385](https://arxiv.org/abs/1804.02385) | Original finite lower-bound-5 construction. |
| [Heule, arXiv:1805.12181](https://arxiv.org/html/1805.12181) | Minkowski sums, exact rotations, SAT refutations, and proof-guided graph reduction. Section 3.5 separates exact geometry from graph coloring. |
| [Heule's CNP-SAT repository](https://github.com/marijnheule/CNP-SAT) | Exact coordinate files `vtx/510.vtx`, `529.vtx`, `S199.vtx`, `T721.vtx`, corresponding edges, CNFs, proofs, and checks. Pin the downloaded commit and hash each input. |
| [Biere, Kauers, Ritirc: Checking distance one proofs](https://fmv.jku.at/dist1sqrtgb/) | Exact polynomial reduction for edge lengths; Singular and independently checked PAC proofs. |
| [Parts, arXiv:2010.12665v2](https://arxiv.org/html/2010.12665v2) | 509 vertices and 2442 edges; smallest unrestricted 5-chromatic example reported by the current source below. |
| [Haugland, arXiv:2608.04542v4](https://arxiv.org/html/2608.04542v4) | Latest inspected revision is August 17, 2026. A structured heptagonal construction yields a spindle-free 2131-vertex 5-chromatic graph. It acknowledges Heule's earlier 1441-vertex spindle-free example; 2131 is **not** the spindle-free size record. |

The ordinary Heule construction uses a 31-point direction set and Minkowski sums, then rotates copies by `theta_i = acos((2i-1)/(2i))`. In particular, `theta_4` has cosine `7/8` and sine `sqrt(15)/8`. Heule reports that `(S199 + S199) union theta_4(S199 + S199)` was relatively difficult to 5-color, but was still 5-colorable. This is a useful reproducibility baseline, not positive evidence for the target. [Heule, Sections 2.1, 4, 5](https://arxiv.org/html/1805.12181).

## Scope traps

* Seven colors for a **restricted map coloring** do not rule out arbitrary 6-colorings. The hypotheses in [Voronov, arXiv:2502.01958](https://arxiv.org/abs/2502.01958) concern map boundaries and finite local structure.
* A graph forbidding every distance in a nonzero interval around 1 is not a unit-distance graph. [Voronov, arXiv:2304.10163](https://arxiv.org/abs/2304.10163) concerns this interval variant.
* Unit-distance graphs drawn in the plane need not be planar graphs. A planar-graph coloring theorem is irrelevant to crossing edges in these realizations.
* A positive numerical residual, a SAT timeout, or failure of a coloring heuristic is not a noncolorability certificate.
* An extra incorrectly inferred edge can create false noncolorability. Missing true unit edges cannot invalidate a lower-bound proof for the remaining graph, although they weaken a search.
* A virtual edge enforcing distinct colors must be expanded into an actual unit-distance gadget with a certificate at **six** colors. A gadget established only for four colors cannot be imported into the target CNF.

## A provable obstruction to a large family of naive searches

[Madore, arXiv:1509.07023](https://arxiv.org/html/1509.07023), Proposition 4.6, proves that `Q(sqrt(3),sqrt(11))^2` is at most 5-colorable. Corollary 3.4 uses a residue field of size congruent to 3 modulo 4; Lemma 4.5 supplies an explicit 5-coloring of the unit-distance graph over `F_11^2`.

The following extension is derived here, rather than quoted as a theorem from the paper. It rules out the entire familiar multiquadratic coordinate family:

**Every unit-distance graph with coordinates in `Q(sqrt(3),sqrt(5),sqrt(11))^2` is 5-colorable.**

For the coordinate ring with rational denominators coprime to 11, the proof is elementary. Reduce coefficients modulo 11 and set

```
sqrt(3) -> 5,  sqrt(5) -> 4,  sqrt(11) -> 0.
```

These assignments satisfy each defining relation: `5^2 = 3`, `4^2 = 5`, and `0^2 = 11` modulo 11. Independence of the eight standard multiquadratic basis monomials makes this a well-defined ring map. An edge's exact equation maps to `(dx)^2 + (dy)^2 = 1` in `F_11`, so its endpoints cannot collapse. Compose this graph homomorphism with the table below.

For arbitrary rational denominators, use the same local-field argument as Madore: `sqrt(3)` and `sqrt(5)` lift from their simple roots modulo 11 to `Q_11`, by Hensel's lemma. Thus the coordinate field embeds into `Q_11(sqrt(11))`, whose residue field is `F_11`. Madore's Corollary 3.4 applies. This paragraph establishes existence of a coloring; the elementary ring map above is the direct finite implementation when denominator conditions hold.

Here is the 121-cell coloring from Madore's Lemma 4.5, transcribed and independently checked in this session using integer arithmetic on all **1452 directed edges / 726 undirected edges**:

```
31021234201
12104123432
21230212340
02313404134
40121340212
34012134121
12340310402
21434023034
42343402340
04121040232
40412303013
```

Rows and columns are indexed by `0,...,10` cyclically. Recheck by enumerating `dx,dy in 0,...,10` with `(dx*dx+dy*dy)%11 == 1`, then comparing every cell with its translate. The 12 directions are `(0,+/-1)`, `(+/-1,0)`, `(+/-3,+/-5)`, and `(+/-5,+/-3)` modulo 11.

Consequences: arbitrary repeated Minkowski sums, translations, and rotations whose coefficients remain in this field cannot yield a 6-chromatic graph, let alone the requested 7-chromatic graph. The `sqrt(15)/8` spindle rotation stays inside the field. Adjoining `sqrt(2)` defeats this particular residue-field obstruction because 2 is a nonsquare modulo 11. That observation does **not** prove a higher chromatic number or the absence of other coloring obstructions.

## Exact form of the heptagonal directions

Haugland's Appendix A gives paths covering the 740-vertex seed; Section 3 explains the unions producing the 2131-vertex graph. The paper identifies nearby points using a floating-point enumeration lemma. That is not a substitute for the user's exact certificate. Section 4 reports at most six candidate equivalence classes of 4-colorings of the generator lattice; it explicitly does not prove that those local colorings extend globally. [Haugland v4](https://arxiv.org/html/2608.04542v4).

The following algebraic reduction was derived and checked symbolically in this session. It allows the generator construction to be performed without tolerance. Represent planar points as complex numbers. Let

```
t = exp(pi*i/21),
Phi_42(t) = t^12 + t^11 - t^9 - t^8 + t^6 - t^4 - t^3 + t + 1 = 0,
7*v = 4*t^11 - 4*t^10 - 6*t^9 + t^8 + t^7 + 5*t^6
      + t^5 - 2*t^4 + 3*t^2 + t + 2.
```

The 84 directions are `u_(2j)=t^j`, `u_(2j+1)=v*t^j`, for `0 <= j < 42`. Complex conjugation is `t -> t^-1 = t^41`. Exact polynomial reduction verifies `v*conjugate(v)=1`; powers of `t` also have norm 1. The intended complex embedding is fixed by `t=exp(pi*i/21)`, and real coordinates are `(Re(z), Im(z))`.

One derivation starts with `s=exp(pi*i/42)` and `i=s^21`. The difference `Q_0-R_0` in Haugland's table equals

```
i*s^7/(s^12-s^-12) - i*s^35/(s^24-s^-24).
```

Reduce modulo `Phi_84(s)=s^24+s^22-s^18-s^16+s^12-s^8-s^6+s^2+1`; only even powers survive. Set `t=s^2` to obtain `v` above. Each point now needs just 12 rational coefficients, and equality is equality of reduced coefficient vectors. An edge is checked by reducing `(z-w)*conjugate(z-w)-1` to zero.

The final spindle rotation can be represented by `rho=(7+i*sqrt(15))/8`. Since `i*sqrt(3)` lies in `Q(t)`, the complex field only needs a further `sqrt(5)` extension. Equivalently, exact real coordinates fit inside `Q(zeta_84)^+(sqrt(5))`. The degree-12 representation above is for complex points; it does not claim that their separate real and imaginary parts both lie in `Q(t)`.

## Field screen for the heptagonal extension

No theorem asserting 5- or 6-colorability of `Q(zeta_84)^+(sqrt(5))^2` was found. This is a limit of the search, not proof that the field is suitable.

A concrete residue calculation shows why Madore's simple `F_11` obstruction does not transfer: on the real cyclotomic field `Q(zeta_84)^+`, the residue degree at 11 is the least positive `f` with `11^f = +/-1 (mod 84)`. The powers modulo 84 are `11,37,71,25,23,1`, so the degree is 6. The residue field has size `11^6`, congruent to 1 modulo 4, and the anisotropy hypothesis of Corollary 3.4 fails.

This computation excludes only that particular application of the theorem. It does not establish any chromatic lower bound. Likewise, modular noncolorability does not imply that a finite-field graph has a Euclidean embedding.

## Concrete next search

1. Reconstruct the Haugland seed from Appendix A in the exact degree-12 representation. Verify all retained edges by conjugate-norm reduction.
2. Add the `sqrt(5)` spindle extension and independently certify a 5-coloring of the published graph as a control. Check its 6-coloring formula explicitly, preserving the coloring witness.
3. For new searches, use exact unions of heptagonal and multiquadratic direction sets, or add a `sqrt(2)` rotation to the multiquadratic family. These choices escape the demonstrated `F_11` coloring map. First calculate cross-copy edges and finite modular coloring obstructions; do not spend larger SAT searches on a field already proved 5-colorable.
4. Search boundary-pair constraints at six colors, retain the exact full unit-distance gadget for each derived constraint, and issue an UNSAT proof for any final candidate. An UNSAT result for five colors would be new but still would not meet the requested certificate.

None of these proposed searches comes with a guarantee that the requested graph exists or will be found.
