# A degree-sixteen real field from two Voronov rotations

The proposed field and rotations below are exact. Their characteristic-five
properties exclude the previously used F25 coloring map, but do not imply
non-five-colorability or non-six-colorability of any Euclidean graph.

## Primary construction and provenance

[Voronov, Neopryatnaya and Dergachev, v4, sections5–6](https://arxiv.org/html/2106.11824v4)
give the L10,2 generators

```
phi0=(sqrt6+sqrt2+i*(sqrt6-sqrt2))/4,
phi1=(sqrt6+i*sqrt3)/3.
M1={0} union {phi0^a phi1^b:0<=a<24, b in{-1,0,1}},
M2=clip(M1+M1,1),  M3=M2+M1.
```

Their Table6 includes rotations determined by `z=2,w=1+phi0^j` for
`j=1,5,7`. These are published five-chromatic constructions, not claimed
six-chromatic examples. The paper's concluding solver difficulty is not
a non-colorability certificate.

The author's cited [repository](https://github.com/vsvor/dist-graphs) was
pinned at `e3714d0a156f6ed151d4521debb2b89ce4f1075c`. Selected files are in
`sources/voronov-dist-graphs/`, with exact URLs, byte counts and SHA256 hashes
in `manifest.json`. They include the Sage notebook, Mathematica M1/M2/M3
coordinate lists, the rotation expressions, and full graph2/graph3 coordinate
ZIPs. No downloaded notebook was executed.

The source notebook's alternate generator is exactly
`a=phi0^4*conjugate(phi1)`, so it gives the same M1 orbit. A restricted
arithmetic parser independently read the author's lists and checked equality
with the exact reconstructed point sets:73,865 and32,257 vertices respectively.
The M3 comparison uses `artifacts/voronov_m32257_points.json`.
See `tools/voronov_source_audit.py` and `research/voronov_source_audit.json`.

For j=1,5,7, `|w|>1`; these points are not in clipped M2. They and2 are in
M3 because M2 contains M1, hence M3 contains M1+M1. This resolves the apparent
M2 wording issue without assuming membership in a clipped set incorrectly.

## Exact rotation formulas

Put `K=Q(sqrt2,sqrt3)`. For `w=1+phi0^j`, write `w=a+i*b` and `s=a²+b²`.
Set

```
h=(s-1)*(9-s),
psi=((s+3)+i*sqrt(h))/(4*w),
Re(psi)=((s+3)*a+sqrt(h)*b)/(4*s),
Im(psi)=(sqrt(h)*a-(s+3)*b)/(4*s).
```

The square root is the positive real one. The identities

```
(s+3)²+h=16*s,
(5-s)²+h=16
```

prove both `|psi|=1` and `|2-psi*w|=1` exactly. Thus each selected rotation
comes with an actual unit edge joining the base and rotated copies.

For the two selected cases:

| j | a | b | s | h |
|---|---|---|---|---|
|1|`(4+sqrt2+sqrt6)/4`|`(sqrt6-sqrt2)/4`|`(4+sqrt2+sqrt6)/2`|`5-sqrt3+3sqrt6+3sqrt2`|
|7|`(4+sqrt2-sqrt6)/4`|`(sqrt2+sqrt6)/4`|`(4+sqrt2-sqrt6)/2`|`5+sqrt3+3sqrt2-3sqrt6`|

Call these radicands h1,h7 and their positive square roots r1,r7. The j=5
radicand is `5+sqrt3+3sqrt6-3sqrt2`; it is not needed for this proposal.

Every selected s is an algebraic unit:

```
s^4-8*s³+20*s²-16*s+1=0,
1/s=-s³+8*s²-20*s+16.
```

Consequently the displayed divisions introduce no denominator divisible
by5. The candidate coordinate ring lies in
`Z[1/6,sqrt2,sqrt3,r1,r7]`.

The rotations were also compared directly against the literal expressions
in the pinned author's rotation file. A restricted parser and exact
cross-multiplication give

```
source rotation row2 = phi0^13 * psi1,
source rotation row3 = phi0^19 * psi7.
```

The source's complicated nested radicand in row3 is exactly `h7/2`.
Since M3 is invariant under phi0, these formula differences produce the
same two rotated M3 point sets.

`tools/voronov_algebra_audit.py` verifies all these identities using
constant-plus-radical arithmetic over the existing exact multiquadratic
implementation. It imports neither the new field nor its verifier. Expected
sixteen-coordinate real and imaginary parts of psi1,psi7 are saved in
`research/voronov_rotation_algebra.json` for independent implementation checks.

## Degree and positive embeddings

The useful simplification is

```
h1*h7=-14+48*sqrt2.
```

The signs under the four real embeddings of K are:

| sign sqrt2 | sign sqrt3 | h1 | h7 | h1*h7 |
|---|---|---|---|---|
|+|+|+|+|+|
|−|+|−|+|−|
|+|−|+|+|+|
|−|−|+|−|−|

A square in a totally real field cannot be negative at any real embedding.
Thus h1,h7 and their product are all nonsquares in K. Their square classes
are independent, and

```
E=K(r1,r7),    [E:Q]=16.
```

The selected embedding is real, with both new radicals positive. For example,
the rational bounds `7/5<sqrt2<3/2`, `17/10<sqrt3<7/4`, and
`12/5<sqrt6<5/2` prove positivity of h1,h7 and negativity at the specified
other embeddings without numerical tolerances. Each h_j has polynomial
`X^4-20X³+1432X-4412`; each new square root has degree8 over Q, with
polynomial `X^8-20X^6+1432X²-4412`. The radicals should therefore not be
described as individual degree-four numbers.

For complex arithmetic, adjoining i doubles the degree to32. Independence
also follows directly over K(i): a member a of K can become square there
only if a or−a is square in K, and the mixed-sign radicands and their
product exclude both possibilities. The proposed real basis is

```
(sqrt2)^e2 * (sqrt3)^e3 * r1^e1 * r7^e7,  each exponent in{0,1}.
index=4*(e1+2*e7)+(e2+2*e3).
```

## Every F25 image fails jointly

Write `F25=F5[w]/(w²-2)`. The four maps from the base coordinate ring have
the following images. A tuple in the audit JSON represents coefficients
of `1,w`.

| sqrt2 | sqrt3 | h1 | h7 |
|---|---|---|---|
|w|2w|2+w, nonsquare|3, square|
|w|−2w|3, square|2+w, nonsquare|
|−w|2w|3, square|2−w, nonsquare|
|−w|−2w|2−w, nonsquare|3, square|

This is exhaustive: the roots of2 and3 are respectively ±w and ±2w.
The program enumerates all25 possible square roots for each image, rather
than inferring square status from floating-point or symbolic heuristics.
Equivalently, the product has image `1±3w`, whose norm to F5 is
`1-2*3²=3`, a nonsquare. A product of two squares cannot have this norm.
Thus no map of the base ring to F25 extends after both roots are adjoined.
The optional choices `i=±2` do not alter this real-coordinate obstruction.

Each individual radical allows two of the four base maps, illustrating why
one rotation alone does not remove every old F25 coloring route. The joint
failure removes that particular upper-bound method only. A different
quotient, a nonlinear rule, or a direct finite coloring may still use five
or six colors.

## Coupling geometry and bounded proposal

The initial proposal is

```
M3 union psi1*M3 union psi7*M3.
```

The three copies intersect only at0. Its vertex count is therefore exactly
`3*32257-2=96769`, not merely bounded by that number. Each rotated copy
has its guaranteed unit link to the base point2, in addition to any other
exactly verified cross-copy edges.

There are no unit edges directly joining nonzero vertices of the two rotated
copies. To see this, write
`psi_j=a_j*(1+i*r_j/(s_j+3))`, with `a_j=(s_j+3)/(4*w_j)` in K(i).
For base points v,u, put `P=a1*v*conjugate(a7*u)`. The coefficients of r1
and r7 in their squared distance force `Im(P)=0`; the coefficient of r1*r7
forces `Re(P)=0`. Independence of the four extension basis terms then gives
`P=0`, hence v=0 or u=0. The two branches are coupled through the entire
shared base graph; that structural observation does not imply a coloring.

If further coupling is needed after testing this graph, adding
`psi1*psi7*M3` creates a fourth copy with guaranteed unit links to both
rotated branches, by rotating the already certified base-branch links.
This is an implementable follow-up within the same field, not a claim that
its chromatic number increases. Any retained edge must still be independently
checked, and any non-six-colorability claim still requires an independently
verified refutation of the unrestricted six-color CNF.
