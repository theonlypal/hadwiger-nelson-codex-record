# A separate multiquadratic coordinate family

This implementation uses exact real coordinates in

\[
K=\mathbb Q(\sqrt2,\sqrt3,\sqrt5,\sqrt{11}).
\]

Index `m` in a 16-entry coefficient array means
`product(sqrt(r[i]) for each set bit i of m)`, where `r=[2,3,5,11]`.
All square roots are positive. Multiplication of basis entries `a,b` gives
the entry `a XOR b`, multiplied by the product of radicands whose bits occur
in `a AND b`. Distinct products of these four primes represent independent
square classes, so these 16 elements form a rational basis of the field.
Consequently canonical rational arrays establish point equality exactly.

An exported vertex has separate arrays `x,y`. Each retained edge satisfies
`(x1-x2)^2+(y1-y2)^2=1` by rational coefficient equality. Vertex selection
inside a circle uses exact sign evaluation: rational intervals enclosing
each positive square root are refined until the sign is determined. No
floating-point comparison decides an algebraic equality or clipping bound.

## Primary-source construction

The starting construction is from Marijn J. H. Heule,
[Computing Small Unit-Distance Graphs with Chromatic Number 5,
arXiv:1805.12181v1](https://arxiv.org/html/1805.12181v1), Sections 2.1 and 4.1.
Using complex multiplication on real coordinate pairs, set

\[
h=(1/2,\sqrt3/2),\qquad q=(\sqrt{33}/6,\sqrt3/6).
\]

The 30 directions are `h^j*q^k` for `0<=j<6` and `-2<=k<=2`, and adjoining
the origin yields V31. Form V151 by keeping the points of V31+V31 whose
squared radius is at most one, then form V1939=V31+V151. The exact computed
cardinalities are 31, 151, and 1,939. The unit rotation

\[
\theta_4=(7/8,\sqrt{15}/8)
\]

produces the baseline union `V1939 union theta4*V1939`, which has 3,877
vertices and 26,814 retained exact unit edges. Its construction is a
baseline, not a claim of non-6-colorability.

For a smaller geometry check, the diamond
`{0,(1,0),h,(1,0)+h}` united with its copy rotated by
`theta3=(5/6,sqrt(11)/6)` gives the 7-vertex Moser spindle.

## Additional rotations

The explicitly checked norm-one multipliers

\[
r_1=(\sqrt2/2,\sqrt2/2),\qquad r_2=(1/3,2\sqrt2/3)
\]

introduce the independent square root of two. The code supports unions of
rotated copies and sums with a rotated six-direction unit star. These are
search families; entering a larger coordinate field does not establish an
increase in chromatic number.

Numerical squared-distance screening is permitted only to propose edges.
Every retained edge is checked exactly, and graph files explicitly decline
to claim induced-graph completeness. A verified non-6-colorable subgraph
would satisfy the target even if candidate screening omitted other edges.

Implementation: `hn/multiquadratic.py`, `hn/mq_construct.py`.

## Escaping the proved reduction at five

The separate `finite_field_screen.md` establishes a five-coloring for the
entire coordinate ring with rational denominators prime to five. Adding
sqrt(2) alone does not escape that bound. The denominator-five experiment
instead takes actual pairs of old centers at squared separations `10/3`,
`20/9`, or `25/9`. Their common unit neighbors are

\[
\frac{a+b}{2}\ \pm\ \frac{i(b-a)}2
\sqrt{\frac{4-s}{s}}.
\]

The respective square-root factors are `sqrt(5)/5`, `2*sqrt(5)/5`, and
`sqrt(11)/5`, all in the implemented field. A full screened pass through the
14,389-vertex first intersection closure found 61,056 new exact centers,
including 58,896 with a denominator divisible by five after canonical
normalization. Candidate enumeration uses a numerical broad phase; all
accepted separations and common-neighbor unit distances are exact.

The finite budget selected 4,268 complete sixfold rotational orbits. The
resulting graph in `runs/mq-pole5-v1/graph.json` has 39,997 vertices and
220,230 exact unit edges. Its six-core retains 17,617 vertices and 132,666
edges, including 3,324 newly introduced pole vertices. An independently
checked native SAT model gives a proper six-coloring of the full graph, so
this graph does not satisfy the requested target.

For 1,800 of those new six-core vertices, the recorded generating separation
is `20/9`. There is a stronger argument for these vertices than inspecting
their canonical rational denominators. Let
`K0=Q(sqrt(2),sqrt(3),sqrt(11))`. This field is unramified at five and its old
coordinates are integral at every prime above five. For `s=10/3` or `20/9`,
the normalized valuation of `s` is one. The two coordinate differences
cannot both have positive valuation, or their squared norm would have
valuation at least two. At least one difference therefore has valuation
zero. The perpendicular offset multiplier has valuation `-1/2`; one new
coordinate has valuation `-1/2`, which an integral midpoint cannot cancel.
These vertices have a genuine pole at every prime above five, so choosing a
different residue embedding does not restore the previous ring coloring.

The exact recorded-separation counts and the scope of that argument are in
`runs/mq-pole5-v1/pole5-valuation.json`. A separate native five-color
diagnostic on the exact five-core completed SAT in 276.733 seconds. Its
model satisfied every unrestricted CNF clause. Reversing the verified
elimination order extended that model to a proper five-coloring of the full
39,997-vertex graph; geometry and the resulting coloring were independently
checked. The complete witness is
`runs/mq-pole5-five-core-v1/full-coloring.json`. This particular candidate
does not improve the known lower bound, despite escaping the stated ring
coloring obstruction.

## Blocking one saved five-coloring

`hn/pole5_adaptive_extension.py` screened exact unit-neighbor centers against
the saved five-coloring of the pole graph. It selected 529 new centers, each
with neighbors of all five saved colors, at least six neighbors in the old
six-core, and both pole and integral neighbors. Thus the old coloring cannot
extend unchanged to any one of these new centers. This is a constraint on
that particular coloring, not a proof that every five-coloring fails.

The resulting graph has 40,526 vertices and 226,143 exact unit edges. All
529 new vertices survive the six-core, which grows to 18,879 vertices and
140,945 edges. Native six-color SAT completed in 9.910 seconds. A separate
five-core diagnostic then completed SAT in 372.628 seconds, and its coloring
was extended and independently checked on every full-graph edge. The full
five-color witness is
`runs/mq-pole5-adaptive-five-core-v1/full-coloring.json`. Both checks reject
this graph as a target certificate.

The old pole graph and this extension have 258 distinct directed unit
vectors. An exact scan found no pair differing by the angle
`theta4=(7/8,sqrt(15)/8)`. A translated, rotated, or reflected copy of the
published Heule realization needs such a pair because it contains the rays
`1` and `theta4` from its origin. This excludes isometric copies of that
specific realization. It does not exclude arbitrary embeddings of the same
abstract graph or all five-chromatic subgraphs.

## Completing the primary 510-vertex source realization

The coordinates and edge list come from the primary
[CNP-SAT repository](https://github.com/marijnheule/CNP-SAT), pinned to commit
`bb414955a6ef5f49f7df2b245b1e778aa67c068a`. Raw files are retained under
`sources/heule-cnp-sat/`. The arithmetic-only parser `hn/heule_import.py`
accepts rational operations and square roots supported by the declared
multiquadratic field; it does not execute source expressions.

Of the 510 source points, 357 already occur in the 40,526-vertex extension.
Adding the missing 153 gives 40,679 vertices and 227,005 exact unit edges.
This preserves the full primary 510-point graph, with 2,504 edges checked
against the pinned edge file and matched through an explicit index map.
The completion adds 294 old-new edges and 568 new-new edges. Only 65 of its
153 new vertices survive the six-core; the full six-core contains 18,944
vertices and 141,403 edges.

The native six-color test completed SAT in 42.620 seconds and checked all
2,012,894 canonical CNF clauses plus the full proper coloring. Results and
the exact source embedding are in `runs/mq-pole5-heule510-v1/`. The graph is
therefore another rejected non-six-colorability candidate.

A bounded local five-color repair allowed changes within distances zero
through three of the new vertices, reaching 29,481 free vertices at radius
three. Every restricted formula was UNSAT under the retained boundary
colors. These are restrictions on one boundary coloring and establish no
chromatic lower bound. The formulas, free-vertex maps, and outcomes are
retained in `runs/mq-pole5-heule510-local-five-v1/`. Five-colorability of this
completed graph has not been settled by those local diagnostics.

The standalone source geometry and an explicit five-coloring were checked
as a calibration. A fresh unrestricted four-color run on its 510 vertices
timed out after 60.136 seconds. No fresh non-four-colorability proof was
obtained, and the source's published chromatic claim is kept distinct from
the finite certificates produced here. None of these diagnostics satisfy
the requested non-six-colorability gate.
