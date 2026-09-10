# Finite-field coloring obstruction screen, 2026-09-10

The graph in this note is **exactly** the graph with vertex set `F_q²` and
adjacency `(x-x')²+(y-y')²=1` in `F_q`. A finite-field realization is not a
Euclidean realization. A coloring of this graph can be pulled back through a
ring homomorphism, so it is useful for ruling out algebraic search families.

## A verified five-coloring of the entire F25 graph

Write `F25=F5[w]/(w²-2)` and write a vertex as
`(x,y)=(a+bw,c+dw)`, with `a,b,c,d` in `F5`. The formula

```
C(x,y) = a+c+d mod 5
```

is a proper five-coloring. This was found by enumerating the 156 rank-one
linear maps up to nonzero scalar; 12 have kernels containing no unit
direction. It was checked against the independently generated graph
`runs/modular-f25-k6/graph.json`: 625 vertices, 7,500 undirected edges, zero
monochromatic edges. The complete coloring vector is in
`runs/finite-quotients-f25/linear-5-coloring.json`.

Here is a short proof that does not depend on a SAT solver. By translation
invariance, suppose a unit direction has color zero. Then `a=-c-d`. The
coefficient of `w` in `x²+y²=1` gives

```
c(d-b)=bd.
```

If `b=0`, then `cd=0`. For `d=0` the real coefficient of the norm is `2c²`;
for `c=0` it is `3d²`. Neither can be 1 in `F5` because 2 and 3 are nonsquares.
If `b != 0`, put `t=d/b` and `u=c/b`. Then `u(t-1)=t`, excluding `t=1`.
The remaining four cases are:

| t | u | (a²+2b²+c²+2d²)/b² |
|---:|---:|---:|
| 0 | 0 | 2 |
| 2 | 2 | 0 |
| 3 | 4 | 0 |
| 4 | 3 | 2 |

Again the norm cannot be 1. Thus every unit direction has nonzero color,
which proves the coloring formula for every vertex pair.

The exhaustive additive-character spectrum is recorded as integer
coefficient lists in `runs/finite-quotients-f25/spectrum.json`. For a character
indexed by `h in F5^4`, the eigenvalue is `sum_{s in S} zeta_5^(h.s)`.
The list gives all 625 characters with multiplicities. Its least eigenvalue
is `-1-3sqrt(5)`. In particular it is less than -5. The floating-point values
in that file are labelled as approximate display values; the integer counts
determine the exact cyclotomic expressions. Comparing the 16 expressions
exactly and applying Hoffman's bound gives
`chi >= 1+24/(1+3sqrt(5)) > 4`; together with the formula above this establishes
`chi(F25²)=5`. The five-color upper bound is the part required for the search
obstruction and has the elementary proof above.

## Consequence for the current multiquadratic search

Let

```
R = Z_(5)[sqrt(2),sqrt(3),sqrt(5),sqrt(11)],
Z_(5) = {m/n in Q : 5 does not divide n}.
```

There is a ring homomorphism to `F25` given by

```
sqrt(2) -> w, sqrt(3) -> 3w, sqrt(5) -> 0, sqrt(11) -> 1.
```

Every image squares to the claimed radicand; all rational denominators in
`Z_(5)` are invertible. Composing this map in each coordinate with `C` gives
a proper five-coloring of **every** unit-distance graph whose coordinates
belong to `R`. Choosing `2w` instead of `3w` for `sqrt(3)` is also valid.
This is a proof about the whole ring, rather than a heuristic about tested
finite graphs. The independent verifier is checking the explicit lifted
models for the existing Euclidean artifacts.

This does **not** establish the same upper bound for the entire number field
`Q(sqrt(2),sqrt(3),sqrt(5),sqrt(11))`: denominators divisible by 5 prevent this
ring reduction. Moreover, `-1` is a square in `F25`, so the anisotropy
hypothesis required to extend a residue coloring to arbitrary field
coordinates does not hold. The ring pullback and the extra hypothesis for
fields are precisely distinguished in [Madore, paragraphs 2.2 and 3.1,
Proposition 3.2](https://arxiv.org/html/1509.07023).

Consequently, adding more exact circle intersections at squared separation
`4/3` or `8/3`, Minkowski sums, or rotations with denominators prime to 5
inside this ring cannot produce a non-six-colorable graph. An actionable
way to evade this particular obstruction within the existing exact MQ
implementation is to introduce a genuine denominator 5, such as the unit
rotation `(3+4i)/5`, and then construct actual unit links between its layers.
Escaping an upper-bound proof is only a necessary change of search family;
it supplies no lower bound on the chromatic number of a new candidate.

## F49 screen

For `F49=F7[w]/(w²-3)`, the exact unit-direction set has 48 elements and the
full graph has 2,401 vertices and 57,624 undirected edges. The helper
enumerated every rank-one and rank-two `F7`-linear quotient map, once per
row space using reduced row echelon form:

* No rank-one map has a loopless quotient.
* Of 2,850 rank-two maps, 1,118 create loops and 1,732 do not.
* The loopless maps give 848 distinct direction sets on 49 quotient vertices.
* Six-color SAT attempts returned `UNSAT_UNVERIFIED` for all 848 quotients.
  No proof certificates were requested for these screening solves.

These outcomes concern only the enumerated quotient coloring ansatz. They
do not prove that the full F49 graph is not six-colorable. Exact character
counts and individually labelled solve outcomes are retained in
`runs/finite-quotients-f49/`; `tools/finitecolor_helper.py` regenerates them.

## Primary literature and a bound that must not be used

[Vinh's 2005 paper](https://arxiv.org/html/math/0510092) defines precisely
this graph, gives degree `q-(-1)^((q-1)/2)`, and constructs an upper bound
`(p^n+p^(n-1))/2` for `q=p^n`. That gives 15 colors for q=25 and 28 for q=49,
so it does not settle this screen. Its Lemma 4 states the nontrivial
eigenvalue bound `|lambda|<=sqrt(q)`. The exact F25 eigenvalue above
disproves that bound as stated. The displayed lower bound in
[Vinh's 2006 Theorem 1](https://arxiv.org/html/math/0606482) uses the same
bound and would incorrectly force at least six colors for q=25.

The eigenvalues' Kloosterman-sum description comes from the original paper
[Medrano, Myers, Stark and Terras, *Finite analogues of Euclidean space*,
1996](https://www.sciencedirect.com/science/article/pii/0377042795002618).
The later [Vinh, *Explicit Ramsey graphs and Erdős distance problems*,
2008](https://emis.de/ft/5652) uses the `2sqrt(q)` scale for nontrivial
eigenvalues. No primary source found in this search supplies a six-color
certificate or certified non-six-colorability for the full F49 graph.
The explicit F25 formula and its elementary proof above resolve the useful
part of the current modular screen directly.

## Can the five-color formula extend to the whole number field?

No such extension has been proved here, and the primary-source search did
not find a theorem giving a five-color bound for the full multiquadratic
field or for the plane over all compass-constructible real numbers. There
are two rigorous obstructions to the suggested simple extensions.

First, no additive map from the whole additive group `K²` to `F5` can be a
proper coloring: every `v` is `5*(v/5)`, so every such map is zero. Even an
additive map on the connected unit-distance component cannot extend the
formula after admitting the rational unit directions `(3/5,±4/5)`. Their
sum, minus the unit vector `(1,0)`, is `(1/5,0)`. Additivity then forces
the color of the unit vector `(1,0)=5*(1/5,0)` to be zero. This argument
does not exclude nonlinear five-colorings.

Second, a finite coloring of the entire completion at 5 cannot be defined
by a Borel function of its 5-adic digits, valuation, and leading residues.
The completion contains `Q5`, where `-1` is a square by Hensel lifting
the simple root 2 of `X²+1` modulo 5. Hence `x²+y²` is isotropic over the
completion. [Bardestani and Mallahi-Karai, Theorem 1.1, *On a generalization
of the Hadwiger–Nelson problem*](https://arxiv.org/html/1507.05300) proves
that a nondegenerate quadratic graph over a characteristic-zero local
field has finite **Borel** chromatic number exactly when its form is
anisotropic. Thus a finite Borel extension on that whole completion is
impossible. This result does not exclude a five-coloring on the countable
number field: it is essential not to replace the Borel or whole-completion
conditions by a claim about arbitrary colorings of `K²`.

For completeness, an exact prime screen finds that the first odd prime
`p=3 mod4` for which all of `2,3,5,11` are squares or zero modulo p is
239; the next three below 1,000 are 359, 431 and 479. In particular, the
small-residue anisotropic reduction that proved the earlier
`Q(sqrt3,sqrt5,sqrt11)` five-color bound at 11 is unavailable for the
expanded field. These are simple residue calculations, not a theorem
excluding every possible valuation-based method.

[Moorhouse's preliminary draft, revised 3 March 2010, Theorem 1.2](https://ericmoorhouse.org/pub/chromatic.pdf)
shows that some number field attains the chromatic number of the real
plane, with a tower of quadratic extensions over an odd-degree field
whose plane is two-colorable. Its theorem does not assert that the
odd-degree starting field is `Q`, and supplies no five-color upper bound
for multiquadratic fields. It therefore does not rule out the denominator-5
search. The draft's finite-field table also gives `chi(F9²)=3`, independently
demonstrating why the incompatible bound in Vinh's 2006 paper must not be
used without checking it.
