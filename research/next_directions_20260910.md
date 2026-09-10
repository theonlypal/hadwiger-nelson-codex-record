# A direction family escaping the current additive five-coloring

This note records a construction proposal and exact exclusion of one class
of colorings. It does **not** establish a six- or seven-chromatic unit-distance
graph. In particular, a finite quotient with large chromatic number is not a
Euclidean realization of that quotient.

## Concrete next family

Let `t=exp(pi*i/21)` and

```
u=(2+t+3t²-2t⁴+t⁵+5t⁶+t⁷+t⁸-6t⁹-4t¹⁰+4t¹¹)/7.
U84 ={t^j,u*t^j:0<=j<42},
U126={t^j,u*t^j,u²*t^j:0<=j<42}.
```

All 126 directions are distinct and have exact complex norm one. The
implementation verifies these facts using rational arithmetic modulo Phi42.
The new orbit has

```
u²=(-2+t+4t²+3t³+4t⁴-4t⁵+t⁶+t⁷-t⁸+t¹⁰-3t¹¹)/7.
```

This orbit escapes the old additive lattice `L=Z[t]+u Z[t]`. For a point
whose coordinates have denominator dividing 7, write `7z=sum b_j t^j`.
The functional `b0-4b11 mod7` vanishes on every generator `t^j,u*t^j`
for `0<=j<12`. Consequently it vanishes throughout L. For u² it equals
`-2-4*(-3)=10=3 mod7`. This is a short exact certificate of nonmembership,
not a numerical separation. The 24 generator checks and the new residue
are stored in `runs/haugland-residue-exact-spectrum/residues.json`.

The proposed finite construction is the exact analogue of a short-path
region, now using U126: take all vertices on paths of at most six unit
steps between selected endpoints. Meet unordered half-path sums and include
all subset sums of every exact matching multiset. This generates every
ordering without introducing approximate equalities. Before using any pair
as a spindle gadget, separately certify its color relation at the relevant
palette size. Merely transplanting a four-color relation into a five- or
six-color argument is invalid.

## Independently checked additive-F5 obstruction

Reduce rational coefficients modulo 5; denominators of U126 are invertible.
There are 63 projective direction residues after identifying nonzero scalar
multiples. There is no F5-linear functional on the 12 coefficient coordinates
which is nonzero on every residue. The same holds in the 48-coordinate
tower, because all these directions have zero coordinates in its other
three parts.

Artifacts are in `runs/tower-odd-generator-minimal/haugland-126/`:

* `directions.json` contains the 48-coordinate residues;
  `directions-12.txt` contains their first twelve coordinates.
* `linear.cnf` has 2,475 variables and 15,424 clauses. Glucose4 produced
  `proof.drat`; the independent `drat-trim` replay in `drat-check.log`
  reports `s VERIFIED`.
* `tools/check_additive_five.c` provides a separate exhaustive check without
  a SAT encoding. All twelve coordinate directions occur, forcing all
  coefficients to be nonzero. Normalize the first coefficient to 1 and
  enumerate the remaining `4^11=4,194,304` assignments.
  `exhaustive-check.log` records zero survivors.

The C checker reads only the finite residue vectors; it imports no
construction arithmetic. It does not by itself certify that the supplied
vectors represent exact unit directions. That connection is checked by the
exact direction construction and residue script.

These results exclude additive coefficient colorings modulo5. They do not
exclude nonlinear periodic colorings, colorings using another modulus,
or arbitrary five-colorings of any generated graph.

## Why characteristic5 sees a new orbit

Modulo5, the cyclotomic polynomial factors as

```
f1=x^6+2x^4+3x^3+2x²+1,
f2=x^6+x^5+3x^4+4x³+3x²+x+1.
```

Both factors are irreducible. Each factor quotient is F15625; conjugation
acts as `a -> a^125`, and its fixed field is F125. The norm-one group has
126 elements. The image of t has order42 in both factors. In the first
factor, u has order9 and `u³=t^14`; in the second it has order18 and
`u³=t^35`. Thus the three cosets in U126 cover the entire norm-one group
in each projection. U84 covers only two cosets.

The two projections should not be confused with the full product ring.
Their different relations for u³ mean that U126 is not simply a single
126-element subgroup of the product. The relation `u^6=t^28` does hold
in both projections. Nor does coverage of the norm-one group alone prove
the full 12-coordinate additive obstruction; the exhaustive check above
establishes that separate fact.

The reduced characteristic3 quotient behaves differently:

```
Phi42 mod3=(x^6+2x^5+x^4+2x³+x²+2x+1)².
```

Quotienting by the irreducible degree-six factor gives F729. Here t has
order14, u has order4, and `u²=-1`. U84 and U126 both cover the same
28-element norm-one group for F729/F27. This statement is about the reduced
quotient: the full ring modulo Phi42 in characteristic3 retains a nilpotent
part and still distinguishes 126 residues.

## Small finite quotient screens and exact spectra

For a direction set S in Fp^6, the Cayley graph's additive-character
eigenvalues are `sum(s in S) exp(2*pi*i*(a dot s)/p)`, indexed by every
`a in Fp^6`. Counting each residue of the dot product computes them exactly.
No numerical eigenvalue calculation is needed.

For p=3, the reduced graph has 729 vertices and degree28. Its full spectrum,
written as `eigenvalue:multiplicity`, is

```
-8:112, -5:84, -2:168, 1:84,
4:168, 7:84, 10:28, 28:1.
```

Hoffman's bound gives `chi >= ceil(1+28/8)=5` for this finite graph.
Thirty-second Glucose runs at k=5 and k=6 both returned UNKNOWN; neither
timeout is a lower bound. The saved k=4 timeout is superseded by the exact
spectral bound.

For each p=5 factor, the finite graph has 15,625 vertices and degree126.
Its minimum eigenvalue is exactly `-9-5sqrt5`. All character-count tuples
and multiplicities are retained, and their ordering is checked by integer
comparisons of `a+b*sqrt5`. Hoffman's bound gives

```
chi >= ceil(1+126/(9+5sqrt5)) = 8.
```

For example the minimum comes from residue counts `(18,22,32,32,22)`;
the corresponding character sum is `-9-5sqrt5`. To check that the displayed
bound exceeds7, the needed strict inequality reduces to `144>125`.
Thus a coloring with at most seven colors of either single factor quotient
cannot provide a universal coloring of U126. This is a finite-field
obstruction to that method; it does not prove a planar chromatic lower bound.

The reduced p=3 rank-three linear quotient screen enumerated all 33,880
row spaces. Of these, 13,664 introduce a loop, and the remainder yield715
distinct Cayley graphs on27 vertices. All715 k=6 solver calls returned
UNSAT, but those calls were not independently proof-checked. Their only
use here is to report that this bounded search found no coloring witness.
The diagnostic is `runs/haugland-residue-20260910-v2/rank3-quotients.json`.

The exact factor, norm-group, lattice, and spectral checks replay with:

```
.venv/bin/python tools/haugland_residue_screen.py --out NEW_DIRECTORY
```

The canonical completed arithmetic run is
`runs/haugland-residue-exact-spectrum/`. The first attempted residue run
failed on a tuple-versus-list FLINT constructor call before producing
results; it remains preserved. The v2 run contains the finite graph and
bounded SAT attempts. Its later spectrum files were exploratory; the
canonical run includes exact sign comparisons and the final archived script.

## What the reviewed primary constructions actually provide

[Haugland v4, 17 August2026](https://arxiv.org/html/2608.04542v4)
uses precisely U84. Its short-path construction creates a non-monochromatic
pair for four colors, then a monochromatic pair, then a spindle proving
five colors necessary. I found no higher multiplicative powers of u used
in that paper. Its Section4 discusses at most six equivalence classes of
four-colorings; these are not six-colorings. The paper expressly leaves
global extension of its local four-colorings unproved. Its numerical
point-merging argument is not assumed by our exact constructors.

[Heule's2018 construction](https://arxiv.org/html/1805.12181)
supports choosing rotations that couple actual radius shells and then
examining forced palette patterns. It also reports a larger doubled
Minkowski-sum construction as five-colorable, despite costly solving. It
does not supply a certified non-five-colorable unit-distance gadget.

[Parts' odd-distance construction](https://arxiv.org/html/2206.12632)
does produce six-chromatic graphs when multiple odd distances are edges.
The paper's virtual-edge motivation concerns four-colorings, not an
established extension to five or six colors. Replacing those odd-distance
edges by unit-distance constraints would require separately certified
gadgets at the required palette size.

[Parts' two-distance construction](https://arxiv.org/html/2010.12656)
forbids both1 and the golden ratio. Its pentagon and spindle mechanisms
are useful examples of palette amplification, but their diagonal edges
are not unit edges. [Voronov's interval result](https://arxiv.org/html/2304.10163)
likewise concerns a positive interval of forbidden distances, not the
single exact distance1.

The additive-character spectral method for finite quadratic-form graphs
has a primary precedent in [Medrano, Myers, Stark and Terras,
*Finite analogues of Euclidean space*](https://www.sciencedirect.com/science/article/pii/0377042795002618).
The spectra used above are freshly enumerated integer character counts,
not an invocation of an asymptotic estimate or a quoted finite chromatic
number. In particular, the characteristic5 quotient is an anisotropic
norm graph over F125; it should not be silently identified with a split
sum-of-two-squares graph over that field.

## A palette limitation on further spindle amplification

If H is five-colorable and u,v are nonadjacent, start with a five-coloring
and recolor both vertices with an unused sixth color. This is a proper
six-coloring with equal endpoint colors. Recoloring just one vertex with
that sixth color also shows that no pair can be forced equal in every
six-coloring. Thus a five-colorable seed cannot supply a nontrivial forced
pair relation for six colors. A spindle route to seven must first obtain
an appropriate six-chromatic intermediate or establish a joint obstruction
not based on such a pair relation.

The current mixed Haugland–Heule seed and its first adaptive extension
both have explicit checked five-colorings; their repeated size growth is
therefore insufficient by itself. U126 removes one demonstrated algebraic
obstruction to the next search, while leaving the required finite geometric
and non-six-colorability certificates entirely outstanding.

## Executed mixed odd-rotation candidate

The next exact candidate retained the full `G=HauglandG3+V31` seed and
formed `G union uG`, rotating about the origin with the same exact u above.
`hn/tower_rotation_search.py` preserves both seed edge sets, adds proposed
cross-copy edges only after exact norm checks, and records copy membership.
This is an induced-graph-free certificate format: the saved edge list is
the graph, and enumeration completeness is not claimed.

`runs/mixed-g3-v31-odd-rotation-v1/` contains the completed execution:

* 125,929 distinct vertices, 991,423 edges, and55 shared vertices;
* 985 strict cross-copy edges, incident to1,309 vertices;
* all495,270 edges from each seed copy preserved, with102 common edges;
* twelve numerical near-unit proposals rejected by exact arithmetic.

The separate tower verifier accepted every coordinate and edge before
native solving. CaDiCaL returned a proper six-coloring in101.986 seconds.
The complete native model, all7,963,402 canonical CNF clauses, and the
edge coloring were independently checked. The model is
`six-color/coloring.json` and its report is `six-color/verification.json`.
The graph SHA256 is
`baff9f6fe24055f81ac01fe66fa7530f8c663a05ead97ad13b5b9c917dcfca9b`.
The completed partial solver proof was losslessly compressed with a checked
decompression hash. No five-color diagnostic was requested or run.

This candidate therefore does not satisfy the required non-six-colorability
condition, despite escaping the previously found additive five-color rule.
