# Upper-color filters for the odd-power T6 spindle

The exact graph `runs/odd-power-t6-spindle-v1/graph.json` has 82,041 vertices
and 604,354 recorded unit edges. The experiments below attempt to obtain
explicit proper six-colorings. Failure of any of these restricted forms
does **not** prove that the graph is not six-colorable.

## Additive five-coloring plus one independent color class

The initial map is the first 24 entries of `linear_functional` in
`runs/mixed-haugland-g3-v31-linear5/result.json`. Additional maps come from
`runs/mixed-five-forced-pair-screen/functional-sampling.json`. Their first
24 entries act on the coefficient order

`a_0,...,a_11,b_0,...,b_11`

for `sum_j a_j*zeta42^j + sqrt(5)*sum_j b_j*zeta42^j`. Rational coefficients
are reduced as numerator times inverse denominator modulo five. Every
coefficient denominator in this graph is invertible modulo five.

For tentative colors `L(v)` in `{0,1,2,3,4}`, let `s_v` indicate that vertex
`v` is moved to color 5. Proper coloring requires the following 2-CNF:

- For every graph edge `uv`, `not s_u or not s_v`.
- For every edge with `L(u)=L(v)`, `s_u or s_v`.

Thus the selected vertices must form an independent set covering every
monochromatic edge. `tools/sixth_color_repair.py` constructs the implication
graph and computes its strongly connected components. A successful model
would be checked on the entire original graph with the independent exact
geometry and coloring verifier.

All 20 tested maps failed this repair. Their monochromatic edge counts range
from 23,565 to 42,961. The contradictions have small finite witnesses:

| Witness | Maps |
| --- | ---: |
| Monochromatic triangle | 2 |
| Monochromatic cycle of length five | 11 |
| Implication paths `s_v -> not s_v -> s_v` | 7 |

For a monochromatic edge, exactly one endpoint must be selected. An odd
cycle makes that alternation impossible. In the other seven cases, the bad
graph is bipartite, but independence constraints between its components
produce contradictory implications. Each saved path has only 6–10
literals.

The separate verifier `tools/verify_sixth_repair_obstructions.py` replayed
all 20 witnesses, checking the exact unit distances and rational modulo-five
reductions on every invoked edge. Its report is
`runs/odd-power-t6-sixth-repair-v2/independent-obstruction-verification.json`.
The run retains each functional, compressed full 2-CNF with hashes, and its
small witness. The earlier `v1` directory preserves an initialization error;
the completed screen is `v2`.

## Old-direction cosets modulo seven

Let `D=U84 union rho*U84`, where `rho` is the existing spindle rotation.
Multiply all 24 coordinate coefficients of `z` and every `d in D` by seven,
then reduce modulo seven. All resulting denominators are invertible. The
vectors `7*d` span a rank-two subspace of `F7^24`.

Let `H` be a row basis of the 22-dimensional annihilator of that span. The
value `h(z)=H(7*z) mod7` is unchanged along old-direction edges. Exactly
13 values occur in this graph: one coset has 34,557 vertices, six have 3,955,
and six have 3,959.

The quotient identifies vertices with the same pair `(h(z),L(z))`. A proper
six-coloring of the quotient would pull back to the original graph. All
20 maps give 65 types, but each quotient has loops, arising from 12–32
original edges whose endpoints share both values. This rejects those
quotient colorings immediately. The matrices, coset assignments, type maps,
and loop witnesses are in `runs/odd-power-t6-coset-repair-v1/`.

As a bounded refinement, each endpoint of those initial loop edges was given
its own singleton type. This creates 89–129 types and removes all loops.
The 20 resulting six-color quotient formulas all returned UNSAT in a total
of 3.722 seconds. Their CNFs and raw DRAT traces are retained under
`runs/odd-power-t6-coset-split-repair-v1/`; those traces were not independently
checked and the status is explicitly `QUOTIENT_UNSAT_UNVERIFIED`.

Even a checked quotient UNSAT proof would establish only failure of that
constant-on-types coloring form. It would not be a lower-bound certificate
for the planar graph. These upper filters produced no full coloring and no
non-six-colorability certificate.

## Bounded TabuCol and local SAT repair

`tools/tabu_color.cpp` is an untrusted candidate generator driven by
`hn/solve_local.py`. It initializes the 32,241 exact overlapping vertices
from the saved smaller graph coloring, greedily colors the others in
descending degree order, and then leaves every vertex free. It maintains
neighbor-color counts and the conflicted-vertex list. Each move samples up
to 64 conflicted vertices and uses tabu tenure with aspiration for a new
global best. The seed, input bytes, compiler command/version, source,
binary, candidate states, and observed conflict counts are retained.

Both runs stopped at the two-million-move cap before the 120-second limit:

| Seed | Initial conflicts | Best conflicts | Moves | Native seconds |
| --- | ---: | ---: | ---: | ---: |
| 0 | 5,809 | 1,005 | 2,000,000 | 3.868 |
| 1 | 5,675 | 1,008 | 2,000,000 | 3.678 |

The wrapper recomputed every conflicting edge of each saved best state.
Both results are `UNKNOWN_IMPROPER_BEST`; neither is a proper coloring or
evidence of a chromatic lower bound. Files are in
`runs/odd-power-t6-tabu-seed0-v1/` and
`runs/odd-power-t6-tabu-seed1-v1/`.

`hn/local_coloring_repair.py --allow-improper-base` then tried the better
state. It initially frees both endpoints of every conflicting edge, so all
remaining fixed boundary vertices already form a proper coloring. Growing
the free region produced the following results:

| Neighborhood radius | Free vertices | Restricted six-color result |
| --- | ---: | --- |
| 0 | 1,888 | UNSAT under fixed boundary colors |
| 1 | 43,073 | UNSAT under fixed boundary colors |
| 2 | 81,322 | Timeout at 20.897 seconds |

The exact restricted formulas and free-vertex maps are retained in
`runs/odd-power-t6-tabu-hybrid-repair-v1/`. These outcomes establish no lower
bound for the original graph. A positive repair, had one been found, would
have been independently checked on the full exact graph.
