# Active research continuation

The user explicitly requires continued work until a finite exact unit-distance
graph has an independently checkable non-6-colorability certificate. The
Codex goal is active. **No target certificate has been found.**

## Established evidence

The original seven graphs and twelve palette-guided Haugland augmentations
have explicit proper6colorings. Increasing graph size alone did not fix their
small cross-field cut. `research/core_audit.md` and
`research/coloring_pattern_analysis.json` describe those checks.

The joint Haugland center experiment `runs/adaptive-g3-joint-centers-v1`
adds3496 exact centers, reaching5627vertices/25653edges.1296new centers survive
the6core, and the cross-field cut excluding the original hub grows2→1298.
It still has an independently checked6coloring. Do not confuse a lack of
single-center support with a lack of useful joint support.

A new exact real field implementation covers
`Q(sqrt2,sqrt3,sqrt5,sqrt11)` in `hn/multiquadratic.py`, with a separate
independent verifier `tools/verify_mq.py`. It reproduces Heule's3877vertex
union and constructs circle-intersection closures. The V1939 closure at
squared separation4/3 has14389vertices/116742edges; its6core has14293vertices.
Glucose timed out at6colors, but nativeCaDiCaL found and independently
verified a6coloring. A54271vertex/468600edge sqrt2star also timed out before
the stronger algebraic obstruction below resolved it.

## A proved family exclusion: the denominator-5 boundary

All currently checked multiquadratic candidates whose rational basis
coefficients have denominators coprime to5 have explicit **5-colorings**.
This is a constructive exclusion, not an inference from solver behavior.

In F25=F5[w]/(w²-2), map the positive real square-root generators by

```
sqrt2 -> w, sqrt3 -> 2w, sqrt5 -> 0, sqrt11 -> 1.
```

These images satisfy every defining relation. Invert denominators modulo5.
For coordinates x=x0+x1*w and y=y0+y1*w, use color
`(x0+y0+y1) mod5`. Independent exhaustive enumeration checks all625
displacements, including all24unit displacements; none of those unit
displacements has color0. The function therefore colors every graph in
`Z_(5)[sqrt2,sqrt3,sqrt5,sqrt11]^2` with5colors.

`runs/mod5-family-elimination/index.json` records independently checked
concrete graph colorings, including the two previously timed-out
intersection graphs, the54271vertex star, anchored rotations, and the
adaptive MQ rounds. Original SAT timeouts remain preserved. Redundant
MQ SAT processes were terminated only after a concrete5coloring was checked.
Their raw proof streams are being losslessly compressed with decompressed
hash validation.

This does **not** prove5colorability of the whole number field: denominator5
poles prevent this reduction. Merely adjoining more ordinary square roots
does not escape F25 when coefficient denominators remain prime to5. The
Haugland cyclotomic42 field also escapes this particular residue map, since
F25 has no primitive seventh root of unity.

## Denominator-five execution

Build exact unit-circle intersections on the V1939 first closure (which
already contains sqrt2), choosing actual squared center separations that
introduce genuine denominator5 poles:

| s | sqrt((4-s)/s) |
|---|---|
|5/4|sqrt55/5|
|5/2|sqrt15/5|
|10/3|sqrt5/5|
|20/9|2sqrt5/5|
|25/9|sqrt11/5|

For centers a,b with squared separation s, the two points are
`(a+b)/2 ± i*(b-a)*sqrt((4-s)/s)/2`. Check both distances exactly and count
uncancelled denominator5 coefficients. Prefer balanced complete60degree
orbits and verify joint6core/cross-field coupling. Arbitrary denominator5
rotated copies tested so far create no useful exclusive cross-copy edges;
do not repeat those disjoint-copy tests.

`runs/mq-pole5-v1` contains the executed construction: 39,997 vertices and
220,230 exact edges. Complete 60-degree orbits retain 25,608 new pole
vertices. The 6-core retains 3,324 pole vertices, so the extension creates
real coupling. Native six-color search found a checked coloring in 82.63
seconds. Its 5-core has 28,297 vertices and 181,518 edges. The separate
`runs/mq-pole5-five-core-v1` diagnostic found a proper 5-coloring in 276.73
seconds; independently reversing the checked peeling order extends it to
the entire graph. This candidate therefore does not improve the lower bound.

## Mixed cyclotomic tower and direction exclusions

`hn/cyclotomic_tower.py` implements the exact degree-48 field
`Q(zeta_42)(sqrt5,sqrt33)`. The independent verifier is
`tools/verify_tower.py`; its 13 additional tests bring the verifier suite
to 46 tests at that stage. `hn/mixed_cyclotomic.py` forms exact Minkowski sums with
Heule's 31-point unit-direction configuration.

| Seed | Vertices | Edges | Native six-color result |
| --- | ---: | ---: | --- |
| Haugland G1 + V31 | 21,034 | 155,179 | SAT, independently verified |
| Haugland G3 + V31 | 62,992 | 495,270 | SAT, independently verified |

The second graph uses 198 signed unit directions. Reducing the 48 rational
basis coefficients modulo 5 gives 96 projectively distinct directions.
`tools/direction_color_screen.py` found a linear functional to F5 that is
nonzero on every one. Its explicit coefficients and full proper 5-coloring
are in `runs/mixed-haugland-g3-v31-linear5`.

This is a direction-set exclusion: it colors every translate of these
directions, but it is not a claim about all unit directions in the field.
Coloring-guided augmentations are in `runs/adaptive-mixed-tower-v2`.
Round 000 has 63,742 vertices and 502,532 edges and also passes independent
geometry and the explicit 5-coloring. Any further augmentation must be
checked against this functional before investing in SAT.

Current work seeks exact new unit directions and cross-copy edges that
escape the saved additive colorings, plus color-guided pole extensions.
An inability to find an additive coloring is only a search filter; it is
never evidence of non-5- or non-6-colorability.

`hn/color_filter.py` now verifies every exported edge against a supplied
five-color functional before allowing an ansatz-based rejection. It checks
exact geometry independently, never marks the goal satisfied, treats poles
as inapplicability, and checks all new edges rather than just the direction
pool. The full test suite currently has 57 passing tests, including the
later phase-driver and Kempe-witness checks.

## New odd-generator orbit

The exact unit directions `U126={u1^k*t^j:0<=k<3,0<=j<42}` include one
42-direction orbit absent from Haugland's published U84 construction.
This direction set admits no additive coefficient functional to F5;
the small ansatz CNF and independently verified DRAT proof are in
`runs/tower-odd-generator-minimal/haugland-126`. This is explicitly NOT
a proof about arbitrary graph colorings.

`hn/odd_power_search.py` constructs finite candidates. The radius-two ball
has 7,939 vertices and 38,682 edges; native solvers independently verified
proper six- and five-colorings (the five-color solve took 112.10 seconds).
The exact endpoint-walk variant enumerates 2,991 step multisets of length
at most six from 0 to i*sqrt3, producing 31,750 path vertices. Applying the
two-stage spindle construction gives 82,041 vertices and 604,354 edges in
`runs/odd-power-t6-spindle-v1`; its unrestricted target SAT run timed out
after 240 seconds. The partial binary proof is retained as a hash-verified
gzip stream. This is UNKNOWN, not a lower-bound result.

Finite residue graphs are being studied only as upper-bound filters:
a checked coloring may lift along an edge-preserving map. A finite
quotient's non-colorability never establishes non-colorability of a lift.

The combined direction experiment `runs/mixed-g3-v31-odd-rotation-v1`
forms `G union u1*G` from the mixed tower seed. It has 125,929 vertices,
991,423 edges, 55 shared vertices, and 985 strict cross-copy edges.
Native six-color SAT completed in 101.99 seconds; geometry, all graph edges,
and all 7,963,402 CNF clauses were independently checked. A phase-guided
five-color diagnostic timed out after 300.24 seconds. Neither a new direction orbit nor
failure of an additive coloring is a chromatic lower bound.

## Current proof-guided work

`tools/cadical_phase.cpp` is a thin API driver linked against the same pinned
CaDiCaL library. Its only clause input is the original DIMACS file. It calls
`phase(literal)` for preferences, never `assume` or `add` for hints. Tests
confirm that a false hint does not change satisfiability, and an emitted
UNSAT proof verifies against the original CNF. Driver binary SHA256:
`4c36bb0431f7daa3043bea79368be853d1bcb947ec3ea82bea670e75e1c57891`.

The `runs/odd-power-t6-phased-v1` run maps all 32,241 vertices of the
old T6 spindle's verified six-coloring into the new 82,041-vertex graph.
Its 193,446 signed phase preferences do not change the canonical CNF. The
run timed out after 600.24 seconds. A second phase run based on an improper
Tabu model timed out after 180.18 seconds. Their proof streams remain
losslessly archived. Algebraic sixth-color repair failed only for its
specified ansatz and did not establish non-colorability.

The candidate is now **eliminated by a proper six-coloring** in
`runs/odd-power-t6-hex-repair-v1`. `tools/hex_seed_repair.py` checks geometric
seven-color starts, removes one color by minimizing each affected vertex's
conflicts, and lets every vertex move in local search. The selected start
had 745 conflicts; 1,528 Tabu moves eliminated them in 0.012 seconds.
The separate exact verifier checked all 82,041 vertices and 604,354 edges,
and checked the complete six-color assignment. All older UNKNOWN results
remain unchanged, but they no longer motivate search on this graph.

The mixed seed has 12,799 vertex pairs with identical coefficient residues
modulo 5. Twenty-six additive colorings span all 26 active coefficient
columns and separate every other pair. Independently replayed Kempe swaps
separate 10,474 of the residue-collision pairs. A phase-guided unrestricted
conditional five-color solve also separated the unique remaining pair at
distance squared 25/3. Two small nonlinear offset colorings separate all
remaining pairs. Independent end-to-end verification covers all
1,983,964,536 distinct vertex pairs in
`runs/mixed-five-nonlinear-offset-completion/independent-complete-pair-verification.json`.
Thus this graph has no forced-equality pair in proper five-colorings.
The native conditional run was stopped after the full positive witness
was checked; its original UNKNOWN and partial proof are preserved.

The pole branch now also includes a 40,679-vertex/227,005-edge candidate
containing every point and edge of the pinned Heule 510 realization. Its
six-coloring is verified. The 529-center predecessor has a verified full
five-coloring; adding points that block one model did not prove an
obstruction to all models. A fresh four-color proof attempt on Heule 510
timed out, so its published lower bound is kept separate from fresh checks.

## Voronov independent-radical construction

Primary sources are arXiv:2106.11824v4 and `vsvor/dist-graphs` revision
`e3714d0a156f6ed151d4521debb2b89ce4f1075c`. The locally generated 73-,
865-, and 32,257-point base sets match all primary-source coordinates
exactly; see `research/voronov_source_audit.json`.

`hn/voronov_field.py` represents the real degree-16 field
Q(sqrt2,sqrt3,a,b), where
`a^2=5+3sqrt2-sqrt3+3sqrt6` and
`b^2=5+3sqrt2+sqrt3-3sqrt6`. Mixed signs in the real embeddings prove
independence of these square classes. The two roots jointly prevent
every extension of the prior F25 coloring map; that is an upper-color
filter exclusion, never a lower-bound certificate.

`tools/verify_voronov.py` independently implements arithmetic and the
exact field-degree check. Arithmetic tests and 62 relevant verifier
regressions pass, and native dispatch passed a complete positive control.

`runs/voronov-joint-1-7-v1` combines the 32,257-point base with two rotated
copies. The graph has 96,769 vertices and 814,008 exact edges, including
480 strict cross-copy edges. Its SHA256 is
`5ff891f60bb37831e9a6ad112227409d07fe4adcbae58f1406e466fd2919942c`.
The two nonzero radical branches have no direct unit edges; both interact
through the shared base. Native six-color SAT timed out after 300.04 seconds,
but `runs/voronov-joint-1-7-hex-repair-v1` subsequently supplied a verified
proper six-coloring. Its local repair needed 2,486 moves in 0.016 seconds.
The original proof stream remains archived. The global additive F5
coefficient ansatz is UNSAT, so any five-coloring must escape that ansatz;
no chromatic lower bound follows.

Current construction work tests nonzero rotation centers and exact new
vertices that block the saved six-color assignment. Every new graph is
first subjected to the checked geometric-start coloring filter before a
long native proof search. Blocking one saved coloring never means that
all colorings are blocked. All 87 current tests pass.

## Solver and certificate tools

`hn/solve.py`: PySATGlucose4, unrestricted canonical CNF, saved models and
run diagnostics. Original timeout runs did not retain their partial proof
streams; the newer native wrapper does.

`hn/solve_external.py`: pinned nativeCaDiCaL3.0.1, source commit
`c60730422e758ef1cebe7aeddf2dda31c996bf04`, binary SHA256
`39143d3955b59965d9e8ec2b7df71206fa4380d6f17c3c8e3254dce653897fd0`.
It validates geometry before execution, preserves complete CNF/log/model and
streamed binaryDRAT, independently checks the full model on SAT, and invokes
the exact geometry/non6DRAT gate on UNSAT. A diagnostic5color UNSAT must
never set the target goal true. The independentDRAT checker is separately
pinned in `sources/README.md`.

For exactMQ or cyclotomic graphs, a target witness must pass the corresponding
independent verifier on the same graph and proof CNF. Never use a finite-field
graph, a virtual edge, a timeout, a5color obstruction, or an intermediate
6chromatic graph as the requested success certificate.
