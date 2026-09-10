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

## Current priority

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
