# Retained graph core audit

Snapshot: 2026-09-10T16:54:28.462826+00:00. 31 graph files; 31 distinct structural hashes.

A graph is 6-colorable iff its 6-core is 6-colorable; removed vertices restore greedily with at most five previously colored neighbors.

Read-only graph audit. Geometry and UNSAT proofs are not verified here. Large degeneracy or a surviving 6-core does not prove non-6-colorability.

| Graph | Vertices / edges | Degeneracy | 6-core vertices / edges | Removed vertices | Saved proper colorings |
|---|---:|---:|---:|---:|---|
| artifacts/haugland_g1.json | 740 / 3985 | 7 | 740 / 3985 | 0 | none found in snapshot |
| artifacts/haugland_g3.json | 2131 / 12530 | 7 | 2131 / 12530 | 0 | 5, 6 |
| artifacts/heule_3877_mq.json | 3877 / 26814 | 9 | 3829 / 26574 | 48 | none found in snapshot |
| artifacts/heule_7753_sqrt2_mq.json | 7753 / 53628 | 9 | 7657 / 53148 | 96 | 6 |
| artifacts/heule_sqrt2_star_mq.json | 54271 / 468600 | 12 | 54271 / 468600 | 0 | none found in snapshot |
| runs/adaptive-g3-joint-centers-v1/graph.json | 5627 / 25653 | 7 | 3427 / 17732 | 2200 | 6 |
| runs/adaptive-g3-v1/round-000/graph.json | 2281 / 13742 | 7 | 2281 / 13742 | 0 | 6 |
| runs/adaptive-g3-v1/round-001/graph.json | 2431 / 15009 | 8 | 2431 / 15009 | 0 | 6 |
| runs/adaptive-g3-v1/round-002/graph.json | 2581 / 16270 | 8 | 2581 / 16270 | 0 | 6 |
| runs/adaptive-g3-v1/round-003/graph.json | 2731 / 17533 | 8 | 2731 / 17533 | 0 | 6 |
| runs/adaptive-g3-v1/round-004/graph.json | 2881 / 18844 | 8 | 2881 / 18844 | 0 | 6 |
| runs/adaptive-g3-v1/round-005/graph.json | 3031 / 20090 | 8 | 3031 / 20090 | 0 | 6 |
| runs/adaptive-g3-v1/round-006/graph.json | 3181 / 21377 | 8 | 3181 / 21377 | 0 | 6 |
| runs/adaptive-g3-v1/round-007/graph.json | 3329 / 22613 | 8 | 3329 / 22613 | 0 | 6 |
| runs/adaptive-g3-v1/round-008/graph.json | 3473 / 23833 | 8 | 3473 / 23833 | 0 | 6 |
| runs/adaptive-g3-v1/round-009/graph.json | 3609 / 25018 | 9 | 3609 / 25018 | 0 | 6 |
| runs/adaptive-g3-v1/round-010/graph.json | 3759 / 26254 | 9 | 3759 / 26254 | 0 | 6 |
| runs/adaptive-g3-v1/round-011/graph.json | 3909 / 27530 | 9 | 3909 / 27530 | 0 | 6 |
| runs/adaptive-mq-intersections-v1/round-000/graph.json | 14539 / 119031 | 10 | 14476 / 118737 | 63 | 6 |
| runs/adaptive-mq-intersections-v1/round-001/graph.json | 14689 / 121272 | 10 | 14644 / 121060 | 45 | none found in snapshot |
| runs/g3-hex-translations/graph.json | 11848 / 89898 | 10 | 11848 / 89898 | 0 | 6 |
| runs/g3-seven-rotations/graph.json | 14911 / 90404 | 8 | 14911 / 90404 | 0 | 6 |
| runs/g3-triangle-completion/graph.json | 13861 / 96917 | 9 | 13113 / 93297 | 748 | 6 |
| runs/mixed-ball2-2layers/graph.json | 14113 / 63294 | 8 | 13777 / 62370 | 336 | 5, 6 |
| runs/mixed-ball2-3layers/graph.json | 31753 / 137340 | 8 | 31249 / 135912 | 504 | 5, 6 |
| runs/mq-anchored-sqrt2-3/graph.json | 15505 / 107259 | 9 | 15313 / 106299 | 192 | 6 |
| runs/mq-intersections-heule-r1/graph.json | 28777 / 233550 | 10 | 28585 / 232638 | 192 | none found in snapshot |
| runs/mq-intersections-v1939-r1/graph.json | 14389 / 116742 | 10 | 14293 / 116286 | 96 | 6 |
| runs/t6-spindle/graph.json | 32241 / 229108 | 9 | 30647 / 222756 | 1594 | 6 |
| runs/verifier-mq-proof-selftest/triangle.json | 3 / 3 | 2 | 0 / 0 | 3 | 6 |
| runs/verifier-proof-selftest/triangle.json | 3 / 3 | 2 | 0 / 0 | 3 | none found in snapshot |

## Field-partition coupling

Partitions use coordinate support, not inferred copy identities. Counts include all exact stored edges.

| Graph | Part sizes | Cross edges / endpoints | 6-core cross edges / endpoints |
|---|---:|---:|---:|
| artifacts/haugland_g3.json | 1066 / 1065 | 60 / 63 | 60 / 63 |
| artifacts/heule_7753_sqrt2_mq.json | 3877 / 3876 | 60 / 61 | 60 / 61 |
| artifacts/heule_sqrt2_star_mq.json | 3877 / 50394 | 23346 / 27211 | 23346 / 27211 |
| runs/adaptive-g3-joint-centers-v1/graph.json | 1066 / 4561 | 3556 / 3680 | 1356 / 1395 |
| runs/adaptive-g3-v1/round-000/graph.json | 1139 / 1142 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-001/graph.json | 1211 / 1220 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-002/graph.json | 1287 / 1294 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-003/graph.json | 1361 / 1370 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-004/graph.json | 1434 / 1447 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-005/graph.json | 1507 / 1524 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-006/graph.json | 1583 / 1598 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-007/graph.json | 1652 / 1677 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-008/graph.json | 1734 / 1739 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-009/graph.json | 1810 / 1799 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-010/graph.json | 1890 / 1869 | 60 / 63 | 60 / 63 |
| runs/adaptive-g3-v1/round-011/graph.json | 1972 / 1937 | 60 / 63 | 60 / 63 |
| runs/adaptive-mq-intersections-v1/round-000/graph.json | 1963 / 12576 | 25075 / 14514 | 24949 / 14451 |
| runs/adaptive-mq-intersections-v1/round-001/graph.json | 1987 / 12702 | 25309 / 14664 | 25219 / 14619 |
| runs/g3-hex-translations/graph.json | 4393 / 7455 | 420 / 438 | 420 / 438 |
| runs/g3-seven-rotations/graph.json | 1066 / 13845 | 492 / 496 | 492 / 496 |
| runs/g3-triangle-completion/graph.json | 6929 / 6932 | 85 / 89 | 81 / 85 |
| runs/mixed-ball2-2layers/graph.json | 3529 / 10584 | 7224 / 7393 | 7140 / 7225 |
| runs/mixed-ball2-3layers/graph.json | 3529 / 28224 | 14364 / 14533 | 14280 / 14365 |
| runs/mq-anchored-sqrt2-3/graph.json | 3877 / 11628 | 51 / 57 | 51 / 57 |
| runs/mq-intersections-heule-r1/graph.json | 3877 / 24900 | 49800 / 28741 | 49416 / 28549 |
| runs/mq-intersections-v1939-r1/graph.json | 1939 / 12450 | 24900 / 14371 | 24708 / 14275 |
| runs/t6-spindle/graph.json | 16121 / 16120 | 168 / 253 | 86 / 89 |

## Observed bottlenecks

The degree-ranked Haugland adaptive rounds raise degeneracy and retain their added vertices in the 6-core, but the radical-field cut remains at 60 edges on 63 endpoints. Of those edges, 58 share one endpoint, leaving only two cross edges after deleting that hub. Growth within the two sides has not strengthened their coupling.

The 7,753-vertex diagonal sqrt2 union has one articulation vertex in its 6-core. Its 60 field-crossing edges all meet one hub, so this rotation joins the old and rotated copies at their common center without adding a distributed coupling.

The unit-circle intersection candidate is different: almost all vertices survive peeling and its field-crossing edges touch almost every core vertex. Its new points have two old neighbors each, so their many new-to-new edges are essential to their survival. The sqrt2 translated-star graph also retains all vertices and distributes its field crossings widely. These are stronger coupling mechanisms by the audited criteria, without implying a chromatic lower bound.

The joint Haugland mixed-center experiment confirms the same effect: adding all 3,496 screened centers creates 5,627 vertices and 25,653 edges. Although every new center initially had two old neighbors, 1,296 survive joint 6-core peeling. The enlarged core has 3,427 vertices and 17,732 edges; its field-crossing edges outside the original hub increase from 2 to 1,298. A saved proper six-coloring still rejects this graph as a solution. The individual-center degree filter would have discarded this entire coupling mechanism.

## Search criterion

Recompute the 6-core after each proposed augmentation. Reject an augmentation that only adds vertices removed by peeling and adds no edge between old retained vertices. For circle-intersection batches, include interactions between new vertices before peeling: degree two to the old graph alone cannot constrain six colors. Rank surviving batches by growth in core cross edges after removing the largest boundary hub, number of distinct boundary vertices, and coverage of multiple independently obtained six-colorings.

A concrete next construction is joint unit-circle intersection closure on the existing Heule core: create centers at squared separations 4/3 and 8/3, add all exact new-to-old and new-to-new unit edges, peel to the 6-core, and retain only components supported by at least six internal neighbors. Then enumerate a fresh six-coloring and score subsequent intersection batches by core connections to its six color classes. Preserve every finite graph and every SAT coloring; no score is a non-colorability certificate.
