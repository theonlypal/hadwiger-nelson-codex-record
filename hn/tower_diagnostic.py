"""Exact structural audit and bounded five-color diagnostic for G3 + V31.

Factor coordinates come from the recorded construction, but decomposition
sums are recomputed with the independent verifier's rational polynomials.
An edge is strictly cross-factor only when NO decomposition shares either
factor. This avoids confusing a different representation of the same point
with a genuine new coupling edge.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import json
from pathlib import Path
import time

from tools import verify_tower as verifier
from .mixed_cyclotomic import lift_haugland, v31_points
from .solve_external import solve_external


def dump(path, value):
    Path(path).write_text(json.dumps(value, separators=(",", ":"))+"\n")


def key(parts):
    return tuple(tuple(part.coeffs()) for part in parts)


def adjacency(n, edges):
    result = [[] for _ in range(n)]
    for u, v in edges:
        result[u].append(v); result[v].append(u)
    return result


def components(adj):
    seen = bytearray(len(adj)); sizes = []
    for root in range(len(adj)):
        if seen[root]:
            continue
        seen[root] = 1; queue = [root]; count = 0
        while queue:
            vertex = queue.pop(); count += 1
            for other in adj[vertex]:
                if not seen[other]:
                    seen[other] = 1; queue.append(other)
        sizes.append(count)
    return sorted(sizes, reverse=True)


def bottlenecks(adj):
    # Iterative Tarjan traversal avoids recursion depth assumptions.
    n = len(adj); discovery = [-1]*n; low = [0]*n; parent = [-1]*n
    children = [0]*n; articulations = set(); bridges = []; tick = 0
    for root in range(n):
        if discovery[root] >= 0:
            continue
        discovery[root] = low[root] = tick; tick += 1
        stack = [(root, iter(adj[root]))]
        while stack:
            vertex, neighbors = stack[-1]
            try:
                other = next(neighbors)
            except StopIteration:
                stack.pop(); previous = parent[vertex]
                if previous < 0:
                    if children[vertex] > 1:
                        articulations.add(vertex)
                else:
                    low[previous] = min(low[previous], low[vertex])
                    if low[vertex] > discovery[previous]:
                        bridges.append([min(previous, vertex), max(previous, vertex)])
                    if parent[previous] >= 0 and low[vertex] >= discovery[previous]:
                        articulations.add(previous)
                continue
            if discovery[other] < 0:
                parent[other] = vertex; children[vertex] += 1
                discovery[other] = low[other] = tick; tick += 1
                stack.append((other, iter(adj[other])))
            elif other != parent[vertex]:
                low[vertex] = min(low[vertex], discovery[other])
    return {"articulation_vertices": sorted(articulations), "bridges": sorted(bridges)}


def peel_core(adj, threshold):
    degree = [len(neighbors) for neighbors in adj]
    active = bytearray([1])*len(adj)
    queue = deque(v for v, d in enumerate(degree) if d < threshold)
    eliminated = []
    while queue:
        vertex = queue.popleft()
        if not active[vertex] or degree[vertex] >= threshold:
            continue
        neighbors = [other for other in adj[vertex] if active[other]]
        if len(neighbors) != degree[vertex]:
            raise ArithmeticError("core degree bookkeeping mismatch")
        active[vertex] = 0
        eliminated.append({"vertex": vertex, "remaining_neighbors": neighbors})
        for other in neighbors:
            degree[other] -= 1
            if degree[other] < threshold:
                queue.append(other)
    kept = [vertex for vertex, is_active in enumerate(active) if is_active]
    return kept, eliminated


def structural_audit(graph, out):
    out = Path(out)
    geometry = verifier.verify_geometry(graph)
    arithmetic = verifier.TowerArithmetic()
    points = [arithmetic.parse(vertex, index) for index, vertex in enumerate(graph["vertices"])]
    index = {key(point): i for i, point in enumerate(points)}
    factor_a_json = [point.to_json() for point in lift_haugland("g3")]
    factor_b_json = [point.to_json() for point in v31_points()]
    factor_a = [arithmetic.parse(point, i) for i, point in enumerate(factor_a_json)]
    factor_b = [arithmetic.parse(point, i) for i, point in enumerate(factor_b_json)]
    decompositions = [[] for _ in points]
    for a, left in enumerate(factor_a):
        for b, right in enumerate(factor_b):
            total = tuple(x+y for x, y in zip(left, right))
            if key(total) not in index:
                raise ArithmeticError("graph omits a point from the claimed G3 + V31 sum")
            decompositions[index[key(total)]].append([a, b])
    if any(not pairs for pairs in decompositions):
        raise ArithmeticError("graph contains a point without a G3 + V31 decomposition")
    dump(out/"factor-coordinates.json", {"field": graph["field"], "g3": factor_a_json, "v31": factor_b_json})
    dump(out/"decompositions.json", decompositions)
    offset_masks = [sum(1 << b for _, b in pairs) for pairs in decompositions]
    seed_sets = [set(a for a, _ in pairs) for pairs in decompositions]
    categories = Counter(); strict_cross = []; inherited = []; cross_degree = [0]*geometry.vertices
    layer_pair_counts = Counter(); inherited_per_offset = [0]*len(factor_b)
    for u, v in graph["edges"]:
        common_offsets = offset_masks[u] & offset_masks[v]
        common_seeds = bool(seed_sets[u] & seed_sets[v])
        if common_offsets and common_seeds:
            categories["inherited_from_both_factors"] += 1
        elif common_offsets:
            categories["inherited_from_g3_only"] += 1
        elif common_seeds:
            categories["inherited_from_v31_only"] += 1
        else:
            categories["strict_cross_factor"] += 1
            strict_cross.append([u, v]); cross_degree[u] += 1; cross_degree[v] += 1
            layer_pairs = {tuple(sorted((b, d))) for _, b in decompositions[u] for _, d in decompositions[v]}
            layer_pair_counts.update(layer_pairs)
        if common_offsets or common_seeds:
            inherited.append([u, v])
        for offset in range(len(factor_b)):
            if common_offsets & (1 << offset):
                inherited_per_offset[offset] += 1
    dump(out/"strict-cross-edges.json", strict_cross)
    adj = adjacency(geometry.vertices, graph["edges"])
    inherited_adj = adjacency(geometry.vertices, inherited)
    stats = {
        "goal_satisfied": False, "geometry": geometry.__dict__,
        "factor_sizes": [len(factor_a), len(factor_b)],
        "decomposition_count": sum(map(len, decompositions)),
        "decompositions_per_vertex_histogram": dict(Counter(map(len, decompositions))),
        "edge_categories": dict(categories),
        "strict_cross_definition": "No exact decomposition pair shares either the G3 seed point or the V31 offset",
        "inherited_definition": "Shared factor coordinate; edge has exact unit length, so inheritance does not depend on potentially incomplete factor edge lists",
        "inherited_g3_edges_per_offset": inherited_per_offset,
        "strict_cross_incident_vertices": sum(degree > 0 for degree in cross_degree),
        "strict_cross_degree_histogram": dict(Counter(cross_degree)),
        "strict_cross_layer_pair_counts": [{"offset_pair": pair, "edge_count": count} for pair, count in sorted(layer_pair_counts.items())],
        "layer_pair_count_note": "One edge may support several offset pairs when a vertex has multiple decompositions",
        "degree_histogram": dict(Counter(map(len, adj))),
        "minimum_degree": min(map(len, adj)), "maximum_degree": max(map(len, adj)),
        "average_degree": 2*geometry.edges/geometry.vertices,
        "full_component_sizes": components(adj), "inherited_only_component_sizes": components(inherited_adj),
        "full_graph_bottlenecks": bottlenecks(adj), "inherited_only_bottlenecks": bottlenecks(inherited_adj),
        "core_note": "Core sizes are structural diagnostics, not chromatic lower-bound certificates",
    }
    for threshold in (5, 6):
        kept, eliminated = peel_core(adj, threshold)
        remap = [-1]*geometry.vertices
        for new, old in enumerate(kept):
            remap[old] = new
        core_edges = [[remap[u], remap[v]] for u, v in graph["edges"] if remap[u] >= 0 and remap[v] >= 0]
        stats[f"core_{threshold}"] = {"vertices": len(kept), "edges": len(core_edges), "removed_vertices": len(eliminated)}
        dump(out/f"core-{threshold}-elimination.json", {"threshold": threshold, "core_to_original": kept,
                                                       "original_to_core": remap, "eliminated_in_order": eliminated})
        if threshold == 5:
            core_graph = {"schema": graph["schema"], "field": graph["field"],
                          "vertices": [graph["vertices"][old] for old in kept], "edges": core_edges,
                          "metadata": {"purpose": "five-core diagnostic only", "source": graph.get("metadata", {}),
                                       "elimination_map": "core-5-elimination.json", "goal_satisfied": False}}
            dump(out/"core-5-graph.json", core_graph)
    dump(out/"structure.json", stats)
    return stats


def extend_coloring(graph, core_coloring, elimination):
    colors = [None]*len(graph["vertices"])
    for core_index, original in enumerate(elimination["core_to_original"]):
        colors[original] = core_coloring["colors"][core_index]
    for record in reversed(elimination["eliminated_in_order"]):
        neighbor_colors = {colors[other] for other in record["remaining_neighbors"]}
        if None in neighbor_colors or len(record["remaining_neighbors"]) >= 5:
            raise ArithmeticError("invalid five-core coloring extension record")
        colors[record["vertex"]] = next(color for color in range(5) if color not in neighbor_colors)
    coloring = {"k": 5, "colors": colors}
    verifier.verify_coloring(graph, coloring, 5)
    return coloring


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", default="artifacts/mixed_haugland_g3_v31_tower.json", type=Path)
    parser.add_argument("--out", default="runs/mixed-g3-v31-diagnostic", type=Path)
    parser.add_argument("--timeout", default=300, type=float)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    graph = json.loads(args.graph.read_text())
    source_hash = verifier.common._sha256(args.graph)
    started = time.monotonic(); stats = structural_audit(graph, args.out)
    stats["source_graph_sha256"] = source_hash; stats["structural_audit_seconds"] = time.monotonic()-started
    dump(args.out/"structure.json", stats)
    print(json.dumps({"stage": "structural_audit_complete", "edge_categories": stats["edge_categories"],
                      "core_5": stats["core_5"], "core_6": stats["core_6"], "minimum_degree": stats["minimum_degree"]}), flush=True)
    result = solve_external(args.out/"core-5-graph.json", args.out/"native-k5", timeout=args.timeout, configuration="sat", k=5)
    report = {"goal_satisfied": False, "source_graph_sha256": source_hash,
              "native_status": result["status"], "scope": "five-color diagnostic, not the non-six target"}
    if result["status"] == "SAT":
        core_coloring = json.loads((args.out/"native-k5/coloring.json").read_text())
        elimination = json.loads((args.out/"core-5-elimination.json").read_text())
        full_coloring = extend_coloring(graph, core_coloring, elimination)
        full_coloring["graph_sha256"] = source_hash
        dump(args.out/"full-five-coloring.json", full_coloring)
        geometry = verifier.verify_geometry(graph)
        report.update({"certificate_status": "exact-geometry-and-extended-proper-5-coloring",
                       "geometry": geometry.__dict__, "coloring_sha256": verifier.common._sha256(args.out/"full-five-coloring.json")})
    elif result["status"] == "NON5_VERIFIED":
        report.update({"certificate_status": "verified-non5-core-only", "chromatic_lower_bound": 6})
    else:
        report["certificate_status"] = "no-colorability-conclusion-from-native-diagnostic"
    if verifier.common._sha256(args.graph) != source_hash:
        raise ArithmeticError("source graph changed during diagnostic")
    dump(args.out/"summary.json", report)
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
