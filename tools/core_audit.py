#!/usr/bin/env python3
"""Read-only combinatorial audit of retained exact-graph JSON artifacts.

This audits graph structure, not geometry. It never turns solver timeouts,
large minimum degree, or missing saved colorings into chromatic claims.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import heapq
import json
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def core_numbers(adjacency):
    """Minimum-degree peeling; running maximum gives each vertex's core."""
    degree = [len(a) for a in adjacency]
    heap = [(d, i) for i, d in enumerate(degree)]
    heapq.heapify(heap)
    removed = [False] * len(adjacency)
    cores = [0] * len(adjacency)
    degeneracy = 0
    while heap:
        d, v = heapq.heappop(heap)
        if removed[v] or degree[v] != d:
            continue
        degeneracy = max(degeneracy, d)
        removed[v] = True
        cores[v] = degeneracy
        for w in adjacency[v]:
            if not removed[w]:
                degree[w] -= 1
                heapq.heappush(heap, (degree[w], w))
    return cores


def components(adjacency, selected=None):
    remaining = set(range(len(adjacency))) if selected is None else set(selected)
    sizes = []
    while remaining:
        start = remaining.pop()
        stack, size = [start], 1
        while stack:
            v = stack.pop()
            for w in adjacency[v]:
                if w in remaining:
                    remaining.remove(w)
                    stack.append(w)
                    size += 1
        sizes.append(size)
    return sorted(sizes, reverse=True)


def cut_structure(adjacency, selected):
    """Iterative Tarjan bridges and articulation vertices in an induced core."""
    selected = set(selected)
    discovery, low, parent, child_count = {}, {}, {}, Counter()
    cuts, bridges = set(), 0
    clock = 0
    for start in sorted(selected):
        if start in discovery:
            continue
        discovery[start] = low[start] = clock
        clock += 1
        parent[start] = None
        stack = [(start, iter(adjacency[start]))]
        while stack:
            v, iterator = stack[-1]
            try:
                w = next(iterator)
            except StopIteration:
                stack.pop()
                p = parent[v]
                if p is None:
                    if child_count[v] > 1:
                        cuts.add(v)
                else:
                    low[p] = min(low[p], low[v])
                    if low[v] > discovery[p]:
                        bridges += 1
                    if parent[p] is not None and low[v] >= discovery[p]:
                        cuts.add(p)
                continue
            if w not in selected or w == parent[v]:
                continue
            if w not in discovery:
                parent[w] = v
                child_count[v] += 1
                discovery[w] = low[w] = clock
                clock += 1
                stack.append((w, iter(adjacency[w])))
            else:
                low[v] = min(low[v], discovery[w])
    return {"articulation_vertices": len(cuts), "bridges": bridges}


def field_partition(graph, adjacency, core):
    if graph.get("schema") == "hn-exact-v1":
        labels = [int(any(Fraction(x) != 0 for x in p["b"])) for p in graph["vertices"]]
        label = "zero radical part / nonzero radical part"
    elif graph.get("schema") == "hn-multiquadratic-v1":
        try:
            bit = 1 << graph["radicands"].index(2)
        except ValueError:
            return None
        labels = [int(any(Fraction(x) != 0 for axis in ("x", "y")
                          for m, x in enumerate(p[axis]) if m & bit))
                  for p in graph["vertices"]]
        label = "sqrt2-free / nonzero sqrt2 component"
    else:
        return None
    cross = [(u, v) for u, v in graph["edges"] if labels[u] != labels[v]]
    cross_degree = Counter(v for edge in cross for v in edge)
    core_cross = [(u, v) for u, v in cross if u in core and v in core]
    core_cross_degree = Counter(v for edge in core_cross for v in edge)
    return {
        "definition": label,
        "partition_vertices": [labels.count(0), labels.count(1)],
        "cross_edges": len(cross),
        "cross_edge_endpoints": len(cross_degree),
        "cross_degree_histogram_nonzero": dict(sorted(Counter(cross_degree.values()).items())),
        "largest_cross_degree": max(cross_degree.values(), default=0),
        "cross_edges_after_removing_one_maximum_endpoint": len(cross)-max(cross_degree.values(), default=0),
        "six_core_partition_vertices": [sum(labels[v] == i for v in core) for i in (0, 1)],
        "six_core_cross_edges": len(core_cross),
        "six_core_cross_edge_endpoints": len({v for edge in core_cross for v in edge}),
        "six_core_cross_edges_after_removing_one_maximum_endpoint": len(core_cross)-max(core_cross_degree.values(), default=0),
        "interpretation_limit": "coordinate-field partition; not a proof of graph-copy membership",
    }


def color_audit(obj, graph, adjacency, cores):
    colors, k = obj.get("colors"), obj.get("k")
    if type(k) is not int or not isinstance(colors, list) or len(colors) != len(adjacency):
        return {"valid": False, "reason": "malformed coloring"}
    if any(type(c) is not int or not 0 <= c < k for c in colors):
        return {"valid": False, "reason": "out-of-range color"}
    if any(colors[u] == colors[v] for u, v in graph["edges"]):
        return {"valid": False, "reason": "monochromatic edge"}
    coverage = Counter(len({colors[w] for w in adjacency[v]}) for v in range(len(adjacency)))
    core_coverage = Counter(len({colors[w] for w in adjacency[v] if cores[w] >= 6})
                            for v in range(len(adjacency)) if cores[v] >= 6)
    return {"valid": True, "k": k,
            "class_sizes": [colors.count(c) for c in range(k)],
            "six_core_class_sizes": [sum(cores[v] >= 6 and colors[v] == c
                                          for v in range(len(adjacency))) for c in range(k)],
            "distinct_neighbor_colors_histogram": dict(sorted(coverage.items())),
            "six_core_distinct_neighbor_colors_histogram": dict(sorted(core_coverage.items())),
            "locally_forced_vertices_given_other_colors": coverage.get(k-1, 0),
            "interpretation_limit": "local color choices in this one saved coloring only"}


def analyze_graph(graph, witnesses):
    n = len(graph["vertices"])
    adjacency = [set() for _ in range(n)]
    for edge in graph["edges"]:
        if len(edge) != 2 or any(type(v) is not int or not 0 <= v < n for v in edge):
            raise ValueError("malformed graph edge")
        u, v = edge
        if u == v or v in adjacency[u]:
            raise ValueError("loop or duplicate edge")
        adjacency[u].add(v)
        adjacency[v].add(u)
    core = core_numbers(adjacency)
    selected = {v for v, c in enumerate(core) if c >= 6}
    core_edges = sum(u in selected and v in selected for u, v in graph["edges"])
    internal_degrees = [sum(w in selected for w in adjacency[v]) for v in selected]
    degree = [len(a) for a in adjacency]
    result = {
        "vertices": n, "edges": len(graph["edges"]),
        "minimum_degree": min(degree, default=0), "maximum_degree": max(degree, default=0),
        "average_degree": sum(degree)/n if n else 0,
        "degeneracy": max(core, default=0),
        "core_number_histogram": dict(sorted(Counter(core).items())),
        "degree_histogram": dict(sorted(Counter(degree).items())),
        "connected_component_sizes": components(adjacency),
        "six_core": {"vertices": len(selected), "edges": core_edges,
                     "vertex_fraction": len(selected)/n if n else 0,
                     "removable_vertices": n-len(selected),
                     "average_degree": 2*core_edges/len(selected) if selected else 0,
                     "minimum_degree": min(internal_degrees, default=0),
                     "component_sizes": components(adjacency, selected),
                     **cut_structure(adjacency, selected)},
        "field_partition": field_partition(graph, adjacency, selected),
        "saved_colorings": [],
        "six_colorability_from_degeneracy": max(core, default=0) <= 5,
    }
    for path, obj in witnesses:
        result["saved_colorings"].append({"path": path, **color_audit(obj, graph, adjacency, core)})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    paths = sorted(p for directory in (root/"artifacts", root/"runs")
                   if directory.exists() for p in directory.rglob("*.json"))
    objects, colorings, skipped = [], {}, []
    for path in paths:
        # Only JSON is read: CNF/DRAT/progress logs are never loaded.
        try:
            data = path.read_bytes()
            obj = json.loads(data)
        except (OSError, json.JSONDecodeError) as exc:
            skipped.append({"path": str(path.relative_to(root)), "reason": type(exc).__name__})
            continue
        if not isinstance(obj, dict):
            continue
        if "colors" in obj and "k" in obj and "graph_sha256" in obj:
            colorings.setdefault(obj["graph_sha256"], []).append((str(path.relative_to(root)), obj))
        elif obj.get("schema") in ("hn-exact-v1", "hn-multiquadratic-v1") and "vertices" in obj and "edges" in obj:
            structural = {k: obj[k] for k in ("schema", "field", "radicands", "vertices", "edges") if k in obj}
            structural_hash = digest(json.dumps(structural, sort_keys=True, separators=(",", ":")).encode())
            objects.append((path, digest(data), structural_hash, obj))
    grouped = {}
    for path, file_hash, structure_hash, obj in objects:
        group = grouped.setdefault(structure_hash, {"graph": obj, "files": [], "colorings": []})
        group["files"].append({"path": str(path.relative_to(root)), "sha256": file_hash})
        group["colorings"].extend(colorings.get(file_hash, ()))
    records = []
    for structure_hash, group in grouped.items():
        print("Auditing " + group["files"][0]["path"], flush=True)
        record = {"structure_sha256": structure_hash, "files": group["files"],
                  "schema": group["graph"]["schema"],
                  **analyze_graph(group["graph"], group["colorings"])}
        records.append(record)
    report = {
        "schema": "hn-core-audit-v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "graphs_found": len(objects), "distinct_structural_hashes": len(records),
        "scope": "retained graph JSON under artifacts/ and runs/ at this snapshot",
        "deduplication": "SHA256 of schema, field or radicands, ordered vertices and edges; excludes metadata",
        "theorem": "A graph is 6-colorable iff its 6-core is 6-colorable; removed vertices restore greedily with at most five previously colored neighbors.",
        "limits": "Read-only graph audit. Geometry and UNSAT proofs are not verified here. Large degeneracy or a surviving 6-core does not prove non-6-colorability.",
        "skipped_unreadable_json": skipped, "graphs": records,
    }
    out = root/"research"/"core_audit.json"
    out.write_text(json.dumps(report, indent=2)+"\n")
    lines = ["# Retained graph core audit", "", f"Snapshot: {report['created_utc']}. "
             f"{len(objects)} graph files; {len(records)} distinct structural hashes.", "",
             report["theorem"], "", report["limits"], "",
             "| Graph | Vertices / edges | Degeneracy | 6-core vertices / edges | Removed vertices | Saved proper colorings |",
             "|---|---:|---:|---:|---:|---|"]
    for row in records:
        c = row["six_core"]
        ks = sorted({x["k"] for x in row["saved_colorings"] if x.get("valid")})
        lines.append(f"| {row['files'][0]['path']} | {row['vertices']} / {row['edges']} | "
                     f"{row['degeneracy']} | {c['vertices']} / {c['edges']} | {c['removable_vertices']} | "
                     f"{', '.join(map(str,ks)) or 'none found in snapshot'} |")
    lines.extend(["", "## Field-partition coupling", "",
                  "Partitions use coordinate support, not inferred copy identities. Counts include all exact stored edges.", "",
                  "| Graph | Part sizes | Cross edges / endpoints | 6-core cross edges / endpoints |",
                  "|---|---:|---:|---:|"])
    for row in records:
        p = row["field_partition"]
        if p is not None and all(p["partition_vertices"]):
            lines.append(f"| {row['files'][0]['path']} | {p['partition_vertices'][0]} / {p['partition_vertices'][1]} | "
                         f"{p['cross_edges']} / {p['cross_edge_endpoints']} | "
                         f"{p['six_core_cross_edges']} / {p['six_core_cross_edge_endpoints']} |")
    lines.extend(["", "## Observed bottlenecks", "",
                  "The degree-ranked Haugland adaptive rounds raise degeneracy and retain their added vertices in the 6-core, but the radical-field cut remains at 60 edges on 63 endpoints. Of those edges, 58 share one endpoint, leaving only two cross edges after deleting that hub. Growth within the two sides has not strengthened their coupling.", "",
                  "The 7,753-vertex diagonal sqrt2 union has one articulation vertex in its 6-core. Its 60 field-crossing edges all meet one hub, so this rotation joins the old and rotated copies at their common center without adding a distributed coupling.", "",
                  "The unit-circle intersection candidate is different: almost all vertices survive peeling and its field-crossing edges touch almost every core vertex. Its new points have two old neighbors each, so their many new-to-new edges are essential to their survival. The sqrt2 translated-star graph also retains all vertices and distributes its field crossings widely. These are stronger coupling mechanisms by the audited criteria, without implying a chromatic lower bound.", "",
                  "The joint Haugland mixed-center experiment confirms the same effect: adding all 3,496 screened centers creates 5,627 vertices and 25,653 edges. Although every new center initially had two old neighbors, 1,296 survive joint 6-core peeling. The enlarged core has 3,427 vertices and 17,732 edges; its field-crossing edges outside the original hub increase from 2 to 1,298. A saved proper six-coloring still rejects this graph as a solution. The individual-center degree filter would have discarded this entire coupling mechanism.", "",
                  "## Search criterion", "",
                  "Recompute the 6-core after each proposed augmentation. Reject an augmentation that only adds vertices removed by peeling and adds no edge between old retained vertices. For circle-intersection batches, include interactions between new vertices before peeling: degree two to the old graph alone cannot constrain six colors. Rank surviving batches by growth in core cross edges after removing the largest boundary hub, number of distinct boundary vertices, and coverage of multiple independently obtained six-colorings.", "",
                  "A concrete next construction is joint unit-circle intersection closure on the existing Heule core: create centers at squared separations 4/3 and 8/3, add all exact new-to-old and new-to-new unit edges, peel to the 6-core, and retain only components supported by at least six internal neighbors. Then enumerate a fresh six-coloring and score subsequent intersection batches by core connections to its six color classes. Preserve every finite graph and every SAT coloring; no score is a non-colorability certificate.", ""])
    (root/"research"/"core_audit.md").write_text("\n".join(lines))
    print(json.dumps({"report": str(out), "graphs": len(records)}), flush=True)


if __name__ == "__main__":
    main()
