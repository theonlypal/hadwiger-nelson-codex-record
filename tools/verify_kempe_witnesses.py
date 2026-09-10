#!/usr/bin/env python3
"""Replay compact five-coloring witnesses without trusting Kempe discovery.

Every baseline is checked on every graph edge. For a compact swap, every
edge incident to a changed vertex is checked under the reconstructed
coloring; other edges retain their checked baseline colors. This proves
the full reconstructed coloring proper without materializing it. Pair
separation and the disjoint partition of refuted/unresolved pairs are
also checked. This does not claim any unresolved pair is forced equal.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import verify_tower as verifier


def replay(graph_path, witness_directory):
    graph_path, out = Path(graph_path), Path(witness_directory)
    raw = graph_path.read_bytes(); graph = json.loads(raw)
    graph_sha = hashlib.sha256(raw).hexdigest()
    refs = json.loads((out/"refutations.json").read_text())
    pairs = refs["pairs"]; entries = refs["refutations_by_pair_index"]
    unresolved = json.loads((out/"unresolved-pairs.json").read_text())
    covered = {}; adjacency = [[] for _ in graph["vertices"]]
    for u, v in graph["edges"]:
        adjacency[u].append(v); adjacency[v].append(u)
    n = len(adjacency)
    if any(len(pair) != 2 or any(type(v) is not int or not 0 <= v < n for v in pair)
           or pair[0] == pair[1] for pair in pairs):
        raise ValueError("invalid candidate pair")
    if len({tuple(sorted(pair)) for pair in pairs}) != len(pairs):
        raise ValueError("duplicate candidate pair")
    model_count = swap_count = incidence_count = 0
    for model in sorted(out.glob("model-*")):
        baseline = json.loads((model/"baseline-coloring.json").read_text())
        if baseline.get("graph_sha256") != graph_sha:
            raise ValueError("baseline graph hash mismatch")
        verifier.verify_coloring(graph, baseline, 5)
        colors = baseline["colors"]; model_count += 1
        baseline_path = str((model/"baseline-coloring.json").relative_to(out))
        for key, ref in entries.items():
            if ref == {"kind": "checked-baseline-coloring", "witness": baseline_path}:
                index = int(key); u, v = pairs[index]
                if colors[u] == colors[v] or index in covered:
                    raise ValueError("baseline does not uniquely establish its claimed separation")
                covered[index] = ref
        record = json.loads((model/"swap-witnesses.json").read_text())
        if record["graph_sha256"] != graph_sha or record["baseline"] != "baseline-coloring.json":
            raise ValueError("swap graph/baseline binding mismatch")
        for index, witness in enumerate(record["witnesses"]):
            a, b = witness["color_pair"]
            if type(a) is not int or type(b) is not int or not 0 <= a < b < 5:
                raise ValueError("invalid swap color pair")
            changed = set(witness["swapped_vertices"])
            if len(changed) != len(witness["swapped_vertices"]) or any(
                    type(u) is not int or not 0 <= u < n or colors[u] not in (a, b) for u in changed):
                raise ValueError("invalid swap vertex list")
            def color(u):
                return (b if colors[u] == a else a) if u in changed else colors[u]
            local_incidence_count = 0
            for u in changed:
                for v in adjacency[u]:
                    if color(u) == color(v):
                        raise ValueError("compact coloring has a monochromatic edge")
                    local_incidence_count += 1
            if local_incidence_count != witness["affected_edge_incidences_checked"]:
                raise ValueError("reported incident edge count differs")
            incidence_count += local_incidence_count
            for key in witness["separated_pair_indices"]:
                if type(key) is not int or not 0 <= key < len(pairs) or key in covered:
                    raise ValueError("invalid or duplicate separated pair index")
                u, v = pairs[key]
                if color(u) == color(v):
                    raise ValueError("compact coloring does not separate its claimed pair")
                covered[key] = {"kind": "checked-compact-Kempe-coloring",
                                "witness": str((model/"swap-witnesses.json").relative_to(out)),
                                "swap_index": index}
            swap_count += 1
    if {str(key): value for key, value in covered.items()} != entries:
        raise ValueError("refutation index differs from replayed evidence")
    expected_unresolved = [pair for index, pair in enumerate(pairs) if index not in covered]
    if unresolved != expected_unresolved:
        raise ValueError("unresolved pairs are not the exact complement")
    return {"goal_satisfied": False, "graph_sha256": graph_sha,
            "models_checked": model_count, "compact_swaps_checked": swap_count,
            "affected_edge_incidences_checked": incidence_count,
            "pairs_disproved_forced_equal": len(covered), "unresolved_pairs": len(unresolved),
            "scope": "coloring and pair-separation witnesses only; geometry not rerun",
            "unresolved_pairs_have_no_forced_equality_claim": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", type=Path)
    parser.add_argument("witness_directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(replay(args.graph, args.witness_directory), indent=2))
