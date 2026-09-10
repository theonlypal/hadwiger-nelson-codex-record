#!/usr/bin/env python3
"""Checked Kempe-swap counterexamples to proposed forced-equal pairs.

A compact witness consists of a checked baseline five-coloring plus a list
of vertices whose two colors are exchanged. All affected edge incidences
are checked; all unaffected edges inherit the checked baseline. Remaining
pairs are unresolved, never asserted forced equal.
"""

from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import verify_tower as verifier

ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT/"artifacts/mixed_haugland_g3_v31_tower.json"
SOURCE = ROOT/"runs/mixed-five-forced-pair-screen"
OUT = ROOT/"runs/mixed-five-kempe-screen"
COLOR_PAIRS = [(a, b) for a in range(5) for b in range(a+1, 5)]
PAIR_INDEX = {pair: index for index, pair in enumerate(COLOR_PAIRS)}


def dump(path, value):
    Path(path).write_text(json.dumps(value, separators=(",", ":"))+"\n")


def find(parent, vertex):
    while parent[vertex] != vertex:
        parent[vertex] = parent[parent[vertex]]
        vertex = parent[vertex]
    return vertex


def bichromatic_components(graph, colors):
    n = len(colors); parents = [list(range(n)) for _ in COLOR_PAIRS]
    sizes = [[1]*n for _ in COLOR_PAIRS]
    for u, v in graph["edges"]:
        pair = tuple(sorted((colors[u], colors[v])))
        if pair[0] == pair[1]:
            raise ArithmeticError("baseline coloring is not proper")
        slot = PAIR_INDEX[pair]; parent = parents[slot]; size = sizes[slot]
        a, b = find(parent, u), find(parent, v)
        if a != b:
            if size[a] < size[b]:
                a, b = b, a
            parent[b] = a; size[a] += size[b]
    members = [defaultdict(list) for _ in COLOR_PAIRS]
    for vertex, color in enumerate(colors):
        for other in range(5):
            if color == other:
                continue
            slot = PAIR_INDEX[tuple(sorted((color, other)))]
            members[slot][find(parents[slot], vertex)].append(vertex)
    return parents, members


def checked_swap(colors, color_pair, component, adj):
    a, b = color_pair; changed = set(component); checked = 0
    def new_color(vertex):
        value = colors[vertex]
        return (b if value == a else a) if vertex in changed else value
    for vertex in component:
        if colors[vertex] not in color_pair:
            raise ArithmeticError("swap component contains another color")
        for other in adj[vertex]:
            if new_color(vertex) == new_color(other):
                raise ArithmeticError("Kempe swap produced a monochromatic edge")
            checked += 1
    return checked


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); deadline = started+150
    raw = GRAPH.read_bytes(); graph = json.loads(raw); graph_sha = hashlib.sha256(raw).hexdigest()
    sampling = json.loads((SOURCE/"functional-sampling.json").read_text())
    if graph_sha != sampling["graph_sha256"]:
        raise ValueError("source graph binding changed")
    groups = json.loads((SOURCE/"inseparable-classes.json").read_text())
    pairs = []
    for group in groups:
        if len(group) != 2:
            raise ValueError("this bounded screen expects the observed size-two collision classes")
        pairs.append(group)
    remaining = set(range(len(pairs))); refutations = {}; reports = []
    adj = [[] for _ in graph["vertices"]]
    for u, v in graph["edges"]:
        adj[u].append(v); adj[v].append(u)
    cache = {}
    def residue(value):
        if value not in cache:
            fraction = Fraction(value)
            if fraction.denominator % 5 == 0:
                raise ValueError("coefficient has denominator divisible by five")
            cache[value] = fraction.numerator*pow(fraction.denominator, -1, 5) % 5
        return cache[value]
    coords = np.array([[residue(value) for part in point["parts"] for value in part] for point in graph["vertices"]], dtype=np.int16)
    focus = tuple(sorted((6, 732)))
    model_number = 0
    for functional_number, functional in enumerate(sampling["functionals"]):
        colors = ((coords@np.array(functional, dtype=np.int16)) % 5).tolist()
        for walk_step in range(3 if functional_number < 3 else 1):
            if not remaining or time.monotonic() >= deadline:
                break
            model_out = OUT/f"model-{model_number:03d}"; model_out.mkdir()
            baseline = {"k": 5, "colors": colors, "graph_sha256": graph_sha}
            verifier.verify_coloring(graph, baseline, 5)
            dump(model_out/"baseline-coloring.json", baseline)
            # A previous deterministic swap can already separate a residue pair.
            separated_before = []
            for index in sorted(remaining):
                u, v = pairs[index]
                if colors[u] != colors[v]:
                    separated_before.append(index)
                    refutations[index] = {"kind": "checked-baseline-coloring", "witness": str((model_out/"baseline-coloring.json").relative_to(OUT))}
            remaining.difference_update(separated_before)
            parents, members = bichromatic_components(graph, colors)
            selected = defaultdict(list)
            for index in sorted(remaining):
                u, v = pairs[index]; color = colors[u]
                options = []
                for other in range(5):
                    if other == color:
                        continue
                    slot = PAIR_INDEX[tuple(sorted((color, other)))]; parent = parents[slot]
                    a, b = find(parent, u), find(parent, v)
                    if a != b:
                        root = a if len(members[slot][a]) <= len(members[slot][b]) else b
                        options.append((len(members[slot][root]), slot, root))
                if options:
                    _, slot, root = min(options)
                    selected[(slot, root)].append(index)
            checked_witnesses = []
            separated_by_swap = []
            for swap_index, ((slot, root), indices) in enumerate(sorted(selected.items())):
                component = members[slot][root]; color_pair = COLOR_PAIRS[slot]
                checked = checked_swap(colors, color_pair, component, adj)
                changed = set(component)
                for index in indices:
                    u, v = pairs[index]
                    if (u in changed) == (v in changed) or colors[u] != colors[v]:
                        raise ArithmeticError("claimed Kempe witness does not separate the pair")
                    refutations[index] = {"kind": "checked-compact-Kempe-coloring", "witness": str((model_out/"swap-witnesses.json").relative_to(OUT)), "swap_index": swap_index}
                    separated_by_swap.append(index)
                checked_witnesses.append({"color_pair": color_pair, "swapped_vertices": component,
                                          "affected_edge_incidences_checked": checked, "separated_pair_indices": indices})
                if any(tuple(sorted(pairs[index])) == focus for index in indices):
                    full = list(colors); a, b = color_pair
                    for vertex in component:
                        full[vertex] = b if full[vertex] == a else a
                    witness = {"k": 5, "colors": full, "graph_sha256": graph_sha}
                    verifier.verify_coloring(graph, witness, 5)
                    dump(OUT/"separating-6-732-coloring.json", witness)
            dump(model_out/"swap-witnesses.json", {"baseline": "baseline-coloring.json", "graph_sha256": graph_sha,
                                                   "witnesses": checked_witnesses})
            remaining.difference_update(separated_by_swap)
            report = {"model": model_number, "functional": functional_number, "walk_step": walk_step,
                      "component_counts": [len(group) for group in members],
                      "largest_components": [max(map(len, group.values()), default=0) for group in members],
                      "pairs_separated_by_baseline": len(separated_before), "pairs_separated_by_checked_swaps": len(separated_by_swap),
                      "compact_witnesses_checked": len(checked_witnesses), "pairs_remaining": len(remaining)}
            reports.append(report); dump(model_out/"result.json", report)
            print(json.dumps(report), flush=True); model_number += 1
            if not remaining:
                break
            # Deterministic Kempe walk: choose a proper subcomponent that
            # touches a remaining pair endpoint; this is a new proper model,
            # not a forced-equality inference.
            walk_options = []
            for index in sorted(remaining)[:128]:
                for vertex in pairs[index]:
                    for other in range(5):
                        if other == colors[vertex]:
                            continue
                        slot = PAIR_INDEX[tuple(sorted((colors[vertex], other)))]; root = find(parents[slot], vertex)
                        if len(members[slot]) > 1:
                            walk_options.append((len(members[slot][root]), slot, root))
            if not walk_options:
                break
            _, slot, root = min(walk_options); component = members[slot][root]
            checked_swap(colors, COLOR_PAIRS[slot], component, adj)
            a, b = COLOR_PAIRS[slot]; colors = list(colors)
            for vertex in component:
                colors[vertex] = b if colors[vertex] == a else a
        if not remaining or time.monotonic() >= deadline:
            break
    dump(OUT/"refutations.json", {"pairs": pairs, "refutations_by_pair_index": refutations})
    unresolved = [pairs[index] for index in sorted(remaining)]
    dump(OUT/"unresolved-pairs.json", unresolved)
    summary = {"goal_satisfied": False, "graph_sha256": graph_sha, "initial_candidate_pairs": len(pairs),
               "pairs_disproved_forced_equal": len(refutations), "unresolved_pairs": len(remaining),
               "models_checked": len(reports), "elapsed_seconds": time.monotonic()-started,
               "focus_pair_6_732_disproved": any(tuple(sorted(pairs[index])) == focus for index in refutations),
               "unresolved_interpretation": "Not separated by this bounded Kempe search; no forced-equality claim"}
    dump(OUT/"summary.json", summary); print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
