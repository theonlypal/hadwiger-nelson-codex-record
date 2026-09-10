#!/usr/bin/env python3
"""Check finite obstructions to a restricted sixth-color repair ansatz.

Given an additive tentative five-color map L, the ansatz selects an
independent vertex set S meeting every L-monochromatic edge and recolors
S with color 6. Each edge gives (!s_u OR !s_v); monochromatic edges also
give (s_u OR s_v). Verified odd cycles or contradictory implication paths
refute ONLY this ansatz for the specified L, never arbitrary 6-colorings.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.verify import ExactField


def verify(directory):
    directory = Path(directory)
    record = json.loads((directory/"result.json").read_text())
    graph_raw = Path(record["graph"]).read_bytes(); graph = json.loads(graph_raw)
    sha = hashlib.sha256(graph_raw).hexdigest()
    if sha != record["graph_sha256"]:
        raise ValueError("graph hash mismatch")
    if (graph["schema"] != "hn-exact-v1" or graph["field"].get("cyclotomic_order") != 42
            or graph["field"].get("radicand") != 5):
        raise ValueError("this verifier requires the specified cyclotomic sqrt5 basis")
    edges = {tuple(sorted(edge)) for edge in graph["edges"]}
    field = ExactField(graph["field"]); parsed = {}; residue = {}; unit_checked = set()
    def point(vertex):
        if type(vertex) is not int or not 0 <= vertex < len(graph["vertices"]):
            raise ValueError("out of range proof vertex")
        if vertex not in parsed:
            parsed[vertex] = field.parse(graph["vertices"][vertex], vertex)
        return parsed[vertex]
    def coefficients(vertex):
        point(vertex)
        if vertex not in residue:
            values = []
            for value in graph["vertices"][vertex]["a"]+graph["vertices"][vertex]["b"]:
                rational = Fraction(value)
                if rational.denominator % 5 == 0:
                    raise ValueError("coefficient denominator is not invertible modulo5")
                values.append(rational.numerator*pow(rational.denominator, -1, 5) % 5)
            residue[vertex] = values
        return residue[vertex]
    def unit_edge(u, v):
        edge = tuple(sorted((u, v)))
        if edge not in edges:
            raise ValueError("proof invokes an absent edge")
        if edge not in unit_checked:
            if not field.unit_difference(point(u), point(v)):
                raise ValueError("proof edge is not exactly unit length")
            unit_checked.add(edge)
    reports = []
    for path in sorted(directory.glob("map-*.json")):
        proof = json.loads(path.read_text()); functional = proof["linear_functional"]
        if len(functional) != 24 or any(type(c) is not int or not 0 <= c < 5 for c in functional):
            raise ValueError("invalid F5 functional")
        def color(vertex):
            return sum(x*y for x, y in zip(coefficients(vertex), functional)) % 5
        if "odd_cycle" in proof:
            cycle = proof["odd_cycle"]; length = len(cycle)-1
            if length < 3 or length % 2 != 1 or cycle[0] != cycle[-1] or len(set(cycle[:-1])) != length:
                raise ValueError("not a simple closed odd cycle")
            for u, v in zip(cycle, cycle[1:]):
                unit_edge(u, v)
                if color(u) != color(v):
                    raise ValueError("odd-cycle edge is not L-monochromatic")
            reports.append({"map_index": proof["map_index"], "kind": "monochromatic-odd-cycle",
                            "cycle_length": length, "verified": True})
        elif "implication_paths" in proof:
            paths = proof["implication_paths"]
            if len(paths) != 2 or any(len(chain) < 2 for chain in paths):
                raise ValueError("expected two implication paths")
            literal = proof["contradiction_vertex"]+1
            if [paths[0][0], paths[0][-1], paths[1][0], paths[1][-1]] != [literal, -literal, -literal, literal]:
                raise ValueError("paths do not establish both directions of a contradiction")
            for chain in paths:
                for left, right in zip(chain, chain[1:]):
                    if type(left) is not int or type(right) is not int or left*right >= 0:
                        raise ValueError("not an implication of an ansatz edge clause")
                    u, v = abs(left)-1, abs(right)-1; unit_edge(u, v)
                    if left < 0 and color(u) != color(v):
                        raise ValueError("positive cover clause requires a monochromatic edge")
            reports.append({"map_index": proof["map_index"], "kind": "contradictory-implication-paths",
                            "path_lengths": [len(chain)-1 for chain in paths], "verified": True})
        else:
            raise ValueError("no independently checkable finite obstruction supplied")
    return {"goal_satisfied": False, "graph_sha256": sha, "maps_checked": len(reports),
            "unique_exact_unit_edges_checked": len(unit_checked), "coordinate_vertices_checked": len(parsed),
            "obstructions": reports, "scope": "Only the specified additive map plus one independent sixth-color class is refuted",
            "arbitrary_six_colorability_status": "not determined by these witnesses"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("directory", type=Path)
    args = parser.parse_args(); print(json.dumps(verify(args.directory), indent=2))
