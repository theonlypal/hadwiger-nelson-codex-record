#!/usr/bin/env python3
"""Finite proof that no distinct pair in one graph is forced equal at k=5.

This is a negative gadget search result, not a plane-coloring result.
Full-rank additive functionals separate every pair with different residues.
Every residual collision pair is covered by replayed compact Kempe
colorings or one of three independently checked full proper colorings.
No quadratic-size pair table is required.
"""

from collections import defaultdict
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import verify_tower as verifier
from tools.verify_kempe_witnesses import replay

ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT/"artifacts/mixed_haugland_g3_v31_tower.json"
FUNCTIONALS = ROOT/"runs/mixed-five-forced-pair-screen/functional-sampling.json"
KEMPE = ROOT/"runs/mixed-five-kempe-screen"
MODELS = [ROOT/"runs/mixed-five-seeded-pair-repair/25-over-3-repaired/separating-coloring.json",
          ROOT/"runs/mixed-five-nonlinear-offset-screen/map-00/separating-coloring.json",
          ROOT/"runs/mixed-five-nonlinear-offset-completion/combined/separating-coloring.json"]


def rank_mod5(rows):
    matrix = [[int(value) % 5 for value in row] for row in rows]
    pivot_row = 0
    for column in range(len(matrix[0])):
        chosen = next((i for i in range(pivot_row, len(matrix)) if matrix[i][column]), None)
        if chosen is None:
            continue
        matrix[pivot_row], matrix[chosen] = matrix[chosen], matrix[pivot_row]
        scale = pow(matrix[pivot_row][column], -1, 5)
        matrix[pivot_row] = [scale*x % 5 for x in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row != pivot_row:
                coefficient = matrix[row][column]
                matrix[row] = [(x-coefficient*y) % 5 for x, y in zip(matrix[row], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def verify():
    raw = GRAPH.read_bytes(); graph = json.loads(raw); sha = hashlib.sha256(raw).hexdigest()
    geometry = verifier.verify_geometry(graph)
    sampling = json.loads(FUNCTIONALS.read_text())
    if sampling["graph_sha256"] != sha:
        raise ValueError("functional source graph differs")
    cache = {}
    def residue(value):
        if value not in cache:
            rational = Fraction(value)
            if rational.denominator % 5 == 0:
                raise ValueError("coordinate denominator divisible by5")
            cache[value] = rational.numerator*pow(rational.denominator, -1, 5) % 5
        return cache[value]
    coordinates = np.array([[residue(value) for part in point["parts"] for value in part]
                            for point in graph["vertices"]], dtype=np.int16)
    active = np.flatnonzero(np.any(coordinates != 0, axis=0)).tolist()
    functionals = sampling["functionals"]
    if not functionals or any(len(row) != 48 or any(type(c) is not int or not 0 <= c < 5 for c in row) for row in functionals):
        raise ValueError("invalid functional list")
    rank = rank_mod5([[row[j] for j in active] for row in functionals])
    if rank != len(active):
        raise ValueError("functionals do not span all active coordinate columns")
    edges = np.asarray(graph["edges"], dtype=np.int64)
    for functional in functionals:
        colors = (coordinates@np.asarray(functional, dtype=np.int16)) % 5
        if np.any(colors[edges[:, 0]] == colors[edges[:, 1]]):
            raise ValueError("sample functional does not color the graph properly")
    classes = defaultdict(list)
    for vertex, residue_vector in enumerate(coordinates):
        classes[residue_vector.tobytes()].append(vertex)
    candidate_pairs = {pair for ids in classes.values() for pair in itertools.combinations(ids, 2)}
    refs = json.loads((KEMPE/"refutations.json").read_text())
    stored_pairs = [tuple(pair) for pair in refs["pairs"]]
    if len(stored_pairs) != len(candidate_pairs) or set(stored_pairs) != candidate_pairs:
        raise ValueError("Kempe candidate table is not exactly all residue-collision pairs")
    kempe = replay(GRAPH, KEMPE)
    if kempe["graph_sha256"] != sha:
        raise ValueError("Kempe replay graph differs")
    covered = {stored_pairs[int(key)] for key in refs["refutations_by_pair_index"]}
    remaining = candidate_pairs-covered
    additional = []
    for path in MODELS:
        witness = json.loads(path.read_text())
        if witness["graph_sha256"] != sha:
            raise ValueError("additional model graph differs")
        checked = verifier.verify_coloring(graph, witness, 5)
        colors = witness["colors"]
        separated = {pair for pair in remaining if colors[pair[0]] != colors[pair[1]]}
        remaining.difference_update(separated)
        additional.append({"path": str(path), "sha256": verifier.common._sha256(path),
                           "new_collision_pairs_separated": len(separated), "coloring_check": checked})
    if remaining:
        raise ValueError(f"{len(remaining)} residue-collision pairs remain unseparated")
    n = len(graph["vertices"])
    return {"goal_satisfied": False, "graph_sha256": sha, "vertices": geometry.vertices,
            "exact_unit_edges_checked": geometry.edges, "total_distinct_vertex_pairs": n*(n-1)//2,
            "active_residue_columns": active, "verified_functionals": len(functionals),
            "functional_rank_mod5": rank, "distinct_residue_classes": len(classes),
            "different_residue_pairs_separated_by_functionals": n*(n-1)//2-len(candidate_pairs),
            "residue_collision_pairs": len(candidate_pairs), "kempe_pairs_separated": len(covered),
            "additional_models": additional, "remaining_pairs": 0,
            "functional_source_sha256": verifier.common._sha256(FUNCTIONALS),
            "kempe_refutations_sha256": verifier.common._sha256(KEMPE/"refutations.json"),
            "conclusion": "For every two distinct vertices of this exact graph, a supplied proper five-coloring assigns different colors to them. No pair is forced equal in all five-colorings.",
            "scope": "This particular 62,992-vertex graph only; no claim about other graphs or all plane colorings"}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
