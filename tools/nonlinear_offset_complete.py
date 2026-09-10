#!/usr/bin/env python3
"""Bounded completion of remaining offset-pair separation witnesses."""

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import threading

from pysat.formula import CNF
from pysat.solvers import Solver

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import verify_tower as verifier

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/"runs/mixed-five-nonlinear-offset-screen/map-00"
OUT = ROOT/"runs/mixed-five-nonlinear-offset-completion"


def dump(path, value):
    Path(path).write_text(json.dumps(value, separators=(",", ":"))+"\n")


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    raw = (ROOT/"artifacts/mixed_haugland_g3_v31_tower.json").read_bytes()
    graph = json.loads(raw); sha = hashlib.sha256(raw).hexdigest()
    record = json.loads((SOURCE/"result.json").read_text())
    decomp = json.loads((ROOT/"runs/mixed-g3-v31-diagnostic/decompositions.json").read_text())
    factors = json.loads((ROOT/"runs/mixed-g3-v31-diagnostic/factor-coordinates.json").read_text())
    functional = record["linear_functional"]
    def lvalue(point):
        total = 0
        for coefficient, value in zip(functional, (value for part in point["parts"] for value in part)):
            rational = Fraction(value)
            if rational.denominator % 5 == 0:
                raise ValueError("coefficient denominator divisible by5")
            total += coefficient*rational.numerator*pow(rational.denominator, -1, 5)
        return total % 5
    avalues = [lvalue(point) for point in factors["g3"]]
    width = len(factors["v31"])
    def pair_clauses(pair):
        left, right = pair; a, b = decomp[left][0]; c, d = decomp[right][0]
        delta = (avalues[c]-avalues[a]) % 5
        return [[-(b*5+x+1), -(d*5+(x-delta) % 5+1)] for x in range(5)]
    formula = CNF(from_file=str(SOURCE/"restricted.cnf"))
    if verifier.common._sha256(SOURCE/"restricted.cnf") != record["restricted_cnf_sha256"]:
        raise ValueError("source restricted formula hash differs")
    if formula.clauses[-5:] != pair_clauses((8479, 8849)):
        raise ValueError("source formula lacks exact final conditional clauses")
    base = formula.clauses[:-5]
    pairs = json.loads((SOURCE/"remaining-unresolved-pairs.json").read_text())
    representatives = {}
    for pair in pairs:
        key = tuple(map(tuple, pair_clauses(pair)))
        representatives.setdefault(key, pair)
    goals = list(representatives.values())
    jobs = [("combined", goals)]+[(f"pair-{i:02d}", [pair]) for i, pair in enumerate(goals)]
    remaining = list(pairs); reports = []
    for name, selected_pairs in jobs:
        if not remaining:
            break
        output = OUT/name; output.mkdir()
        clauses = base+[clause for pair in selected_pairs for clause in pair_clauses(pair)]
        cnf = CNF(from_clauses=clauses); cnf.to_file(str(output/"restricted.cnf"))
        with Solver(name="g4", bootstrap_with=clauses, with_proof=True) as solver:
            solver.set_phases([b*5+c+1 if c == record["offset_colors"][b] else -(b*5+c+1)
                               for b in range(width) for c in range(5)])
            timer = threading.Timer(5, solver.interrupt); timer.start()
            try:
                answer = solver.solve_limited(expect_interrupt=True)
            finally:
                timer.cancel()
            model = set(solver.get_model()) if answer is True else None
            proof = solver.get_proof(); stats = solver.accum_stats()
        (output/"proof.drat").write_text("\n".join(proof)+("\n" if proof else ""))
        report = {"goal_satisfied": False, "graph_sha256": sha, "target_pairs": selected_pairs,
                  "restricted_cnf_sha256": verifier.common._sha256(output/"restricted.cnf"),
                  "proof_sha256": verifier.common._sha256(output/"proof.drat"), "statistics": stats,
                  "forced_equality_established": False}
        if answer is True:
            offsets = [next(c for c in range(5) if b*5+c+1 in model) for b in range(width)]
            colors = [(avalues[entries[0][0]]+offsets[entries[0][1]]) % 5 for entries in decomp]
            witness = {"k": 5, "colors": colors, "graph_sha256": sha}
            verifier.verify_coloring(graph, witness, 5)
            if any(colors[u] == colors[v] for u, v in selected_pairs):
                raise ValueError("full coloring does not separate a target pair")
            dump(output/"separating-coloring.json", witness)
            separated = [pair for pair in remaining if colors[pair[0]] != colors[pair[1]]]
            remaining = [pair for pair in remaining if colors[pair[0]] == colors[pair[1]]]
            dump(output/"additional-refuted-pairs.json", separated)
            report.update({"status": "SAT_FULL_COLORING_VERIFIED", "offset_colors": offsets,
                           "pairs_disproved_forced_equal": len(separated), "remaining_pairs": len(remaining)})
        else:
            report["status"] = "RESTRICTED_UNSAT_UNVERIFIED" if answer is False else "UNKNOWN_RESTRICTED"
        dump(output/"result.json", report); reports.append(report)
        print(json.dumps(report), flush=True)
    dump(OUT/"remaining-unresolved-pairs.json", remaining)
    dump(OUT/"result.json", {"goal_satisfied": False, "graph_sha256": sha, "initial_pairs": len(pairs),
         "offset_constraint_types": len(representatives), "remaining_pairs": len(remaining),
         "additional_pairs_disproved_forced_equal": len(pairs)-len(remaining),
         "models_verified": sum(report["status"] == "SAT_FULL_COLORING_VERIFIED" for report in reports),
         "forced_equality_established": False})


if __name__ == "__main__":
    main()
