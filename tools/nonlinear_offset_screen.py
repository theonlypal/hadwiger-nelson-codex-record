#!/usr/bin/env python3
"""Small coloring ansatz: additive on G3, arbitrary labels on V31 offsets.

For a vertex with decomposition a+b, propose color L(a)+q_b modulo5.
All alternative exact decompositions and all graph edges constrain q.
Only an explicitly verified full graph coloring refutes forced equality;
UNSAT of these small formulas excludes just the specified ansatz and L.
"""

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import threading
import time

import numpy as np
from pysat.formula import CNF
from pysat.solvers import Solver

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import verify_tower as verifier

ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT/"artifacts/mixed_haugland_g3_v31_tower.json"
STRUCTURE = ROOT/"runs/mixed-g3-v31-diagnostic"
SOURCE = ROOT/"runs/mixed-five-forced-pair-screen/functional-sampling.json"
OUT = ROOT/"runs/mixed-five-nonlinear-offset-screen"
PAIR = (8479, 8849)


def dump(path, value):
    Path(path).write_text(json.dumps(value, separators=(",", ":"))+"\n")


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    raw = GRAPH.read_bytes(); graph = json.loads(raw); sha = hashlib.sha256(raw).hexdigest()
    decompositions = json.loads((STRUCTURE/"decompositions.json").read_text())
    factors = json.loads((STRUCTURE/"factor-coordinates.json").read_text())
    source = json.loads(SOURCE.read_text())
    if source["graph_sha256"] != sha or len(decompositions) != len(graph["vertices"]):
        raise ValueError("source graph binding differs")
    cache = {}
    def residue(value):
        if value not in cache:
            rational = Fraction(value)
            if rational.denominator % 5 == 0:
                raise ValueError("noninvertible coefficient denominator")
            cache[value] = rational.numerator*pow(rational.denominator, -1, 5) % 5
        return cache[value]
    def matrix(points):
        return np.asarray([[residue(value) for part in point["parts"] for value in part] for point in points], dtype=np.int16)
    coordinates = matrix(graph["vertices"]); acoords = matrix(factors["g3"]); bcoords = matrix(factors["v31"])
    first = np.asarray([entries[0] for entries in decompositions], dtype=np.int32)
    ai, bi = first[:, 0], first[:, 1]
    if np.any((acoords[ai]+bcoords[bi]-coordinates) % 5):
        raise ValueError("recorded factor decomposition differs modulo5")
    # Every proposed vertex color is also checked against all alternative
    # decompositions. Exact decomposition validity was separately audited.
    for vertex, entries in enumerate(decompositions):
        for a, b in entries:
            if np.any((acoords[a]+bcoords[b]-coordinates[vertex]) % 5):
                raise ValueError("alternative decomposition differs modulo5")
    edges = np.asarray(graph["edges"], dtype=np.int32)
    u, v = edges[:, 0], edges[:, 1]; width = len(bcoords)
    encoded_pair = (bi[u]*width+bi[v])*5
    reports = []; started = time.monotonic()
    checker = ROOT/"sources/drat-trim/drat-trim"
    config = {"goal_satisfied": False, "graph_sha256": sha, "pair": PAIR,
              "decompositions_sha256": verifier.common._sha256(STRUCTURE/"decompositions.json"),
              "factor_coordinates_sha256": verifier.common._sha256(STRUCTURE/"factor-coordinates.json"),
              "functional_file_sha256": verifier.common._sha256(SOURCE),
              "ansatz": "color(a+b)=L(a)+q_b modulo5; q_b arbitrary per V31 offset",
              "scope": "restricted coloring ansatz; only a full checked SAT witness refutes actual forced equality"}
    dump(OUT/"config.json", config)
    for number, functional in enumerate(source["functionals"]):
        output = OUT/f"map-{number:02d}"; output.mkdir()
        avalues = (acoords@np.asarray(functional, dtype=np.int16)) % 5
        constants = avalues[ai]
        forbidden = set()
        for code in np.unique(encoded_pair+(constants[v]-constants[u]) % 5).tolist():
            offset_pair, difference = divmod(code, 5); left, right = divmod(offset_pair, width)
            for color in range(5):
                forbidden.add(tuple(sorted((-(left*5+color+1), -(right*5+(color-difference) % 5+1)))))
        for vertex, entries in enumerate(decompositions):
            first_a, first_b = entries[0]
            for a, b in entries[1:]:
                delta = int(avalues[a]-avalues[first_a]) % 5
                for x in range(5):
                    for y in range(5):
                        if (x-y-delta) % 5:
                            forbidden.add(tuple(sorted((-(first_b*5+x+1), -(b*5+y+1)))))
        clauses = []
        for offset in range(width):
            clauses.append([offset*5+c+1 for c in range(5)])
            clauses.extend([-(offset*5+c+1), -(offset*5+d+1)] for c in range(5) for d in range(c+1, 5))
        clauses.extend(list(clause) for clause in sorted(forbidden))
        base_offsets = ((bcoords@np.asarray(functional, dtype=np.int16)) % 5).tolist()
        baseline_assignment = {b*5+c+1: c == base_offsets[b] for b in range(width) for c in range(5)}
        if any(not any(baseline_assignment[abs(lit)] == (lit > 0) for lit in clause) for clause in clauses):
            raise ValueError("known additive full coloring fails the ansatz encoding")
        left, right = PAIR; left_b, right_b = int(bi[left]), int(bi[right])
        delta = int(constants[right]-constants[left]) % 5
        pair_clauses = [[-(left_b*5+c+1), -(right_b*5+(c-delta) % 5+1)] for c in range(5)]
        clauses.extend(pair_clauses)
        formula = CNF(from_clauses=clauses); cnf = output/"restricted.cnf"; formula.to_file(str(cnf))
        with Solver(name="g4", bootstrap_with=clauses, with_proof=True) as solver:
            solver.set_phases([var if value else -var for var, value in baseline_assignment.items()])
            timer = threading.Timer(3, solver.interrupt); timer.start()
            try:
                answer = solver.solve_limited(expect_interrupt=True)
            finally:
                timer.cancel()
            model = set(solver.get_model()) if answer is True else None
            proof_lines = solver.get_proof(); stats = solver.accum_stats()
        proof = output/"proof.drat"; proof.write_text("\n".join(proof_lines)+("\n" if proof_lines else ""))
        report = {"map_index": number, "goal_satisfied": False, "linear_functional": functional,
                  "restricted_cnf_sha256": verifier.common._sha256(cnf), "proof_sha256": verifier.common._sha256(proof),
                  "statistics": stats, "forced_equality_established": False}
        if answer is True:
            offsets = [next(c for c in range(5) if b*5+c+1 in model) for b in range(width)]
            colors = ((constants+np.asarray(offsets, dtype=np.int16)[bi]) % 5).tolist()
            coloring = {"k": 5, "colors": colors, "graph_sha256": sha}
            verifier.verify_coloring(graph, coloring, 5)
            if colors[left] == colors[right]:
                raise ValueError("ansatz model fails desired pair separation")
            dump(output/"separating-coloring.json", coloring)
            report.update({"status": "SAT_FULL_COLORING_VERIFIED", "offset_colors": offsets,
                           "forced_equality_disproved": True, "pair_colors": [colors[left], colors[right]]})
        elif answer is False:
            check = subprocess.run([str(checker), str(cnf), str(proof)], capture_output=True, timeout=10)
            (output/"checker.log").write_bytes(check.stdout+check.stderr)
            verified = check.returncode == 0 and b"s VERIFIED" in check.stdout.splitlines()
            report.update({"status": "RESTRICTED_ANSATZ_UNSAT_VERIFIED" if verified else "RESTRICTED_ANSATZ_UNSAT_UNVERIFIED",
                           "checker_sha256": verifier.common._sha256(checker)})
        else:
            report["status"] = "UNKNOWN_RESTRICTED_ANSATZ"
        dump(output/"result.json", report); reports.append(report)
        print(json.dumps({"map_index": number, "status": report["status"]}), flush=True)
        if answer is True:
            break
    dump(OUT/"result.json", {**config, "maps_tested": len(reports), "elapsed_seconds": time.monotonic()-started,
         "statuses": [report["status"] for report in reports], "actual_forced_equality_established": False,
         "actual_forced_equality_disproved": any(report.get("forced_equality_disproved", False) for report in reports)})


if __name__ == "__main__":
    main()
