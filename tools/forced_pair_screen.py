#!/usr/bin/env python3
"""Bounded screening for possible forced-equal pairs in five-colorings.

Different colors in any verified coloring refute forced equality. Equality
in every sampled additive coloring is only a candidate filter. Conditional
native tests add exactly the five clauses demanding different endpoint
colors; only an independently verified DRAT refutation can establish the
conditional forced-equality claim. No outcome here satisfies the non-6 goal.
"""

from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import itertools
import json
from math import isqrt
from pathlib import Path
import random
import subprocess
import sys
import threading
import time

import numpy as np
from flint import fmpq, fmpq_poly
from pysat.formula import CNF
from pysat.solvers import Solver

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import verify_tower as verifier
from hn.solve_external import parse_model, verify_full_cnf_model, compress_completed_run


ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT/"artifacts/mixed_haugland_g3_v31_tower.json"
LINEAR = ROOT/"runs/mixed-haugland-g3-v31-linear5"
OUT = ROOT/"runs/mixed-five-forced-pair-screen"


def dump(path, value):
    Path(path).write_text(json.dumps(value, separators=(",", ":"))+"\n")


def row_basis(rows, width):
    basis = {}
    for original in rows:
        row = [int(value) % 5 for value in original]
        for pivot in sorted(basis):
            factor = row[pivot]
            if factor:
                row = [(x-factor*y) % 5 for x, y in zip(row, basis[pivot])]
        pivot = next((j for j, value in enumerate(row) if value), None)
        if pivot is not None:
            inverse = pow(row[pivot], -1, 5)
            basis[pivot] = [(inverse*x) % 5 for x in row]
    return [basis[pivot] for pivot in sorted(basis)]


def norm(parts, arithmetic):
    conjugates = [arithmetic.conjugate(part) for part in parts]
    output = [fmpq_poly([]) for _ in range(4)]
    for i in range(4):
        for j in range(4):
            output[i ^ j] += verifier.SQUARES[i & j]*parts[i]*conjugates[j]
    return tuple(part % verifier.PHI42 for part in output)


def multiply(left, right):
    output = [fmpq_poly([]) for _ in range(4)]
    for i in range(4):
        for j in range(4):
            output[i ^ j] += verifier.SQUARES[i & j]*left[i]*right[j]
    return tuple(part % verifier.PHI42 for part in output)


def encoded(parts):
    return {"parts": [[str(part[j]) for j in range(12)] for part in parts]}


def spindle_rotation(distance_squared, delta, arithmetic):
    if any(distance_squared[i] for i in range(1, 4)) or distance_squared[0].degree() > 0:
        return {"status": "nonrational-distance-not-analyzed"}
    s = Fraction(str(distance_squared[0][0]))
    if s <= 0:
        raise ArithmeticError("distinct points have nonpositive rational squared distance")
    if s < Fraction(1, 4):
        return {"status": "too-close-for-unit-spindle", "squared_distance_rational": str(s)}
    t = fmpq_poly([0, 1]); b = (2*t**7-1) % verifier.PHI42
    z7 = t**6 % verifier.PHI42
    gauss7 = sum((z7**j for j in (1, 2, 4)), fmpq_poly([]))-sum((z7**j for j in (3, 5, 6)), fmpq_poly([]))
    root21 = (-b*gauss7) % verifier.PHI42
    if root21*root21 % verifier.PHI42 != 21:
        raise ArithmeticError("independent sqrt21 identity failed")
    root_options = {
        1: (fmpq_poly([1]), fmpq_poly([]), fmpq_poly([]), fmpq_poly([])),
        5: (fmpq_poly([]), fmpq_poly([1]), fmpq_poly([]), fmpq_poly([])),
        21: (root21, fmpq_poly([]), fmpq_poly([]), fmpq_poly([])),
        33: (fmpq_poly([]), fmpq_poly([]), fmpq_poly([1]), fmpq_poly([])),
        77: (fmpq_poly([]), fmpq_poly([]), root21/3, fmpq_poly([])),
        105: (fmpq_poly([]), root21, fmpq_poly([]), fmpq_poly([])),
        165: (fmpq_poly([]), fmpq_poly([]), fmpq_poly([]), fmpq_poly([1])),
        385: (fmpq_poly([]), fmpq_poly([]), fmpq_poly([]), root21/3),
    }
    radicand = 3*(4*s-1)
    chosen = None
    if radicand == 0:
        chosen = (1, tuple(fmpq_poly([]) for _ in range(4)))
    else:
        for squareclass, root in root_options.items():
            ratio = radicand/squareclass
            numerator, denominator = isqrt(ratio.numerator), isqrt(ratio.denominator)
            if numerator*numerator == ratio.numerator and denominator*denominator == ratio.denominator:
                scale = fmpq(numerator, denominator)
                chosen = squareclass, tuple(scale*part for part in root)
                break
    report = {"squared_distance_rational": str(s), "required_real_radicand": str(radicand),
              "rotation_cosine": str(1-Fraction(1, 2)/s)}
    if chosen is None:
        report["status"] = "no-root-in-listed-known-quadratic-subfields"
        return report
    squareclass, root = chosen
    denominator = fmpq(str(6*s))
    rotation = [b*part/denominator % verifier.PHI42 for part in root]
    rotation[0] += fmpq(str(1-Fraction(1, 2)/s))
    difference = list(rotation); difference[0] -= 1
    if not arithmetic.is_unit(rotation) or not arithmetic.is_unit(multiply(difference, delta)):
        raise ArithmeticError("spindle rotation failed independent exact identities")
    report.update({"status": "exact-rotation-fits-current-tower", "known_squareclass": squareclass,
                   "rotation": encoded(rotation), "rotation_unit_verified": True, "rotated_endpoint_gap_unit_verified": True})
    return report


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); deadline = started+180
    raw = GRAPH.read_bytes(); graph = json.loads(raw); graph_sha = hashlib.sha256(raw).hexdigest()
    baseline = json.loads((LINEAR/"result.json").read_text())
    if graph_sha != baseline["graph_sha256"]:
        raise ValueError("linear screen graph binding mismatch")
    residues = {}
    def residue(value):
        if value not in residues:
            rational = Fraction(value)
            if rational.denominator % 5 == 0:
                raise ValueError("nonintegral coefficient at 5")
            residues[value] = rational.numerator*pow(rational.denominator, -1, 5) % 5
        return residues[value]
    coords = np.array([[residue(value) for part in point["parts"] for value in part] for point in graph["vertices"]], dtype=np.int16)
    edges = np.asarray(graph["edges"], dtype=np.int64)
    active = np.flatnonzero(np.any(coords != 0, axis=0)).tolist()
    directions = np.array(baseline["directions"], dtype=np.int16)
    formula = CNF(from_file=str(LINEAR/"linear.cnf"))
    rng = random.Random(73013); functionals = []; attempts = []
    def accept(coefficients):
        values = np.asarray(coefficients, dtype=np.int16)
        if np.any((directions@values) % 5 == 0):
            raise ValueError("sampled functional violates a direction constraint")
        colors = (coords@values) % 5
        if np.any(colors[edges[:, 0]] == colors[edges[:, 1]]):
            raise ValueError("sampled functional is not a proper graph coloring")
        if coefficients not in functionals:
            functionals.append(coefficients)
    accept(baseline["linear_functional"])
    sample_deadline = min(deadline-75, time.monotonic()+40)
    with Solver(name="g4", bootstrap_with=formula.clauses) as solver:
        for attempt in range(96):
            remaining = sample_deadline-time.monotonic()
            if remaining <= 0:
                break
            if attempt % 3 == 0:
                assumptions = []
            else:
                indices = rng.sample([j for j in active if j != 0], min(2, len(active)-1))
                assumptions = [j*5+rng.randrange(5)+1 for j in indices]
            phases = [j*5+rng.randrange(5)+1 for j in active]
            solver.set_phases(phases)
            timer = threading.Timer(min(remaining, 0.4), solver.interrupt); timer.start()
            try:
                answer = solver.solve_limited(assumptions=assumptions, expect_interrupt=True)
            finally:
                timer.cancel(); solver.clear_interrupt()
            attempts.append({"assumptions": assumptions, "status": "SAT" if answer is True else "UNSAT-sampling-assumptions" if answer is False else "UNKNOWN"})
            if answer is True:
                model = set(solver.get_model())
                coefficients = [next(c for c in range(5) if j*5+c+1 in model) for j in range(48)]
                accept(coefficients)
                solver.add_clause([-(j*5+coefficients[j]+1) for j in active])
                if len(row_basis([[row[j] for j in active] for row in functionals], len(active))) == len(active):
                    break
    basis = row_basis(functionals, 48)
    signatures = (coords@np.array(basis, dtype=np.int16).T) % 5
    signature_classes = defaultdict(list); residue_classes = defaultdict(list)
    for vertex in range(len(coords)):
        signature_classes[signatures[vertex].tobytes()].append(vertex)
        residue_classes[coords[vertex].tobytes()].append(vertex)
    nontrivial = [ids for ids in signature_classes.values() if len(ids) > 1]
    residue_nontrivial = [ids for ids in residue_classes.values() if len(ids) > 1]
    sample_report = {
        "goal_satisfied": False, "graph_sha256": graph_sha, "functionals": functionals,
        "sampled_functionals_verified": len(functionals), "active_coordinate_columns": active,
        "functional_span_rank": len(basis), "active_functional_span_rank": len(row_basis([[row[j] for j in active] for row in functionals], len(active))),
        "attempts": attempts, "basis_mod5": basis,
        "signature_classes": len(signature_classes), "nontrivial_signature_classes": len(nontrivial),
        "pairs_not_separated_by_samples": sum(len(ids)*(len(ids)-1)//2 for ids in nontrivial),
        "identical_residue_classes": len(residue_classes), "nontrivial_identical_residue_classes": len(residue_nontrivial),
        "identical_residue_pairs": sum(len(ids)*(len(ids)-1)//2 for ids in residue_nontrivial),
        "signature_class_size_histogram": dict(Counter(map(len, signature_classes.values()))),
        "interpretation": "Every separated pair is not forced equal; equality in samples never establishes forced equality",
    }
    dump(OUT/"functional-sampling.json", sample_report)
    dump(OUT/"inseparable-classes.json", nontrivial)
    print(json.dumps({"stage": "functionals_sampled", "samples": len(functionals), "rank": len(basis), "active_columns": len(active), "remaining_pairs": sample_report["pairs_not_separated_by_samples"], "identical_residue_pairs": sample_report["identical_residue_pairs"]}), flush=True)
    arithmetic = verifier.TowerArithmetic(); parsed = {}
    def point(index):
        if index not in parsed:
            parsed[index] = arithmetic.parse(graph["vertices"][index], index)
        return parsed[index]
    candidates = []
    for group in sorted(nontrivial, key=lambda ids: (-len(ids), ids[0])):
        for u, v in itertools.combinations(group[:8], 2):
            delta = tuple(a-b for a, b in zip(point(v), point(u)))
            distance_squared = norm(delta, arithmetic)
            rotation = spindle_rotation(distance_squared, delta, arithmetic)
            candidates.append({"vertices": [u, v], "identical_residues": bool(np.array_equal(coords[u], coords[v])),
                               "squared_distance": encoded(distance_squared), "spindle": rotation,
                               "forced_equality_established": False})
            if len(candidates) >= 160:
                break
        if len(candidates) >= 160:
            break
    candidates.sort(key=lambda p: (p["spindle"]["status"] != "exact-rotation-fits-current-tower", p["spindle"]["status"] == "nonrational-distance-not-analyzed", p["vertices"]))
    dump(OUT/"candidate-pairs.json", candidates)
    results = []
    if candidates and deadline-time.monotonic() > 12:
        base_cnf = verifier.coloring_cnf(graph, 5)
        header, body = base_cnf.split(b"\n", 1); parts = header.decode().split()
        for index, candidate in enumerate(candidates[:3]):
            remaining = deadline-time.monotonic()
            if remaining < 12:
                break
            limit = min(15, int(remaining-8)); run_out = OUT/f"conditional-{index:02d}"; run_out.mkdir()
            u, v = candidate["vertices"]
            additions = [[-(u*5+c+1), -(v*5+c+1)] for c in range(5)]
            cnf = (f"p cnf {parts[2]} {int(parts[3])+5}\n".encode()+body+
                   b"".join((" ".join(map(str, clause))+" 0\n").encode() for clause in additions))
            cnf_path = run_out/"conditional.cnf"; cnf_path.write_bytes(cnf)
            proof = run_out/"proof.drat"; solution = run_out/"solution.out"
            solver = ROOT/"sources/cadical/build/cadical"; checker = ROOT/"sources/drat-trim/drat-trim"
            command = [str(solver), "--sat", "--no-colors", "-t", str(limit), "-w", str(solution), str(cnf_path), str(proof)]
            run_started = time.monotonic()
            with (run_out/"solver.stdout.log").open("wb") as stdout, (run_out/"solver.stderr.log").open("wb") as stderr:
                process = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=limit+5)
            result = {"goal_satisfied": False, "purpose": "conditional five-color pair inequality test only",
                      "vertices": [u, v], "added_inequality_clauses": additions, "base_graph_sha256": graph_sha,
                      "base_five_cnf_sha256": hashlib.sha256(base_cnf).hexdigest(), "conditional_cnf_sha256": hashlib.sha256(cnf).hexdigest(),
                      "native_exit_code": process.returncode, "elapsed_seconds": time.monotonic()-run_started,
                      "command": command, "solver_sha256": verifier.common._sha256(solver),
                      "proof_sha256": verifier.common._sha256(proof), "proof_bytes": proof.stat().st_size,
                      "forced_equality_established": False}
            if process.returncode == 10:
                assignment = parse_model(solution, 5*len(graph["vertices"]))
                verify_full_cnf_model(cnf_path, assignment)
                colors = [next(c for c in range(5) if assignment[vertex*5+c+1]) for vertex in range(len(coords))]
                coloring = {"k": 5, "colors": colors, "graph_sha256": graph_sha}
                verifier.verify_coloring(graph, coloring, 5)
                if colors[u] == colors[v]:
                    raise ArithmeticError("conditional SAT witness does not separate its pair")
                dump(run_out/"separating-coloring.json", coloring)
                result.update({"status": "SAT", "forced_equality_disproved": True})
            elif process.returncode == 20:
                result["status"] = "CONDITIONAL_UNSAT_UNVERIFIED"
                if cnf_path.read_bytes() != cnf:
                    raise ArithmeticError("conditional CNF binding changed")
                check_limit = max(1, min(45, deadline-time.monotonic()))
                try:
                    checked = subprocess.run([str(checker), str(cnf_path), str(proof)], capture_output=True, text=True, timeout=check_limit)
                    (run_out/"checker.log").write_text(checked.stdout+checked.stderr)
                    if checked.returncode == 0 and any(line.strip() == "s VERIFIED" for line in checked.stdout.splitlines()):
                        result.update({"status": "CONDITIONAL_UNSAT_VERIFIED", "forced_equality_established": True,
                                       "checker_sha256": verifier.common._sha256(checker)})
                except subprocess.TimeoutExpired:
                    result["checker_timeout"] = True
            else:
                result["status"] = "UNKNOWN_NATIVE"
            dump(run_out/"result.json", result)
            if result["status"] in ("SAT", "UNKNOWN_NATIVE"):
                result = compress_completed_run(run_out)
            results.append(result)
            print(json.dumps({"stage": "conditional_test", "pair": [u, v], "status": result["status"]}), flush=True)
    summary = {"goal_satisfied": False, "graph_sha256": graph_sha, "sampling": {key: value for key, value in sample_report.items() if key not in ("functionals", "attempts", "basis_mod5")},
               "candidate_pairs_analyzed": len(candidates), "rotation_fits_current_tower": sum(p["spindle"]["status"] == "exact-rotation-fits-current-tower" for p in candidates),
               "conditional_tests": [{"vertices": result["vertices"], "status": result["status"], "forced_equality_established": result["forced_equality_established"]} for result in results],
               "elapsed_seconds": time.monotonic()-started,
               "scope": "Unseparated sampled pairs remain candidates, never forced-equality claims without checked conditional UNSAT"}
    dump(OUT/"summary.json", summary); print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
