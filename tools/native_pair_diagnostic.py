#!/usr/bin/env python3
"""Phase-guided native conditional coloring diagnostics with an exact gate.

Only the graph coloring clauses and k endpoint-inequality clauses enter
the formula. A checked UNSAT proof establishes forced equality in every
k-coloring of this graph, never unrestricted non-k-colorability and never
the Hadwiger--Nelson target. Phase files contain search preferences only.
"""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hn.solve_external import parse_model, verify_full_cnf_model
from tools import verify_tower as verifier
from tools.seeded_pair_repair import archive_completed_proof

ROOT = Path(__file__).resolve().parents[1]


def dump(path, value):
    Path(path).write_text(json.dumps(value, separators=(",", ":"))+"\n")


def conditional_cnf(graph, pair, k=5):
    if type(k) is not int or k < 1:
        raise ValueError("color count must be a positive integer")
    if len(pair) != 2 or any(type(v) is not int or not 0 <= v < len(graph["vertices"]) for v in pair) or pair[0] == pair[1]:
        raise ValueError("invalid endpoint pair")
    base = verifier.coloring_cnf(graph, k)
    header, body = base.split(b"\n", 1); fields = header.decode().split()
    u, v = pair
    additions = [[-(u*k+c+1), -(v*k+c+1)] for c in range(k)]
    expected = (f"p cnf {fields[2]} {int(fields[3])+k}\n".encode()+body+
                b"".join((" ".join(map(str, clause))+" 0\n").encode() for clause in additions))
    return base, expected, additions


def check_conditional_proof(graph, pair, cnf, proof, checker, k=5, timeout=180):
    """Check exact augmented bytes, then DRAT on private immutable copies."""
    _, expected, additions = conditional_cnf(graph, pair, k)
    cnf, proof, checker = Path(cnf), Path(proof), Path(checker)
    if cnf.read_bytes() != expected:
        raise ValueError("conditional CNF differs from exact independent regeneration")
    proof_sha = verifier.common._sha256(proof)
    checker_sha = verifier.common._sha256(checker)
    with tempfile.TemporaryDirectory(prefix="hn-conditional-check-") as temporary:
        temporary = Path(temporary)
        sealed_cnf, sealed_proof = temporary/"formula.cnf", temporary/"proof.drat"
        sealed_cnf.write_bytes(expected); shutil.copyfile(proof, sealed_proof)
        if verifier.common._sha256(sealed_proof) != proof_sha:
            raise ValueError("proof copy differs from original")
        sealed_cnf.chmod(0o444); sealed_proof.chmod(0o444)
        result = subprocess.run([str(checker.resolve()), str(sealed_cnf), str(sealed_proof)],
                                capture_output=True, timeout=timeout)
        if sealed_cnf.read_bytes() != expected or verifier.common._sha256(sealed_proof) != proof_sha:
            raise ValueError("sealed proof inputs changed during checking")
    if cnf.read_bytes() != expected or verifier.common._sha256(proof) != proof_sha:
        raise ValueError("original proof inputs changed during checking")
    if verifier.common._sha256(checker) != checker_sha:
        raise ValueError("checker changed during checking")
    verified = result.returncode == 0 and b"s VERIFIED" in result.stdout.splitlines()
    return {"goal_satisfied": False, "conditional_unsat_verified": verified,
            "forced_equality_established": verified, "color_count": k,
            "pair": list(pair), "added_inequality_clauses": additions,
            "conditional_cnf_sha256": hashlib.sha256(expected).hexdigest(),
            "proof_sha256": proof_sha, "checker_sha256": checker_sha,
            "checker_exit_code": result.returncode,
            "scope": "forced pair equality in all k-colorings only"}, result.stdout, result.stderr


def run(graph_path, phase_coloring_path, out, pair=(8479, 8849), timeout=600):
    graph_path, phase_coloring_path = Path(graph_path).resolve(), Path(phase_coloring_path).resolve()
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=False)
    graph_raw = graph_path.read_bytes(); graph = json.loads(graph_raw)
    graph_sha = hashlib.sha256(graph_raw).hexdigest()
    phase_raw = phase_coloring_path.read_bytes(); phase_coloring = json.loads(phase_raw)
    if phase_coloring.get("graph_sha256") != graph_sha:
        raise ValueError("phase coloring graph hash differs")
    geometry = asdict(verifier.verify_geometry(graph))
    verifier.verify_coloring(graph, phase_coloring, 5)
    colors = phase_coloring["colors"]
    (out/"initial-phase-coloring.json").write_bytes(phase_raw)
    phase_path = out/"phases.txt"
    phase_path.write_text(" ".join(str(vertex*5+c+1 if colors[vertex] == c else -(vertex*5+c+1))
                                  for vertex in range(len(colors)) for c in range(5))+" 0\n")
    base, expected, additions = conditional_cnf(graph, pair, 5)
    cnf, proof, solution = out/"conditional.cnf", out/"proof.drat", out/"solution.out"
    cnf.write_bytes(expected)
    solver = ROOT/"sources/cadical/build/cadical-phase"
    source = ROOT/"tools/cadical_phase.cpp"
    checker = ROOT/"sources/drat-trim/drat-trim"
    command = [str(solver), "--sat", "--seed=0", "--no-colors", "-w", str(solution),
               str(cnf), str(proof), "--phase-file", str(phase_path)]
    config = {"goal_satisfied": False, "graph": str(graph_path), "graph_sha256": graph_sha,
              "geometry": geometry, "k": 5, "pair": list(pair), "timeout_seconds": timeout,
              "phase_source": str(phase_coloring_path), "phase_source_sha256": hashlib.sha256(phase_raw).hexdigest(),
              "phase_file_sha256": verifier.common._sha256(phase_path), "phase_coloring_verified": True,
              "phase_scope": "search preferences only; no fixed variables or assumptions",
              "base_cnf_sha256": hashlib.sha256(base).hexdigest(),
              "conditional_cnf_sha256": hashlib.sha256(expected).hexdigest(),
              "added_inequality_clauses": additions, "command": command,
              "solver_sha256": verifier.common._sha256(solver), "driver_source_sha256": verifier.common._sha256(source),
              "checker_sha256": verifier.common._sha256(checker), "proof_format": "binary DRAT",
              "scope": "conditional pair inequality; not an unrestricted chromatic lower bound"}
    dump(out/"config.json", config)
    started = time.monotonic(); timed_out = False
    with (out/"solver.stdout.log").open("wb") as stdout, (out/"solver.stderr.log").open("wb") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        dump(out/"process.json", {"pid": process.pid, "started_unix_seconds": time.time()})
        try:
            process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True; process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill(); process.wait()
    result = {**config, "native_exit_code": process.returncode, "timed_out": timed_out,
              "elapsed_seconds": time.monotonic()-started, "forced_equality_established": False,
              "proof_sha256": verifier.common._sha256(proof), "proof_bytes": proof.stat().st_size}
    dump(out/"native-outcome.json", result)
    if graph_path.read_bytes() != graph_raw or phase_coloring_path.read_bytes() != phase_raw:
        raise ValueError("source graph or source phase coloring changed during run")
    if (out/"initial-phase-coloring.json").read_bytes() != phase_raw or cnf.read_bytes() != expected:
        raise ValueError("sealed formula or phase coloring changed during run")
    for path, expected_sha in ((phase_path, config["phase_file_sha256"]), (solver, config["solver_sha256"]),
                               (source, config["driver_source_sha256"]), (checker, config["checker_sha256"])):
        if verifier.common._sha256(path) != expected_sha:
            raise ValueError(f"run input changed: {path}")
    if process.returncode == 10:
        assignment = parse_model(solution, len(colors)*5)
        clause_count = verify_full_cnf_model(cnf, assignment)
        colors = [next(c for c in range(5) if assignment[vertex*5+c+1]) for vertex in range(len(colors))]
        witness = {"k": 5, "colors": colors, "graph_sha256": graph_sha}
        verifier.verify_coloring(graph, witness, 5)
        if colors[pair[0]] == colors[pair[1]]:
            raise ValueError("native SAT model does not separate the endpoints")
        dump(out/"separating-coloring.json", witness)
        result.update({"status": "SAT_VERIFIED", "forced_equality_disproved": True,
                       "full_cnf_clauses_checked": clause_count, "pair_colors": [colors[v] for v in pair]})
    elif process.returncode == 20:
        result["status"] = "CONDITIONAL_UNSAT_UNVERIFIED"
        dump(out/"result.json", result)
        try:
            certificate, stdout, stderr = check_conditional_proof(graph, pair, cnf, proof, checker, k=5)
            (out/"checker.stdout.log").write_bytes(stdout); (out/"checker.stderr.log").write_bytes(stderr)
            dump(out/"conditional-proof-verification.json", certificate)
            result.update({"status": "CONDITIONAL_UNSAT_VERIFIED" if certificate["conditional_unsat_verified"] else "CONDITIONAL_UNSAT_UNVERIFIED",
                           "forced_equality_established": certificate["forced_equality_established"]})
        except subprocess.TimeoutExpired:
            result["checker_timed_out"] = True
    else:
        result["status"] = "UNKNOWN"
    dump(out/"result.json", result)
    if result["status"] in ("SAT_VERIFIED", "UNKNOWN"):
        result = archive_completed_proof(out)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("graph", type=Path)
    parser.add_argument("phase_coloring", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pair", type=int, nargs=2, default=[8479, 8849])
    parser.add_argument("--timeout", type=float, default=600)
    args = parser.parse_args()
    print(json.dumps(run(args.graph, args.phase_coloring, args.out, tuple(args.pair), args.timeout), indent=2))
