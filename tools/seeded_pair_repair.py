#!/usr/bin/env python3
"""One bounded, fully witnessed conditional five-coloring repair.

The CNF is unrestricted graph 5-coloring plus five pair-inequality clauses.
The initial phases are a known full proper coloring selected to minimize
the conflict cost of recoloring an endpoint. Neither phases nor UNKNOWN imply
anything about forced equality. A checked conditional UNSAT would only
establish the pair gadget claim, never the non-6 target.
"""

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
import sys
import threading
import time

from pysat.formula import CNF
from pysat.solvers import Solver

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hn.solve_external import verify_full_cnf_model
from tools import verify_tower as verifier

ROOT = Path(__file__).resolve().parents[1]
GRAPH = ROOT/"artifacts/mixed_haugland_g3_v31_tower.json"
KEMPE = ROOT/"runs/mixed-five-kempe-screen"
OUT = ROOT/"runs/mixed-five-seeded-pair-repair"
PAIR = (979, 2711)


def dump(path, value):
    Path(path).write_text(json.dumps(value, separators=(",", ":"))+"\n")


def archive_completed_proof(out):
    """Keep the exact completed proof bytes in a checked lossless archive."""
    out = Path(out); result = json.loads((out/"result.json").read_text())
    if result["status"] not in ("SAT_VERIFIED", "UNKNOWN"):
        raise ValueError("only completed diagnostic outcomes can be archived")
    if "proof_storage" in result or not result["proof_bytes"]:
        return result
    proof = out/"proof.drat"; archive = out/"proof.drat.gz"
    expected_sha = result["proof_sha256"]
    if verifier.common._sha256(proof) != expected_sha:
        raise ValueError("proof changed before archival")
    with archive.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=6, fileobj=raw, mtime=0) as compressed:
            with proof.open("rb") as original:
                shutil.copyfileobj(original, compressed, 1024*1024)
    recovered = hashlib.sha256(); size = 0
    with gzip.open(archive, "rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            recovered.update(block); size += len(block)
    if recovered.hexdigest() != expected_sha or size != result["proof_bytes"]:
        raise ValueError("archive did not reproduce exact proof bytes")
    if verifier.common._sha256(proof) != expected_sha:
        raise ValueError("proof changed during archival")
    result["proof_storage"] = {"compression": "gzip", "path": str(archive),
        "archive_sha256": verifier.common._sha256(archive), "archive_bytes": archive.stat().st_size,
        "uncompressed_sha256": expected_sha, "uncompressed_bytes": size,
        "decompression_verified": True}
    dump(out/"result.json", result)
    proof.unlink()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair", type=int, nargs=2, default=PAIR)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--timeout", type=float, default=30)
    args = parser.parse_args()
    out = args.out.resolve(); pair = tuple(args.pair)
    out.mkdir(parents=True, exist_ok=False)
    raw = GRAPH.read_bytes(); graph = json.loads(raw)
    graph_sha = hashlib.sha256(raw).hexdigest()
    adjacency = [[] for _ in graph["vertices"]]
    for u, v in graph["edges"]:
        adjacency[u].append(v); adjacency[v].append(u)
    seeds = []
    for path in sorted(KEMPE.glob("model-*/baseline-coloring.json")):
        model = json.loads(path.read_text()); verifier.verify_coloring(graph, model, 5)
        if model.get("graph_sha256") != graph_sha:
            raise ValueError("seed graph hash mismatch")
        for endpoint in pair:
            counts = Counter(model["colors"][v] for v in adjacency[endpoint])
            for color in range(5):
                if color != model["colors"][endpoint]:
                    seeds.append((counts[color], str(path), endpoint, color))
    conflicts, seed_path, endpoint, color = min(seeds)
    baseline = json.loads(Path(seed_path).read_text())
    phases_colors = list(baseline["colors"])
    dump(out/"initial-phase-colors.json", {"k": 5, "colors": phases_colors,
          "graph_sha256": graph_sha, "is_proper_coloring": True,
          "purpose": "SAT phase seed only; no variables frozen"})
    base = verifier.coloring_cnf(graph, 5); header, body = base.split(b"\n", 1)
    fields = header.decode().split(); u, v = pair
    additions = [[-(u*5+c+1), -(v*5+c+1)] for c in range(5)]
    data = (f"p cnf {fields[2]} {int(fields[3])+5}\n".encode()+body+
            b"".join((" ".join(map(str, clause))+" 0\n").encode() for clause in additions))
    cnf_path = out/"conditional.cnf"; cnf_path.write_bytes(data)
    config = {"goal_satisfied": False, "graph_sha256": graph_sha, "vertices": pair,
              "solver": "PySAT Glucose4", "timeout_seconds": args.timeout,
              "seed_coloring": seed_path, "preferred_endpoint_recolor": [endpoint, color],
              "preferred_endpoint_recolor_conflicts": conflicts,
              "initial_monochromatic_edges": 0, "phase_variables_frozen": False,
              "base_cnf_sha256": hashlib.sha256(base).hexdigest(),
              "conditional_cnf_sha256": hashlib.sha256(data).hexdigest(),
              "added_inequality_clauses": additions,
              "proof_scope": "conditional 5-coloring gadget only"}
    dump(out/"config.json", config)
    formula = CNF(from_file=str(cnf_path))
    phases = [var if c == phases_colors[vertex] else -var
              for vertex in range(len(phases_colors)) for c in range(5)
              for var in (vertex*5+c+1,)]
    started = time.monotonic()
    with Solver(name="g4", bootstrap_with=formula.clauses, with_proof=True) as solver:
        solver.set_phases(phases)
        timer = threading.Timer(args.timeout, solver.interrupt); timer.start()
        try:
            answer = solver.solve_limited(expect_interrupt=True)
        finally:
            timer.cancel()
        elapsed = time.monotonic()-started
        native_model = solver.get_model() if answer is True else None
        proof_lines = solver.get_proof()
        stats = solver.accum_stats()
    proof = out/"proof.drat"
    proof.write_text("\n".join(proof_lines)+("\n" if proof_lines else ""))
    result = {**config, "elapsed_seconds": elapsed, "solver_statistics": stats,
              "forced_equality_established": False, "proof_sha256": verifier.common._sha256(proof),
              "proof_bytes": proof.stat().st_size}
    dump(out/"native-outcome.json", {"solver_return_value": answer,
         "elapsed_seconds": elapsed, "statistics": stats})
    if answer is True:
        dump(out/"raw-solver-model.json", native_model)
        variable_count = len(graph["vertices"])*5
        if len(native_model) != variable_count or {abs(lit) for lit in native_model} != set(range(1, variable_count+1)):
            raise ValueError("partial or malformed solver model")
        assignment = [None]*(variable_count+1)
        for literal in native_model:
            assignment[abs(literal)] = literal > 0
        verify_full_cnf_model(cnf_path, assignment)
        colors = [next(c for c in range(5) if assignment[vertex*5+c+1])
                  for vertex in range(len(graph["vertices"]))]
        witness = {"k": 5, "colors": colors, "graph_sha256": graph_sha}
        verifier.verify_coloring(graph, witness, 5)
        if colors[u] == colors[v]:
            raise ValueError("SAT model fails conditional inequality")
        dump(out/"separating-coloring.json", witness)
        remaining = json.loads((KEMPE/"unresolved-pairs.json").read_text())
        separated = [pair for pair in remaining if colors[pair[0]] != colors[pair[1]]]
        dump(out/"additional-refuted-pairs.json", separated)
        result.update({"status": "SAT_VERIFIED", "forced_equality_disproved": True,
                       "additional_pairs_disproved_forced_equal": len(separated),
                       "pair_colors": [colors[u], colors[v]]})
    elif answer is False:
        checker = ROOT/"sources/drat-trim/drat-trim"
        checker_sha = verifier.common._sha256(checker)
        checked = subprocess.run([str(checker), str(cnf_path), str(proof)],
                                 capture_output=True, timeout=20)
        (out/"checker.stdout.log").write_bytes(checked.stdout)
        (out/"checker.stderr.log").write_bytes(checked.stderr)
        if cnf_path.read_bytes() != data or verifier.common._sha256(proof) != result["proof_sha256"]:
            raise ValueError("conditional proof inputs changed during checking")
        if verifier.common._sha256(checker) != checker_sha:
            raise ValueError("checker changed during checking")
        verified = checked.returncode == 0 and b"s VERIFIED" in checked.stdout.splitlines()
        if verified:
            verifier.verify_geometry(graph)
        result.update({"status": "CONDITIONAL_UNSAT_VERIFIED" if verified else "CONDITIONAL_UNSAT_UNVERIFIED",
                       "forced_equality_established": verified, "checker_sha256": checker_sha})
    else:
        result["status"] = "UNKNOWN"
    dump(out/"result.json", result)
    if result["status"] in ("SAT_VERIFIED", "UNKNOWN"):
        result = archive_completed_proof(out)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
