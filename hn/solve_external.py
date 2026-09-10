"""Native CaDiCaL runs with exact geometry/model/proof verification.

Every run keeps the original canonical CNF, unedited stdout/stderr, native
solution file, and streamed binary DRAT proof, including partial timeout
proofs. Exit 10 must yield a complete, independently checked model; exit 20
must pass independent geometric and DRAT checks. Only a verified six-color
refutation satisfies the target; optional five-color runs are diagnostics.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import platform
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

from tools import verify as common
from tools import verify_mq
from tools import verify_tower


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOLVER = ROOT / "sources/cadical/build/cadical"
DEFAULT_CHECKER = ROOT / "sources/drat-trim/drat-trim"


def _dump(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+"\n")


def compress_completed_run(out: str | Path) -> dict:
    """Losslessly compress a completed diagnostic proof; preserve raw hash.

    A native exit code in the final result is required. The decompressed
    bytes and the still-present original are checked before the original is
    removed. This function never operates on an active run's proof.
    """
    out = Path(out).resolve()
    result_path = out / "result.json"
    result = json.loads(result_path.read_text())
    if result.get("status") not in (
        "SAT", "UNKNOWN", "UNKNOWN_TIMEOUT", "UNKNOWN_NATIVE", "STOPPED_ALGEBRAIC_5_COLORING",
        "ERROR_NATIVE_EXIT", "NON5_VERIFIED"
    ) or "native_exit_code" not in result:
        raise ValueError("proof compression requires a completed native diagnostic run")
    if "proof_storage" in result:
        return result
    proof = out / "proof.drat"
    if not proof.exists():
        return result
    original_sha = result["proof_sha256"]
    if common._sha256(proof) != original_sha:
        raise ValueError("completed proof hash changed before compression")
    compressed = out / "proof.drat.gz"
    with compressed.open("xb") as output:
        with gzip.GzipFile(filename="", mode="wb", fileobj=output, compresslevel=6, mtime=0) as zipped:
            with proof.open("rb") as original:
                shutil.copyfileobj(original, zipped, length=1024*1024)
    recovered = hashlib.sha256()
    recovered_bytes = 0
    with gzip.open(compressed, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            recovered.update(chunk)
            recovered_bytes += len(chunk)
    if recovered.hexdigest() != original_sha or recovered_bytes != result["proof_bytes"]:
        raise ValueError("compressed proof does not reproduce the original bytes")
    if common._sha256(proof) != original_sha:
        raise ValueError("original proof changed during compression")
    result["proof_storage"] = {
        "path": str(compressed), "compression": "gzip", "decompression_verified": True,
        "uncompressed_sha256": original_sha, "uncompressed_bytes": recovered_bytes,
        "compressed_sha256": common._sha256(compressed), "compressed_bytes": compressed.stat().st_size,
    }
    # Persist the verified compressed artifact's identity before unlinking.
    _dump(result_path, result)
    proof.unlink()
    return result


def parse_model(path: Path, variable_count: int) -> list[bool]:
    """Reject partial, contradictory, out-of-range, or unterminated models."""
    assignment = [None] * (variable_count+1)
    terminated = False
    statuses = []
    for line in path.read_text().splitlines():
        if line.startswith("s "):
            statuses.append(line.strip())
        elif line.startswith("v "):
            for token in line[2:].split():
                literal = int(token)
                if terminated:
                    raise ValueError("model has literals after its terminating zero")
                if literal == 0:
                    terminated = True
                    continue
                variable = abs(literal)
                if variable > variable_count:
                    raise ValueError("model contains an out-of-range variable")
                if assignment[variable] is not None:
                    raise ValueError("model assigns a variable more than once")
                assignment[variable] = literal > 0
    if statuses != ["s SATISFIABLE"] or not terminated:
        raise ValueError("model has no unique SAT status and terminating zero")
    if any(value is None for value in assignment[1:]):
        raise ValueError("model is incomplete")
    return assignment


def verify_full_cnf_model(cnf: Path, assignment: list[bool]) -> int:
    checked = 0
    with cnf.open() as stream:
        header = stream.readline().split()
        if len(header) != 4 or header[:2] != ["p", "cnf"] or int(header[2]) != len(assignment)-1:
            raise ValueError("unexpected canonical CNF header")
        expected_clauses = int(header[3])
        for line in stream:
            clause = list(map(int, line.split()))
            if not clause or clause[-1] != 0 or 0 in clause[:-1]:
                raise ValueError("unexpected canonical CNF clause")
            if not any(assignment[abs(literal)] == (literal > 0) for literal in clause[:-1]):
                raise ValueError(f"native model falsifies CNF clause {checked}")
            checked += 1
    if checked != expected_clauses:
        raise ValueError("CNF clause count does not match its header")
    return checked


def solve_external(graph_path, out, timeout=300.0, solver=DEFAULT_SOLVER,
                   checker=DEFAULT_CHECKER, seed=0, configuration="default", k=6,
                   phase_file=None):
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be positive and finite")
    if configuration not in ("default", "sat", "unsat"):
        raise ValueError("unsupported CaDiCaL configuration")
    if type(k) is not int or k not in (5, 6):
        raise ValueError("native color count must be 5 (diagnostic) or 6 (target)")
    graph_path = Path(graph_path).resolve(strict=True)
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    solver = Path(solver).resolve(strict=True)
    checker = Path(checker).resolve(strict=True)
    raw_graph = graph_path.read_bytes()
    graph = json.loads(raw_graph)
    verifier = {
        "hn-exact-v1": common,
        "hn-multiquadratic-v1": verify_mq,
        "hn-cyclotomic-tower-v1": verify_tower,
    }.get(graph.get("schema"))
    if verifier is None:
        raise ValueError("unsupported exact geometry schema")
    # Validate geometry before committing any compute to a candidate.
    geometry = verifier.verify_geometry(graph)
    cnf = out / "coloring.cnf"
    canonical = common.coloring_cnf(graph, k)
    cnf.write_bytes(canonical)
    proof = out / "proof.drat"
    solution = out / "solution.out"
    stdout_path, stderr_path = out / "solver.stdout.log", out / "solver.stderr.log"
    command = [str(solver), "--no-colors", f"--seed={seed}"]
    if configuration != "default":
        command.append("--"+configuration)
    command.extend(["-w", str(solution), str(cnf), str(proof)])
    phase_path = Path(phase_file).resolve(strict=True) if phase_file is not None else None
    if phase_path is not None:
        command.extend(["--phase-file", str(phase_path)])
    version = subprocess.run([str(solver), "--version"], text=True, capture_output=True, check=True).stdout.strip()
    provenance = {
        "schema": "hn-native-sat-run-v1", "goal_satisfied": False,
        "success_certificate": False, "solver": "cadical", "solver_version": version,
        "solver_path": str(solver), "solver_sha256": common._sha256(solver),
        "graph": str(graph_path), "graph_sha256": hashlib.sha256(raw_graph).hexdigest(),
        "cnf_sha256": hashlib.sha256(canonical).hexdigest(), "geometry": geometry.__dict__,
        "vertices": geometry.vertices, "edges": geometry.edges, "k": k,
        "command": command, "timeout_seconds": timeout, "seed": seed,
        "configuration": configuration, "proof_format": "binary DRAT",
        "python": sys.version, "platform": platform.platform(),
    }
    if phase_path is not None:
        provenance['phase_hint_file'] = str(phase_path)
        provenance['phase_hint_sha256'] = common._sha256(phase_path)
        provenance['phase_hint_scope'] = 'Search preferences only; canonical CNF unchanged'
    # Record source identity separately from the executable hash where available.
    source = solver.parent.parent
    source_revision = subprocess.run(["git", "-C", str(source), "rev-parse", "HEAD"], capture_output=True, text=True)
    if source_revision.returncode == 0:
        provenance["solver_source_revision"] = source_revision.stdout.strip()
    _dump(out/"config.json", provenance)
    print(json.dumps({"stage": "native_solver_started", "out": str(out), **provenance}), flush=True)
    result = dict(provenance)
    started = time.monotonic()
    timed_out = False
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            process.send_signal(signal.SIGTERM)
            try:
                returncode = process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                returncode = process.wait()
    result.update({"elapsed_seconds": time.monotonic()-started, "native_exit_code": returncode,
                   "timed_out": timed_out, "stdout_sha256": common._sha256(stdout_path),
                   "stderr_sha256": common._sha256(stderr_path)})
    if proof.exists():
        result.update({"proof_sha256": common._sha256(proof), "proof_bytes": proof.stat().st_size})
    try:
        if common._sha256(graph_path) != provenance["graph_sha256"]:
            raise ValueError("source graph changed during native solve")
        if cnf.read_bytes() != canonical:
            raise ValueError("canonical CNF changed during native solve")
        if common._sha256(solver) != provenance["solver_sha256"]:
            raise ValueError("native solver executable changed during solve")
        if phase_path is not None and common._sha256(phase_path) != provenance['phase_hint_sha256']:
            raise ValueError('Phase hint file changed during solve')
        if timed_out:
            result["status"] = "UNKNOWN_TIMEOUT"
            result["partial_proof_retained"] = proof.exists()
        elif returncode == 10:
            assignment = parse_model(solution, k*geometry.vertices)
            checked = verify_full_cnf_model(cnf, assignment)
            colors = []
            for vertex in range(geometry.vertices):
                selected = [color for color in range(k) if assignment[k*vertex+color+1]]
                if len(selected) != 1:
                    raise ValueError("model does not select exactly one color per vertex")
                colors.append(selected[0])
            coloring = {"k": k, "colors": colors, "graph_sha256": provenance["graph_sha256"]}
            verifier.verify_coloring(graph, coloring, k)
            _dump(out/"coloring.json", coloring)
            verification = {"goal_satisfied": False, "certificate_status": f"geometry-and-valid-{k}-coloring",
                            "geometry": geometry.__dict__, "graph_sha256": provenance["graph_sha256"],
                            "coloring_sha256": common._sha256(out/"coloring.json"),
                            "complete_native_model_checked": True, "cnf_clauses_checked": checked}
            _dump(out/"verification.json", verification)
            result.update({"status": "SAT", "proper_coloring_checked": True,
                           "complete_native_model_checked": True, "cnf_clauses_checked": checked})
        elif returncode == 20:
            # A native UNSAT result is only a trigger for independent checking.
            result["status"] = "UNSAT_UNVERIFIED"
            _dump(out/"result.json", result)
            if k == 6:
                verified = verifier.verify_non6_certificate(graph, cnf, proof, checker, k=6)
            else:
                diagnostic_geometry = verifier.verify_geometry(graph)
                diagnostic_proof = common._verify_coloring_proof(graph, cnf, proof, checker, k=5)
                verified = {"goal_satisfied": False, "certificate_status": "exact-unit-distance-and-non-5-colorable",
                            "geometry": diagnostic_geometry.__dict__, "chromatic_lower_bound": 6,
                            **diagnostic_proof}
            _dump(out/"verification.json", verified)
            result.update({"status": "UNSAT_INDEPENDENTLY_VERIFIED" if k == 6 else "NON5_VERIFIED",
                           "goal_satisfied": k == 6, "success_certificate": k == 6})
        elif returncode == 0:
            result["status"] = "UNKNOWN_NATIVE"
            result["partial_proof_retained"] = proof.exists()
        else:
            result["status"] = "ERROR_NATIVE_EXIT"
    except (ValueError, OSError) as exc:
        result.update({"status": "ERROR_VERIFICATION", "error": str(exc),
                       "goal_satisfied": False, "success_certificate": False})
        _dump(out/"result.json", result)
        print(json.dumps(result), flush=True)
        raise
    _dump(out/"result.json", result)
    if result["status"] in ("SAT", "UNKNOWN_TIMEOUT", "UNKNOWN_NATIVE", "NON5_VERIFIED"):
        result = compress_completed_run(out)
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--solver", type=Path, default=DEFAULT_SOLVER)
    parser.add_argument("--checker", type=Path, default=DEFAULT_CHECKER)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--configuration", choices=("default", "sat", "unsat"), default="default")
    parser.add_argument("--k", type=int, choices=(5, 6), default=6)
    args = parser.parse_args()
    solve_external(args.graph, args.out, args.timeout, args.solver, args.checker, args.seed, args.configuration, k=args.k)


if __name__ == "__main__":
    main()
