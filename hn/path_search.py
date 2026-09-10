"""Extend Haugland's spindle assembly from G1 to the full exact T6 graph."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

from .construct import t_points, g2_points, g3_points
from .search import export
from .solve import solve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("runs/t6-spindle"))
    parser.add_argument("--timeout", type=float, default=60)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    config = {
        "family": "exact-T6-spindle",
        "definition": "T6 contains every point on an at-most-six-step A-B path using H's 84 unit directions",
        "A": "0", "B": "2*zeta42^7-1",
        "first_assembly": "zeta42^-7*T6-1 union zeta42^7*T6+1",
        "second_assembly": "S union rho*(S+1)-1",
        "rho": "(7+(2*zeta42^7-1)*sqrt(5))/8",
        "edge_acceptance": "exact squared norm equals 1",
        "candidate_enumeration": "NumPy squared-distance broad phase abs(d2-1)<1e-8; no induced completeness claim",
        "k_sequence": [6, 5], "timeout_per_solver_seconds": args.timeout,
        "argv": sys.argv,
        "source": "https://arxiv.org/html/2608.04542v1",
        "program_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (args.out / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    start = time.monotonic()
    seed = t_points(6)
    assembly = g2_points(seed)
    points = g3_points(assembly)
    print(json.dumps({"stage": "path_points", "T6": len(seed),
                      "first_assembly": len(assembly), "second_assembly": len(points),
                      "elapsed_seconds": time.monotonic()-start}), flush=True)
    graph_path = args.out / "graph.json"
    export(points, graph_path, config)
    result = solve(graph_path, 6, args.out / "six-color", args.timeout)
    if result["status"] == "SAT":
        solve(graph_path, 5, args.out / "five-color", args.timeout)


if __name__ == "__main__":
    main()
