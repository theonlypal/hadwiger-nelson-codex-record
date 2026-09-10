#!/usr/bin/env python3
"""Independent exact checker for hn-multiquadratic-v1 unit-distance graphs.

For the strictly supported radicands [2, 3, 5, 11], coefficient mask m means
sqrt(product of primes whose bits are set in m), in the positive real
embedding. Distinct prime valuations prove these four rational squareclasses
are independent. The resulting multiquadratic extension has degree 2**4 = 16,
so this is a unique rational basis. Products of two basis elements satisfy
B_i*B_j = product(primes in i & j) * B_(i xor j).

Coordinates and edge norms are evaluated with integer arithmetic after
clearing denominators. No construction module is imported. The combinatorial
CNF and external DRAT checks are shared with verify.py; both public success
gates verify the appropriate exact geometry before examining a proof.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from math import gcd, lcm, prod
from pathlib import Path
import re
import sys
from typing import Any

try:
    from . import verify as common
except ImportError:
    import verify as common


VerificationError = common.VerificationError
RADICANDS = [2, 3, 5, 11]
DEGREE = 16
BASIS_SQUARES = tuple(
    prod(prime for bit, prime in enumerate(RADICANDS) if mask & (1 << bit))
    for mask in range(DEGREE)
)
RATIONAL = re.compile(r"[+-]?[0-9]+(?:/[+]?[0-9]+)?\Z")


@dataclass(frozen=True)
class Point:
    x: tuple[int, ...]
    y: tuple[int, ...]
    denominator: int


@dataclass(frozen=True)
class GeometryResult:
    vertices: int
    edges: int
    radicands: tuple[int, ...] = (2, 3, 5, 11)
    field_degree: int = DEGREE


def _point(vertex: Any, index: int) -> Point:
    if not isinstance(vertex, dict):
        raise VerificationError(f"vertex {index} must be an object")
    coordinates = []
    for axis in ("x", "y"):
        values = vertex.get(axis)
        if not isinstance(values, list) or len(values) != DEGREE:
            raise VerificationError(f"vertex {index}.{axis} must contain 16 rational strings")
        for value in values:
            if not isinstance(value, str) or RATIONAL.fullmatch(value) is None:
                raise VerificationError(f"vertex {index}.{axis} has a non-rational string")
            try:
                coordinates.append(Fraction(value))
            except (ValueError, ZeroDivisionError) as exc:
                raise VerificationError(f"vertex {index}.{axis} has an invalid rational") from exc
    denominator = lcm(*(coefficient.denominator for coefficient in coordinates))
    numerators = tuple(
        coefficient.numerator * (denominator // coefficient.denominator)
        for coefficient in coordinates
    )
    return Point(numerators[:DEGREE], numerators[DEGREE:], denominator)


def _difference(left: Point, right: Point) -> Point:
    denominator = lcm(left.denominator, right.denominator)
    ls, rs = denominator // left.denominator, denominator // right.denominator
    numerators = tuple(ls*a-rs*b for a, b in zip(left.x+left.y, right.x+right.y))
    divisor = gcd(denominator, *numerators)
    if divisor != 1:
        numerators = tuple(value // divisor for value in numerators)
        denominator //= divisor
    # Squared norms are unchanged by simultaneous sign reversal. This is an
    # exact cache key normalization, not an approximation or geometric test.
    first = next((value for value in numerators if value), 0)
    if first < 0:
        numerators = tuple(-value for value in numerators)
    return Point(numerators[:DEGREE], numerators[DEGREE:], denominator)


def _unit(point: Point) -> bool:
    norm = [0] * DEGREE
    for coordinate in (point.x, point.y):
        nonzero = [(mask, coefficient) for mask, coefficient in enumerate(coordinate) if coefficient]
        for offset, (mask, coefficient) in enumerate(nonzero):
            norm[0] += BASIS_SQUARES[mask] * coefficient * coefficient
            for other_mask, other_coefficient in nonzero[offset+1:]:
                norm[mask ^ other_mask] += (
                    2 * BASIS_SQUARES[mask & other_mask] * coefficient * other_coefficient
                )
    return norm[0] == point.denominator**2 and all(value == 0 for value in norm[1:])


def verify_geometry(graph: Any) -> GeometryResult:
    if not isinstance(graph, dict) or graph.get("schema") != "hn-multiquadratic-v1":
        raise VerificationError("multiquadratic geometry requires schema hn-multiquadratic-v1")
    radicands = graph.get("radicands")
    if not isinstance(radicands, list) or radicands != RADICANDS or any(type(p) is not int for p in radicands):
        raise VerificationError("radicands must be exactly [2, 3, 5, 11] in this order")
    n, edges = common._structure(graph)
    points = [_point(vertex, index) for index, vertex in enumerate(graph["vertices"])]
    seen: dict[Point, int] = {}
    for index, point in enumerate(points):
        if point in seen:
            raise VerificationError(f"vertices {seen[point]} and {index} are the same exact point")
        seen[point] = index
    verified_differences = set()
    for u, v in edges:
        difference = _difference(points[u], points[v])
        if difference not in verified_differences:
            if not _unit(difference):
                raise VerificationError(f"edge ({u}, {v}) does not have exact squared length 1")
            verified_differences.add(difference)
    return GeometryResult(n, len(edges))


coloring_cnf = common.coloring_cnf
verify_coloring = common.verify_coloring


def verify_non6_certificate(
    graph: Any,
    cnf_path: str | Path,
    proof_path: str | Path,
    checker_path: str | Path,
    k: int = 6,
    timeout: float | None = None,
) -> dict[str, Any]:
    if type(k) is not int or k != 6:
        raise VerificationError("the target certificate must establish non-6-colorability")
    geometry = verify_geometry(graph)
    proof = common._verify_non6_proof(graph, cnf_path, proof_path, checker_path, k, timeout)
    return {
        "goal_satisfied": True,
        "certificate_status": "exact-unit-distance-and-non-6-colorable",
        "geometry": geometry.__dict__,
        **proof,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", required=True, type=Path)
    parser.add_argument("--k", type=int, default=6)
    parser.add_argument("--coloring", type=Path)
    parser.add_argument("--write-cnf", type=Path)
    parser.add_argument("--cnf", type=Path)
    parser.add_argument("--proof", type=Path)
    parser.add_argument("--checker", type=Path)
    parser.add_argument("--checker-timeout", type=float)
    args = parser.parse_args(argv)
    try:
        graph = json.loads(args.graph.read_text())
        proof_arguments = (args.cnf, args.proof, args.checker)
        if any(proof_arguments) and not all(proof_arguments):
            raise VerificationError("a non-6 certificate requires --cnf, --proof, and --checker together")
        if all(proof_arguments):
            report = verify_non6_certificate(graph, args.cnf, args.proof, args.checker, args.k, args.checker_timeout)
        else:
            report = {
                "goal_satisfied": False,
                "certificate_status": "geometry-only",
                "geometry": verify_geometry(graph).__dict__,
            }
        if args.coloring:
            report["coloring"] = verify_coloring(graph, json.loads(args.coloring.read_text()), args.k)
            if report["goal_satisfied"]:
                raise VerificationError("a valid 6-coloring contradicts the claimed non-6 certificate")
            report["certificate_status"] = "geometry-and-valid-coloring"
        if args.write_cnf:
            canonical = coloring_cnf(graph, args.k)
            args.write_cnf.write_bytes(canonical)
            report["cnf_sha256"] = hashlib.sha256(canonical).hexdigest()
        report["graph_sha256"] = common._sha256(args.graph)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    except (VerificationError, OSError, ValueError) as exc:
        print(json.dumps({"goal_satisfied": False, "error": str(exc)}, sort_keys=True), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
