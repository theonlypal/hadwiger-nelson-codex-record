#!/usr/bin/env python3
"""Independent exact checker for Q(zeta_42)(sqrt(5),sqrt(33)).

In hn-cyclotomic-tower-v1, parts indexed by masks 0,1,2,3 multiply
1,sqrt(5),sqrt(33),sqrt(165), respectively. Each part has twelve rational
coefficients in increasing powers of zeta_42=exp(2*pi*i/42). The square roots
are positive real numbers. Complex conjugation fixes them and inverts zeta.

Q(zeta_42) can ramify only at 2,3,7. Q(sqrt(5)) ramifies at 5, so adjoining
sqrt(5) doubles the degree. The resulting field can ramify only at 2,3,5,7;
Q(sqrt(33)), of discriminant 132, ramifies at 11. Adjoining sqrt(33) therefore
doubles the degree again. This proves uniqueness of the four-part basis and
total degree 48. Ramification cannot appear in a subfield without appearing
in its overfield. All norms below use exact rational polynomial arithmetic.

This checker imports no hn arithmetic. Only graph structure, coloring CNFs,
coloring checks and external proof checking are shared with verify.py.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any

from flint import fmpq, fmpq_poly

try:
    from . import verify as common
except ImportError:
    import verify as common


VerificationError = common.VerificationError
PHI42 = fmpq_poly([1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1])
SQUARES = (1, 5, 33, 165)
RATIONAL = re.compile(r"[+-]?[0-9]+(?:/[+]?[0-9]+)?\Z")


@dataclass(frozen=True)
class GeometryResult:
    vertices: int
    edges: int
    cyclotomic_order: int = 42
    radicands: tuple[int, ...] = (5, 33)
    field_degree: int = 48


class TowerArithmetic:
    def __init__(self):
        t = fmpq_poly([0, 1])
        inverse = t**41 % PHI42
        self.conjugate_powers = [fmpq_poly([1])]
        for _ in range(1, 12):
            self.conjugate_powers.append(self.conjugate_powers[-1]*inverse % PHI42)

    def parse(self, vertex, index):
        if not isinstance(vertex, dict) or not isinstance(vertex.get("parts"), list) or len(vertex["parts"]) != 4:
            raise VerificationError(f"vertex {index} must have four parts")
        result = []
        for part, values in enumerate(vertex["parts"]):
            if not isinstance(values, list) or len(values) != 12:
                raise VerificationError(f"vertex {index} part {part} must contain 12 rational strings")
            coefficients = []
            for value in values:
                if not isinstance(value, str) or RATIONAL.fullmatch(value) is None:
                    raise VerificationError(f"vertex {index} part {part} has a non-rational string")
                try:
                    coefficients.append(fmpq(value))
                except (ValueError, ZeroDivisionError) as exc:
                    raise VerificationError(f"vertex {index} part {part} has an invalid rational") from exc
            result.append(fmpq_poly(coefficients))
        return tuple(result)

    def conjugate(self, value):
        result = fmpq_poly([])
        for power, coefficient in enumerate(value.coeffs()):
            if coefficient:
                result += coefficient*self.conjugate_powers[power]
        return result

    def is_unit(self, parts):
        conjugates = tuple(self.conjugate(part) for part in parts)
        norm = [fmpq_poly([]) for _ in range(4)]
        for i in range(4):
            if parts[i]:
                norm[0] += SQUARES[i]*parts[i]*conjugates[i]
                for j in range(i+1, 4):
                    if parts[j]:
                        norm[i ^ j] += SQUARES[i & j]*(parts[i]*conjugates[j]+parts[j]*conjugates[i])
        norm = [part % PHI42 for part in norm]
        return norm[0] == 1 and all(part == 0 for part in norm[1:])


def verify_geometry(graph: Any) -> GeometryResult:
    if not isinstance(graph, dict) or graph.get("schema") != "hn-cyclotomic-tower-v1":
        raise VerificationError("tower geometry requires schema hn-cyclotomic-tower-v1")
    field = graph.get("field")
    if not isinstance(field, dict) or type(field.get("cyclotomic_order")) is not int or field["cyclotomic_order"] != 42:
        raise VerificationError("tower cyclotomic_order must be 42")
    radicands = field.get("radicands")
    if not isinstance(radicands, list) or radicands != [5, 33] or any(type(r) is not int for r in radicands):
        raise VerificationError("tower radicands must be exactly [5, 33] in that order")
    n, edges = common._structure(graph)
    arithmetic = TowerArithmetic()
    points = [arithmetic.parse(vertex, index) for index, vertex in enumerate(graph["vertices"])]
    seen = {}
    for index, point in enumerate(points):
        key = tuple(tuple(part.coeffs()) for part in point)
        if key in seen:
            raise VerificationError(f"vertices {seen[key]} and {index} are the same exact point")
        seen[key] = index
    certified_differences = set()
    for u, v in edges:
        difference = tuple(left-right for left, right in zip(points[u], points[v]))
        first = next(coefficient for part in difference for coefficient in part.coeffs() if coefficient)
        if first < 0:
            difference = tuple(-part for part in difference)
        key = tuple(tuple(part.coeffs()) for part in difference)
        if key not in certified_differences:
            if not arithmetic.is_unit(difference):
                raise VerificationError(f"edge ({u}, {v}) does not have exact squared length 1")
            certified_differences.add(key)
    return GeometryResult(n, len(edges))


coloring_cnf = common.coloring_cnf
verify_coloring = common.verify_coloring


def verify_non6_certificate(graph, cnf_path, proof_path, checker_path, k=6, timeout=None):
    if type(k) is not int or k != 6:
        raise VerificationError("the target certificate must establish non-6-colorability")
    geometry = verify_geometry(graph)
    proof = common._verify_non6_proof(graph, cnf_path, proof_path, checker_path, k, timeout)
    return {"goal_satisfied": True, "certificate_status": "exact-unit-distance-and-non-6-colorable",
            "geometry": geometry.__dict__, **proof}


def main(argv=None):
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
            report = {"goal_satisfied": False, "certificate_status": "geometry-only",
                      "geometry": verify_geometry(graph).__dict__}
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
    except (VerificationError, ValueError, OSError) as exc:
        print(json.dumps({"goal_satisfied": False, "error": str(exc)}, sort_keys=True), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
