#!/usr/bin/env python3
"""Independent exact geometry and certificate checker for the Voronov field.

K0=Q(sqrt2,sqrt3). Let h1=5+3sqrt2-sqrt3+3sqrt6 and
h7=5+3sqrt2+sqrt3-3sqrt6. Coordinate basis index4*m+j is
(sqrt(h1) if m&1)*(sqrt(h7) if m&2)*B_j, with
B=(1,sqrt2,sqrt3,sqrt6). All square roots use the positive real embedding.

Prime squareclasses give [K0:Q]=4. Exact rational interval signs below
show h1,h7,h1*h7 each negative at a real embedding and positive at another,
so none is a square in K0. Their squareclasses are independent, giving
degree16 and a unique rational basis. Geometry uses integer arithmetic
after clearing denominators; no construction arithmetic is imported.

Only exact geometry plus independently checked unrestricted non-6 DRAT
can satisfy the public target gate. Coloring and geometry alone cannot.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from math import gcd, lcm
from pathlib import Path
import re
import sys
from typing import Any

try:
    from . import verify as common
except ImportError:
    import verify as common


VerificationError = common.VerificationError
FIELD = {"base_radicands": [2, 3],
         "extension_squares": [[5, 3, -1, 3], [5, 3, 1, -3]]}
BASE_SQUARES = (1, 2, 3, 6)
H1 = (5, 3, -1, 3)
H7 = (5, 3, 1, -3)
DEGREE = 16
RATIONAL = re.compile(r"[+-]?[0-9]+(?:/[+]?[0-9]+)?\Z")


def base_multiply(left, right):
    output = [0]*4
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            output[i ^ j] += BASE_SQUARES[i & j]*a*b
    return tuple(output)


EXTENSION_SQUARES = ((1, 0, 0, 0), H1, H7, base_multiply(H1, H7))


def _basis_product(i, j):
    extension_i, base_i = divmod(i, 4)
    extension_j, base_j = divmod(j, 4)
    extension = extension_i ^ extension_j
    base = base_i ^ base_j
    scale = BASE_SQUARES[base_i & base_j]
    square = EXTENSION_SQUARES[extension_i & extension_j]
    return tuple((4*extension+(base ^ mask), scale*coefficient*BASE_SQUARES[base & mask])
                 for mask, coefficient in enumerate(square) if coefficient)


PRODUCTS = tuple(tuple(_basis_product(i, j) for j in range(DEGREE)) for i in range(DEGREE))


def field_basis_certificate():
    """Exact sign bounds establish the field degree and physical roots."""
    bounds = ((Fraction(1414, 1000), Fraction(1415, 1000)),
              (Fraction(1732, 1000), Fraction(1733, 1000)))
    for radicand, (lower, upper) in zip((2, 3), bounds):
        if not 0 < lower < upper or not lower*lower < radicand < upper*upper:
            raise ArithmeticError("invalid rational square-root enclosure")
    if EXTENSION_SQUARES[3] != (-14, 48, 0, 0):
        raise ArithmeticError("h1*h7 exact identity differs")
    intervals = []; signs = []
    for sign2, sign3 in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        positive = ((Fraction(1), Fraction(1)), bounds[0], bounds[1],
                    (bounds[0][0]*bounds[1][0], bounds[0][1]*bounds[1][1]))
        embedding_signs = (1, sign2, sign3, sign2*sign3)
        embedding_intervals = []; embedding_result = []
        for polynomial in (H1, H7, EXTENSION_SQUARES[3]):
            lower = upper = Fraction(0)
            for coefficient, sign, (a, b) in zip(polynomial, embedding_signs, positive):
                first, second = coefficient*sign*a, coefficient*sign*b
                lower += min(first, second); upper += max(first, second)
            if lower > 0:
                embedding_result.append(1)
            elif upper < 0:
                embedding_result.append(-1)
            else:
                raise ArithmeticError("sign enclosure is inconclusive")
            embedding_intervals.append([str(lower), str(upper)])
        signs.append(embedding_result)
        intervals.append({"sqrt2_sign": sign2, "sqrt3_sign": sign3,
                          "intervals_h1_h7_product": embedding_intervals})
    if signs[0][:2] != [1, 1] or any({row[j] for row in signs} != {-1, 1} for j in range(3)):
        raise ArithmeticError("required independent squareclasses not established")
    return {"degree": 16, "h1_h7_product": list(EXTENSION_SQUARES[3]),
            "embedding_signs_h1_h7_product": signs, "exact_sign_intervals": intervals,
            "squareclass_argument": "Each of h1,h7,h1*h7 is negative at a real embedding, hence nonsquare in K0; the two extension squareclasses are independent."}


FIELD_BASIS_CERTIFICATE = field_basis_certificate()


@dataclass(frozen=True)
class Point:
    x: tuple[int, ...]
    y: tuple[int, ...]
    denominator: int


@dataclass(frozen=True)
class GeometryResult:
    vertices: int
    edges: int
    field_degree: int = 16


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
    numerators = tuple(coefficient.numerator*(denominator//coefficient.denominator) for coefficient in coordinates)
    return Point(numerators[:16], numerators[16:], denominator)


def _difference(left: Point, right: Point) -> Point:
    denominator = lcm(left.denominator, right.denominator)
    ls, rs = denominator//left.denominator, denominator//right.denominator
    numerators = tuple(ls*a-rs*b for a, b in zip(left.x+left.y, right.x+right.y))
    divisor = gcd(denominator, *numerators)
    denominator //= divisor; numerators = tuple(value//divisor for value in numerators)
    first = next((value for value in numerators if value), 0)
    if first < 0:
        numerators = tuple(-value for value in numerators)
    return Point(numerators[:16], numerators[16:], denominator)


def _unit(point: Point) -> bool:
    norm = [0]*16
    for coordinate in (point.x, point.y):
        nonzero = [(i, value) for i, value in enumerate(coordinate) if value]
        for position, (i, a) in enumerate(nonzero):
            for j, b in nonzero[position:]:
                scale = a*b*(1 if i == j else 2)
                for target, coefficient in PRODUCTS[i][j]:
                    norm[target] += scale*coefficient
    return norm[0] == point.denominator**2 and all(value == 0 for value in norm[1:])


def verify_geometry(graph: Any) -> GeometryResult:
    if not isinstance(graph, dict) or graph.get("schema") != "hn-voronov-field-v1":
        raise VerificationError("Voronov geometry requires schema hn-voronov-field-v1")
    specification = graph.get("field")
    if not isinstance(specification, dict) or specification != FIELD:
        raise VerificationError("field metadata must exactly match the supported Voronov basis")
    if any(type(value) is not int for value in specification["base_radicands"]
           + [value for square in specification["extension_squares"] for value in square]):
        raise VerificationError("field radicands and coefficients must be integers")
    n, edges = common._structure(graph)
    points = [_point(vertex, index) for index, vertex in enumerate(graph["vertices"])]
    seen = {}
    for index, point in enumerate(points):
        if point in seen:
            raise VerificationError(f"vertices {seen[point]} and {index} are the same exact point")
        seen[point] = index
    verified = set()
    for u, v in edges:
        difference = _difference(points[u], points[v])
        if difference not in verified:
            if not _unit(difference):
                raise VerificationError(f"edge ({u}, {v}) does not have exact squared length 1")
            verified.add(difference)
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
        graph = json.loads(args.graph.read_text()); proof_arguments = (args.cnf, args.proof, args.checker)
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
            canonical = coloring_cnf(graph, args.k); args.write_cnf.write_bytes(canonical)
            report["cnf_sha256"] = hashlib.sha256(canonical).hexdigest()
        report["graph_sha256"] = common._sha256(args.graph)
        print(json.dumps(report, indent=2, sort_keys=True)); return 0
    except (VerificationError, OSError, ValueError) as exc:
        print(json.dumps({"goal_satisfied": False, "error": str(exc)}), file=sys.stderr); return 1


if __name__ == "__main__":
    sys.exit(main())
