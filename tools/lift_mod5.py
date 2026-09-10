#!/usr/bin/env python3
"""Exact five-color certificates for 5-integral multiquadratic graphs.

In F_25 = F_5[w]/(w^2-2), send sqrt(2),sqrt(3),sqrt(5),sqrt(11)
to w,2w,0,1 respectively. For x=x0+x1*w and y=y0+y1*w use
color x0+y0+y1 (mod 5). Exhaustive finite arithmetic verifies that no
unit displacement has color zero. The map applies to rational coefficients
whose reduced denominators are coprime to 5. No claim about coefficients
with denominator divisible by 5 is made by this finite ring certificate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

try:
    from . import verify_mq as verifier
except ImportError:
    import verify_mq as verifier


def add(left, right):
    return ((left[0]+right[0]) % 5, (left[1]+right[1]) % 5)


def multiply(left, right):
    return ((left[0]*right[0]+2*left[1]*right[1]) % 5,
            (left[0]*right[1]+left[1]*right[0]) % 5)


IMAGES = ((0, 1), (0, 2), (0, 0), (1, 0))
BASIS_IMAGES = []
for mask in range(16):
    value = (1, 0)
    for bit, image in enumerate(IMAGES):
        if mask & (1 << bit):
            value = multiply(value, image)
    BASIS_IMAGES.append(value)


def finite_certificate():
    for radicand, image in zip((2, 3, 5, 11), IMAGES):
        if multiply(image, image) != (radicand % 5, 0):
            raise ArithmeticError("radical image violates its defining relation")
    directions = []
    for x0 in range(5):
        for x1 in range(5):
            for y0 in range(5):
                for y1 in range(5):
                    x, y = (x0, x1), (y0, y1)
                    if add(multiply(x, x), multiply(y, y)) == (1, 0):
                        color = (x0+y0+y1) % 5
                        if color == 0:
                            raise ArithmeticError("finite unit displacement has zero color")
                        directions.append({"x": x, "y": y, "color_difference": color})
    if len(directions) != 24:
        raise ArithmeticError("unexpected finite unit-direction count")
    return {
        "goal_satisfied": False,
        "purpose": "five-coloring obstruction for a coefficient ring",
        "field": "F_5[w]/(w^2-2)", "radical_images": IMAGES,
        "color": "x0+y0+y1 modulo 5",
        "enumerated_displacements": 625, "unit_directions": directions,
        "unit_direction_count": 24, "zero_color_unit_directions": 0,
        "scope": "Q(sqrt2,sqrt3,sqrt5,sqrt11) coefficients with reduced denominators coprime to 5",
    }


def reduce_coordinate(coefficients):
    value = (0, 0)
    for coefficient, image in zip(coefficients, BASIS_IMAGES):
        rational = Fraction(coefficient)
        if rational.denominator % 5 == 0:
            raise ValueError("coefficient denominator divisible by 5; this ring map does not apply")
        residue = rational.numerator * pow(rational.denominator, -1, 5) % 5
        value = add(value, (residue*image[0] % 5, residue*image[1] % 5))
    return value


def lift_graph(graph_path, output_directory):
    graph_path = Path(graph_path).resolve(strict=True)
    output_directory = Path(output_directory)
    raw = graph_path.read_bytes()
    graph = json.loads(raw)
    geometry = verifier.verify_geometry(graph)
    colors = []
    for vertex in graph["vertices"]:
        x, y = reduce_coordinate(vertex["x"]), reduce_coordinate(vertex["y"])
        colors.append((x[0]+y[0]+y[1]) % 5)
    coloring = {"k": 5, "colors": colors, "graph_sha256": hashlib.sha256(raw).hexdigest(),
                "method": "exact F25 ring map with color x0+y0+y1 mod5"}
    verifier.verify_coloring(graph, coloring, 5)
    prefix = "graph-"+coloring["graph_sha256"][:16]
    coloring_path = output_directory / (prefix+".coloring.json")
    coloring_path.write_text(json.dumps(coloring, separators=(",", ":"))+"\n")
    report = {"goal_satisfied": False, "certificate_status": "geometry-and-explicit-valid-5-coloring",
              "graph": str(graph_path), "graph_sha256": coloring["graph_sha256"],
              "geometry": geometry.__dict__, "k": 5, "all_denominators_coprime_to_5": True,
              "coloring": str(coloring_path.resolve()), "coloring_sha256": verifier.common._sha256(coloring_path),
              "all_graph_edges_checked": True}
    (output_directory/(prefix+".verification.json")).write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", action="append", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    finite = finite_certificate()
    (args.out/"finite-field-certificate.json").write_text(json.dumps(finite, indent=2, sort_keys=True)+"\n")
    reports = []
    for path in args.graph:
        report = lift_graph(path, args.out)
        reports.append(report)
        print(json.dumps(report), flush=True)
        (args.out/"index.json").write_text(json.dumps(reports, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
