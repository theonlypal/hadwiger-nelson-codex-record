"""Exact Moser/de Grey/Heule seeds and extensions containing sqrt(2)."""
from __future__ import annotations

import argparse
from functools import lru_cache
import itertools
import json
from pathlib import Path
import time

from .multiquadratic import MQ, Point, RADICANDS, ZERO, ONE, SQRT2, SQRT3, SQRT5, SQRT11, ORIGIN, EAST

SOURCE = "https://arxiv.org/html/1805.12181v1"
HEX_ROTATION = Point(ONE/2, SQRT3/2)
HALF_THETA3 = Point(SQRT3*SQRT11/6, SQRT3/6)
THETA3 = Point(5*ONE/6, SQRT11/6)
THETA4 = Point(7*ONE/8, SQRT3*SQRT5/8)
SQRT2_ROTATION = Point(SQRT2/2, SQRT2/2)
SQRT2_RATIONAL_ROTATION = Point(ONE/3, 2*SQRT2/3)
ROTATIONS = {"hex": HEX_ROTATION, "half-theta3": HALF_THETA3,
             "theta3": THETA3, "theta4": THETA4,
             "sqrt2-diagonal": SQRT2_ROTATION,
             "sqrt2-rational": SQRT2_RATIONAL_ROTATION}


def canonical_points(points):
    return tuple(sorted(set(points), key=lambda p: (p.x.denominator, p.x.coefficients,
                                                    p.y.denominator, p.y.coefficients)))


@lru_cache(maxsize=1)
def unit_directions_v31():
    directions = canonical_points(HEX_ROTATION**j * HALF_THETA3**k
                                  for j in range(6) for k in range(-2, 3))
    assert len(directions) == 30 and all(d.is_unit() for d in directions)
    return directions


def v31_points():
    return (ORIGIN,) + unit_directions_v31()


@lru_cache(maxsize=1)
def v151_points():
    sums = {a+b for a in v31_points() for b in v31_points()}
    clipped = canonical_points(p for p in sums if p.norm_squared() <= ONE)
    assert len(clipped) == 151
    return clipped


@lru_cache(maxsize=1)
def v1939_points():
    points = canonical_points(a+b for a in v31_points() for b in v151_points())
    assert len(points) == 1939
    return points


def moser_spindle_points():
    diamond = (ORIGIN, EAST, HEX_ROTATION, EAST+HEX_ROTATION)
    points = canonical_points(itertools.chain(diamond, (THETA3*p for p in diamond)))
    assert len(points) == 7
    return points


def heule_union_points():
    points = v1939_points()
    return canonical_points(itertools.chain(points, (THETA4*p for p in points)))


def rotated_union(points, rotation, copies=2):
    if not rotation.is_unit():
        raise ValueError("Rotation multiplier must have exact norm 1")
    return canonical_points(rotation**j*p for j in range(copies) for p in points)


def translated_unit_star(points, rotation=SQRT2_ROTATION):
    offsets = (ORIGIN,) + tuple(rotation*HEX_ROTATION**j for j in range(6))
    return canonical_points(p+u for p in points for u in offsets)


def graph_from_points(points, *, directions=None, metadata=None):
    from .search import all_candidate_edges
    points = canonical_points(points)
    started = time.monotonic()
    if directions is None:
        edges = all_candidate_edges(points)
    else:
        directions = canonical_points(directions)
        if not all(p.is_unit() for p in directions):
            raise ValueError("An edge search direction is not unit length")
        index = {p: i for i, p in enumerate(points)}
        pairs = set()
        for i, p in enumerate(points):
            for direction in directions:
                j = index.get(p+direction)
                if j is not None and j != i:
                    pairs.add((min(i, j), max(i, j)))
        edges = [list(e) for e in sorted(pairs)]
    return {"schema": "hn-multiquadratic-v1", "radicands": list(RADICANDS),
            "vertices": [p.to_json() for p in points], "edges": edges,
            "metadata": dict(metadata or {}, source=SOURCE,
                             edge_rule="exact squared Euclidean distance 1",
                             induced_graph_claimed=False,
                             edge_enumeration_seconds=time.monotonic()-started)}


def export(points, out, metadata=None, directions=None):
    graph = graph_from_points(points, directions=directions, metadata=metadata)
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(graph, separators=(",", ":")) + "\n")
    print(json.dumps({"stage": "graph_written", "path": str(path),
                      "vertices": len(graph["vertices"]), "edges": len(graph["edges"])}), flush=True)
    return graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", choices=("moser", "v31", "v151", "v1939", "heule"))
    parser.add_argument("--rotation", choices=tuple(ROTATIONS))
    parser.add_argument("--copies", type=int, default=2)
    parser.add_argument("--translated-star", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    constructors = {"moser": moser_spindle_points, "v31": v31_points,
                    "v151": v151_points, "v1939": v1939_points, "heule": heule_union_points}
    points = constructors[args.seed]()
    if args.rotation:
        points = rotated_union(points, ROTATIONS[args.rotation], args.copies)
    if args.translated_star:
        points = translated_unit_star(points)
    export(points, args.output, {"seed": args.seed, "rotation": args.rotation,
                                "copies": args.copies, "translated_star": args.translated_star,
                                "target": "non-6-colorability; no success asserted"})


if __name__ == "__main__":
    main()
