"""Exact reconstructions and finite search constructors.

Haugland's Appendix A is used only as a list of candidate coordinates. Every
exported edge is checked by exact algebraic arithmetic. The original paper's
floating-point path-merging lemma is not an assumption of this reconstruction.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path

from .field import Element, FIELD_METADATA, ONE, ZERO, ZETA, SQRT5

ROOT = Path(__file__).resolve().parents[1]
SOURCE_URL = "https://arxiv.org/html/2608.04542v1"
SOURCE_PATH = ROOT / "sources" / "haugland_2608.04542v1.html"


@lru_cache(maxsize=1)
def unit_directions():
    """u[0] ... u[83], exactly in the ordering of Haugland's Table 1."""
    # u1=(Q0-R0), reduced modulo Phi_42, then u_(2j)=t^j,
    # u_(2j+1)=u1*t^j. No arctangent or floating-point angle is used.
    odd = Element((2, 1, 3, 0, -2, 1, 5, 1, 1, -6, -4, 4), 7)
    directions = []
    power = ONE
    for j in range(42):
        directions.extend((power, odd * power))
        power *= ZETA
    assert len(set(directions)) == 84
    assert all(v.is_unit() for v in directions)
    return tuple(directions)


def h_points():
    """H's 21 vertices, listed P0...P6,Q0...Q6,R0...R6."""
    p = Element((-2, -3, 0, 5, -1, 0, -1, 0, 3, 4, 0, 1), 7)
    q = Element((-3, 0, 2, 5, 0, -4, 0, 3, 6, 0, -5, -1), 7)
    r = Element((-5, -1, -1, 5, 2, -5, -5, 2, 5, 6, -1, -5), 7)
    powers = tuple(ZETA ** (6*j) for j in range(7))
    return tuple(z * power for z in (p, q, r) for power in powers)


def load_appendix_paths():
    """Read the acquired HTML, retaining only the Appendix A integer rows."""
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(SOURCE_PATH.read_text(), "html.parser")
    paths = []
    for table in soup.find_all("table"):
        for row in table.find_all("tr"):
            cells = [c.get_text(" ", strip=True) for c in row.find_all("td")]
            if len(cells) in (5, 6) and all(s.isdecimal() for s in cells):
                path = tuple(int(s) for s in cells)
                if all(0 <= j < 84 for j in path):
                    paths.append(path)
    if not paths:
        raise ValueError("No Appendix A paths found")
    return tuple(paths)


def load_g1_points():
    points = {ZERO: None}
    directions = unit_directions()
    target = 2 * ZETA**7 - ONE  # i*sqrt(3)
    for path in load_appendix_paths():
        position = ZERO
        for j in path:
            position += directions[j]
            points[position] = None
        if position != target:
            raise ArithmeticError(f"Appendix path does not end at B: {path}")
    return tuple(points)


def g2_points(points=None):
    if points is None:
        points = load_g1_points()
    return tuple(dict.fromkeys(
        itertools.chain((ZETA**(-7) * z - ONE for z in points),
                        (ZETA**7 * z + ONE for z in points))))


def spindle_rotation():
    # i*sqrt(15) = (i*sqrt(3))*sqrt(5), in the selected complex embedding.
    rotation = (7 * ONE + (2 * ZETA**7 - ONE) * SQRT5) / 8
    assert rotation.is_unit()
    return rotation


def g3_points(points=None):
    if points is None:
        points = g2_points()
    rotation = spindle_rotation()
    return tuple(dict.fromkeys(
        itertools.chain(points, (rotation * (z + ONE) - ONE for z in points))))


@lru_cache(maxsize=2)
def t_points(n=6):
    """All points on <=n-step A--B paths, for n in {5,6}, exactly.

    Meet two tables of unordered half-path sums of length <=3. Each matching
    multiset represents every ordering of those steps; all its subset sums
    are precisely the vertices visited in all those orderings. This removes
    the tolerance assumptions used for T5/T6 in the source paper.
    """
    if n not in (5, 6):
        raise ValueError("This bounded constructor supports n=5 or n=6")
    directions = unit_directions()
    table = {ZERO: [()]}
    for length in (1, 2, 3):
        for indices in itertools.combinations_with_replacement(range(84), length):
            position = sum((directions[j] for j in indices), ZERO)
            table.setdefault(position, []).append(indices)
    target = 2 * ZETA**7 - ONE
    paths = set()
    for position, left_paths in table.items():
        right_paths = table.get(target-position, ())
        for left in left_paths:
            for right in right_paths:
                if len(left)+len(right) <= n:
                    paths.add(tuple(sorted(left+right)))
    points = {ZERO}
    for indices in sorted(paths):
        subsets = {ZERO}
        for j in indices:
            subsets |= {z+directions[j] for z in subsets}
        points.update(subsets)
    # Stable order: canonical integer tuple and common denominator.
    return tuple(sorted(points, key=lambda z: (z.denominator, z.coefficients)))


def graph_from_points(points, *, directions=None, metadata=None):
    """Export distinct points and exact unit edges.

    With directions given, perform exact membership lookup for those edges.
    Otherwise use a floating-point spatial prefilter then validate *every*
    retained edge with exact norm equality. Omitting an edge is harmless for
    a non-colorability witness, so no induced-graph claim is made. Geometric
    correctness never depends on the prefilter threshold.
    """
    points = tuple(dict.fromkeys(points))
    edges = set()
    rejected = 0
    if directions is not None:
        directions = tuple(dict.fromkeys(directions))
        if not all(z.is_unit() for z in directions):
            raise ValueError("Non-unit search direction")
        index = {z: j for j, z in enumerate(points)}
        for i, z in enumerate(points):
            for u in directions:
                j = index.get(z + u)
                if j is not None and j != i:
                    edges.add((min(i, j), max(i, j)))
        edge_search = "exact lookup using certified unit directions"
    else:
        approximate = [z.approximate() for z in points]
        buckets = {}
        for i, z in enumerate(approximate):
            bx, by = math.floor(z.real), math.floor(z.imag)
            # Radius 2 is deliberate: a candidate exactly on a cell boundary
            # is not lost because a double rounds to the adjacent cell.
            for dx in range(-2, 3):
                for dy in range(-2, 3):
                    for j in buckets.get((bx+dx, by+dy), ()):
                        delta = z - approximate[j]
                        if abs(delta.real**2 + delta.imag**2 - 1) < 1e-6:
                            if (points[i] - points[j]).is_unit():
                                edges.add((j, i))
                            else:
                                rejected += 1
            buckets.setdefault((bx, by), []).append(i)
        edge_search = "approximate candidates with exact norm verification"
    detail = dict(metadata or {})
    detail.update({"edge_search": edge_search,
                   "induced_graph_claimed": False,
                   "nonunit_candidates_rejected": rejected,
                   "vertex_count": len(points), "edge_count": len(edges)})
    return {"schema": "hn-exact-v1", "field": dict(FIELD_METADATA),
            "vertices": [z.to_json() for z in points],
            "edges": [list(e) for e in sorted(edges)], "metadata": detail}


def graph_points(graph):
    return tuple(Element.from_json(obj) for obj in graph["vertices"])


def write_graph(graph, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(graph, sort_keys=True, separators=(",", ":")) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph", choices=("h", "g1", "g2", "g3", "t5", "t6"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    constructors = {"h": h_points, "g1": load_g1_points,
                    "g2": g2_points, "g3": g3_points,
                    "t5": lambda: t_points(5), "t6": lambda: t_points(6)}
    metadata = {"name": "Haugland " + args.graph.upper(),
                "source": SOURCE_URL,
                "source_sha256": hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest(),
                "source_license": "CC BY 4.0",
                "status": "candidate; no chromatic conclusion asserted"}
    graph = graph_from_points(constructors[args.graph](), metadata=metadata)
    write_graph(graph, args.output)
    print(json.dumps({"path": str(args.output), "vertices": len(graph["vertices"]),
                      "edges": len(graph["edges"])}))


if __name__ == "__main__":
    main()
