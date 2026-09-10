"""Exact unit-circle intersection closure introducing sqrt(2).

For center separation squared s in {4/3,8/3}, common unit neighbors are
(a+b)/2 +/- i*(b-a)*sqrt((4-s)/s)/2. Both points are retained, with parent
provenance. Approximation proposes center pairs; exact rational norm equality
and four exact unit-distance checks decide acceptance.
"""

from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .multiquadratic import MQ, Point, ONE, SQRT2
from .mq_construct import canonical_points, heule_union_points, v1939_points, export
from .solve import solve

SEPARATIONS = {Fraction(4, 3): SQRT2 / 2,
               Fraction(8, 3): SQRT2 / 4}


def common_unit_neighbors(a, b, separation):
    separation = Fraction(separation)
    if separation not in SEPARATIONS:
        raise ValueError("Unsupported squared separation")
    d = b - a
    if d.norm_squared() != MQ(separation):
        raise ValueError("Centers do not have the claimed exact separation")
    midpoint = (a + b) / 2
    offset = Point(-d.y, d.x) * SEPARATIONS[separation]
    neighbors = (midpoint + offset, midpoint - offset)
    for c in neighbors:
        if not (c-a).is_unit() or not (c-b).is_unit():
            raise ArithmeticError("Circle-intersection norm identity failed")
    return neighbors


def intersect_round(points, separations, max_vertices):
    points = tuple(points)
    original_n = len(points)
    index = {p: j for j, p in enumerate(points)}
    output = list(points)
    provenance = []
    pair_counts = {str(s): 0 for s in separations}
    rejected_approximate_pairs = 0
    coords = np.array([[p.x.approximate(), p.y.approximate()] for p in points])
    for a, p in enumerate(points):
        differences = coords[a+1:] - coords[a]
        norms = np.einsum('ij,ij->i', differences, differences)
        for s in separations:
            for offset in np.flatnonzero(np.abs(norms-float(s)) < 1e-8):
                b = a+1+int(offset)
                if (points[b]-p).norm_squared() != MQ(s):
                    rejected_approximate_pairs += 1
                    continue
                pair_counts[str(s)] += 1
                neighbors = common_unit_neighbors(p, points[b], s)
                if len(output) + len(set(neighbors) - index.keys()) > max_vertices:
                    return output, provenance, {
                        "truncated_at_vertex_cap": True,
                        "original_vertices": original_n,
                        "new_vertices": len(output)-original_n,
                        "center_pairs": pair_counts,
                        "approximate_pairs_rejected": rejected_approximate_pairs,
                    }
                for sign, c in zip((1, -1), neighbors):
                    if c not in index:
                        index[c] = len(output)
                        output.append(c)
                        provenance.append({"new_generation_index": index[c],
                                           "center_generation_indices": [a, b],
                                           "separation_squared": str(s), "sign": sign})
    return output, provenance, {"truncated_at_vertex_cap": False,
                               "original_vertices": original_n,
                               "new_vertices": len(output)-original_n,
                               "center_pairs": pair_counts,
                               "approximate_pairs_rejected": rejected_approximate_pairs}


def core_indices(graph, degree=6):
    neighbors = [set() for _ in graph['vertices']]
    for u, v in graph['edges']:
        neighbors[u].add(v)
        neighbors[v].add(u)
    queue = deque(i for i, ns in enumerate(neighbors) if len(ns) < degree)
    removed = set()
    while queue:
        u = queue.popleft()
        if u in removed:
            continue
        removed.add(u)
        for v in neighbors[u]:
            neighbors[v].discard(u)
            if len(neighbors[v]) < degree:
                queue.append(v)
    return set(range(len(neighbors))) - removed


def run(seed_name, out, rounds=1, timeout=60, max_vertices=30000,
        separations=(Fraction(4, 3), Fraction(8, 3))):
    if seed_name not in ('v1939', 'heule') or rounds not in (1, 2):
        raise ValueError("Unsupported seed or number of rounds")
    separations = tuple(Fraction(s) for s in separations)
    if any(s not in SEPARATIONS for s in separations):
        raise ValueError("Unsupported separation")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).read_bytes()
    (out/'program.py').write_bytes(source)
    started = time.monotonic()
    points = list(v1939_points() if seed_name == 'v1939' else heule_union_points())
    seed_set = set(points)
    metadata = {'family': 'exact-common-unit-neighbor-closure', 'seed': seed_name,
                'requested_rounds': rounds, 'separations_squared': list(map(str, separations)),
                'max_vertices': max_vertices, 'solver_timeout_seconds': min(timeout, 60),
                'success_certificate': False,
                'program_sha256': hashlib.sha256(source).hexdigest()}
    (out/'configuration.json').write_text(json.dumps(metadata, indent=2)+'\n')
    round_reports, provenance = [], []
    for r in range(1, rounds+1):
        points, additions, report = intersect_round(points, separations, max_vertices)
        for addition in additions:
            addition['round'] = r
        provenance.extend(additions)
        round_reports.append(report)
        print(json.dumps({'stage': 'intersection_round', 'round': r, **report}), flush=True)
        if report['truncated_at_vertex_cap'] or report['new_vertices'] == 0:
            break
    metadata['rounds'] = round_reports
    graph = export(points, out/'graph.json', metadata)
    canonical = canonical_points(points)
    canonical_index = {p: j for j, p in enumerate(canonical)}
    remapping = [canonical_index[p] for p in points]
    provenance_object = {
        'schema': 'hn-circle-intersection-provenance-v1',
        'graph_sha256': hashlib.sha256((out/'graph.json').read_bytes()).hexdigest(),
        'seed_count': len(seed_set), 'generation_to_graph_index': remapping,
        'new_points': provenance,
    }
    (out/'provenance.json').write_text(json.dumps(provenance_object, separators=(',', ':'))+'\n')
    core = core_indices(graph, 6)
    new_indices = {j for j, p in enumerate(canonical) if p not in seed_set}
    report = {
        'schema': 'hn-circle-intersection-run-v1',
        'vertices': len(canonical), 'edges': len(graph['edges']),
        'new_vertices': len(new_indices), 'six_core_vertices': len(core),
        'new_vertices_in_six_core': len(core & new_indices),
        'six_core_edges': sum(u in core and v in core for u, v in graph['edges']),
        'edges_touching_new_vertices': sum(u in new_indices or v in new_indices
                                          for u, v in graph['edges']),
        'rounds': round_reports,
        'success_certificate': False,
    }
    print(json.dumps({'stage': 'core_analysis', **report}), flush=True)
    statuses = {}
    for k in (6, 5):
        result = solve(out/'graph.json', k, out/f'color-{k}', min(timeout, 60))
        statuses[str(k)] = result['status']
    report.update({'coloring_status': statuses, 'elapsed_seconds': time.monotonic()-started})
    (out/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', choices=('v1939', 'heule'), required=True)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--rounds', type=int, choices=(1, 2), default=1)
    parser.add_argument('--timeout', type=float, default=60)
    parser.add_argument('--max-vertices', type=int, default=30000)
    parser.add_argument('--separations', nargs='+', default=['4/3', '8/3'])
    args = parser.parse_args()
    if args.timeout <= 0 or args.max_vertices <= 0:
        parser.error('Timeout and vertex cap must be positive')
    run(args.seed, args.out, args.rounds, args.timeout, args.max_vertices,
        tuple(map(Fraction, args.separations)))


if __name__ == '__main__':
    main()
