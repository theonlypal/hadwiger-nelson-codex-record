"""Combine two exact Voronov rotations with the same 32,257-point base.

Each displayed unit edge is checked exactly. The two radical generators
have independent square classes, and their joint field escapes the known
F25 ring map. Neither observation asserts any chromatic lower bound.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from .voronov_base import m32257_points
from .voronov_field import Point, lift_base_point, PSI1, PSI7
from .search import all_candidate_edges
from .solve_external import solve_external

FIELD = {'base_radicands': [2, 3],
         'extension_squares': [[5, 3, -1, 3], [5, 3, 1, -3]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=300)
    parser.add_argument('--rotations', choices=['1', '7', '1,7'], default='1,7')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    source = Path(__file__).read_bytes()
    (args.out/'program.py').write_bytes(source)
    base = tuple(lift_base_point(p) for p in m32257_points())
    rotations = {'1': PSI1, '7': PSI7}
    copies = [set(base)] + [set(rotations[j]*p for p in base) for j in args.rotations.split(',')]
    points = sorted(set.union(*copies), key=lambda p: (p.x.denominator, p.x.coefficients,
                                                     p.y.denominator, p.y.coefficients))
    memberships = [sum(1 << j for j, copy in enumerate(copies) if p in copy) for p in points]
    metadata = {'family': 'voronov-joint-independent-radical-rotations',
                'source': 'https://arxiv.org/html/2106.11824v4',
                'source_repository': 'https://github.com/vsvor/dist-graphs',
                'source_revision': 'e3714d0a156f6ed151d4521debb2b89ce4f1075c',
                'rotations': args.rotations.split(','), 'base_vertices': len(base),
                'copy_sizes': [len(c) for c in copies], 'distinct_vertices': len(points),
                'program_sha256': hashlib.sha256(source).hexdigest(),
                'field_program_sha256': hashlib.sha256(Path(__file__).with_name('voronov_field.py').read_bytes()).hexdigest(),
                'point_generation_seconds': time.monotonic()-started,
                'induced_graph_claimed': False, 'success_certificate': False}
    print(json.dumps({'stage': 'points_ready', **metadata}), flush=True)
    edges = all_candidate_edges(points)
    strict = [edge for edge in edges if not (memberships[edge[0]] & memberships[edge[1]])]
    metadata.update(edges=len(edges), strict_cross_copy_edges=len(strict),
                    edge_rule='Numerical candidate screen; exact algebraic norm equality accepts every edge',
                    construction_seconds=time.monotonic()-started)
    graph = {'schema': 'hn-voronov-field-v1', 'field': FIELD,
             'vertices': [p.to_json() for p in points], 'edges': edges, 'metadata': metadata}
    graph_path = args.out/'graph.json'
    graph_path.write_text(json.dumps(graph, separators=(',', ':'))+'\n')
    (args.out/'membership.json').write_text(json.dumps(memberships)+'\n')
    (args.out/'cross-edges.json').write_text(json.dumps(strict)+'\n')
    print(json.dumps({'stage': 'graph_written', **metadata}), flush=True)
    solve_external(graph_path, args.out/'six-color', args.timeout, configuration='sat')


if __name__ == '__main__':
    main()
