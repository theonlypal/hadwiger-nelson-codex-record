"""Exact Minkowski sums coupling Haugland seeds to Heule's V31 seed.

Coordinates lie in Q(zeta_42)(sqrt5,sqrt33), with the schema documented
in hn.cyclotomic_tower. No chromatic conclusion is inferred from the
construction or from numerical candidate enumeration.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time

from .cyclotomic_tower import Element, FIELD_METADATA, ONE, ZERO, ZETA, SQRT5, SQRT33

B = 2*ZETA**7-ONE  # i*sqrt(3) in the selected embedding
HEX_ROTATION = ZETA**7
HALF_THETA3 = (SQRT33+B)/6
THETA3 = (5+B*SQRT33/3)/6
THETA4 = (7+B*SQRT5)/8


def canonical_points(points):
    return tuple(sorted(set(points), key=lambda z: (z.denominator, z.coefficients)))


@lru_cache(maxsize=1)
def v31_points():
    for rotation in (HEX_ROTATION, HALF_THETA3, THETA3, THETA4):
        if not rotation.is_unit():
            raise ArithmeticError('A Heule rotation failed the exact unit identity')
    if HALF_THETA3**2 != THETA3:
        raise ArithmeticError('Heule half-angle identity failed')
    directions = canonical_points(HEX_ROTATION**j*HALF_THETA3**k
                                  for j in range(6) for k in range(-2, 3))
    if len(directions) != 30 or not all(z.is_unit() for z in directions):
        raise ArithmeticError('V31 unit-direction construction failed')
    return (ZERO,)+directions


@lru_cache(maxsize=2)
def lift_haugland(seed='g1'):
    from .construct import load_g1_points, g3_points
    if seed not in ('g1', 'g3'):
        raise ValueError('Supported Haugland seeds are g1 and g3')
    original = load_g1_points() if seed == 'g1' else g3_points()
    return canonical_points(Element.from_legacy(z) for z in original)


def mixed_points(seed='g1', max_vertices=66000):
    points = canonical_points(z+w for z in lift_haugland(seed) for w in v31_points())
    if len(points) > max_vertices:
        raise ValueError(f'Minkowski sum has {len(points)} vertices, exceeding cap {max_vertices}; no truncation performed')
    return points


def graph_from_points(points, *, directions=None, metadata=None):
    from .search import all_candidate_edges, known_direction_edges
    points = canonical_points(points)
    started = time.monotonic()
    edges = all_candidate_edges(points) if directions is None else known_direction_edges(points, directions)
    return {'schema': 'hn-cyclotomic-tower-v1', 'field': dict(FIELD_METADATA),
            'vertices': [z.to_json() for z in points], 'edges': edges,
            'metadata': dict(metadata or {},
                             edge_rule='Exact norm equality; candidate enumeration may omit unit pairs',
                             induced_graph_claimed=False,
                             edge_enumeration_seconds=time.monotonic()-started)}


def export(points, out, metadata=None, directions=None):
    graph = graph_from_points(points, directions=directions, metadata=metadata)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(graph, separators=(',', ':'))+'\n')
    print(json.dumps({'stage': 'graph_written', 'path': str(out),
                      'vertices': len(graph['vertices']), 'edges': len(graph['edges'])}), flush=True)
    return graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', choices=('g1', 'g3'), default='g1')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-vertices', type=int, default=66000)
    parser.add_argument('--lift-only', action='store_true')
    parser.add_argument('--v31-only', action='store_true')
    args = parser.parse_args()
    if args.max_vertices <= 0:
        parser.error('Vertex cap must be positive')
    if args.lift_only and args.v31_only:
        parser.error('Choose at most one control construction')
    started = time.monotonic()
    if args.v31_only:
        points = v31_points()
    elif args.lift_only:
        points = lift_haugland(args.seed)
    else:
        points = mixed_points(args.seed, args.max_vertices)
    from .construct import SOURCE_PATH
    metadata = {
        'family': 'heule-v31-control' if args.v31_only else 'haugland-lift-control' if args.lift_only else 'haugland-minkowski-heule-v31',
        'haugland_seed': args.seed, 'haugland_seed_vertices': len(lift_haugland(args.seed)),
        'v31_vertices': len(v31_points()), 'vertex_cap': args.max_vertices,
        'haugland_source_sha256': hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest(),
        'haugland_source': 'https://arxiv.org/html/2608.04542v1',
        'heule_source': 'https://arxiv.org/html/1805.12181v1',
        'arithmetic_program_sha256': hashlib.sha256(Path(__file__).with_name('cyclotomic_tower.py').read_bytes()).hexdigest(),
        'construction_program_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'point_generation_seconds': time.monotonic()-started,
        'success_certificate': False,
    }
    print(json.dumps({'stage': 'edge_enumeration', 'vertices': len(points), 'metadata': metadata}), flush=True)
    export(points, args.output, metadata)


if __name__ == '__main__':
    main()
