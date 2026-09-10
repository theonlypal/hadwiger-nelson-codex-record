"""Finite candidates adding multiplicative powers of Haugland's odd unit.

The extra direction is u1^2, still an exact norm-one field element. The
direction-only modulo-five filter is not a graph-colorability certificate.
Every exported candidate receives its own unrestricted target SAT test.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

from .construct import unit_directions, g2_points, g3_points, spindle_rotation
from .field import ZERO, ONE, ZETA
from .search import export
from .solve_external import solve_external


def directions(powers=3, spindle_layers=1):
    odd = unit_directions()[1]
    rho = spindle_rotation()
    result = {odd**k * ZETA**j * rho**l
              for k in range(powers) for j in range(42)
              for l in range(spindle_layers)}
    if not all(u.is_unit() for u in result):
        raise ArithmeticError('Non-unit direction')
    return tuple(sorted(result, key=lambda z: (z.denominator, z.coefficients)))


def endpoint_points(pool, max_paths=1000000):
    """All vertices on <=6-step walks from 0 to i*sqrt3 in this pool.

    Unordered triples suffice: every split of a <=6-step multiset appears
    in the two tables. All subset sums recover every ordering of every
    accepted multiset, so this construction is exact and complete for the
    specified finite walk template.
    """
    table = {ZERO: [()]}
    for i, a in enumerate(pool):
        table.setdefault(a, []).append((i,))
        for j in range(i, len(pool)):
            ab = a+pool[j]
            table.setdefault(ab, []).append((i, j))
            for k in range(j, len(pool)):
                table.setdefault(ab+pool[k], []).append((i, j, k))
    print(json.dumps({'stage': 'half_path_table', 'directions': len(pool),
                      'distinct_half_sums': len(table)}), flush=True)
    target = 2*ZETA**7-ONE
    paths = set()
    for position, lefts in table.items():
        for left in lefts:
            for right in table.get(target-position, ()):
                paths.add(tuple(sorted(left+right)))
                if len(paths) > max_paths:
                    raise ValueError('Endpoint-path cap exceeded; no truncated construction exported')
    print(json.dumps({'stage': 'endpoint_multisets', 'count': len(paths)}), flush=True)
    points = {ZERO}
    for indices in sorted(paths):
        subsets = {ZERO}
        endpoint = ZERO
        for j in indices:
            endpoint += pool[j]
            subsets |= {z+pool[j] for z in subsets}
        if endpoint != target:
            raise ArithmeticError('Wrong exact endpoint in accepted path')
        points.update(subsets)
    print(json.dumps({'stage': 'endpoint_path_vertices', 'vertices': len(points)}), flush=True)
    return points


def construct(family, powers=3, spindle_layers=1, max_vertices=100000):
    pool = directions(powers, spindle_layers)
    if family == 'radius-two':
        points = {ZERO, *pool}
        for a, b in itertools.combinations_with_replacement(pool, 2):
            points.add(a+b)
    elif family == 'g3-minkowski':
        odd = unit_directions()[1]
        offsets = (ZERO, *(odd**(powers-1)*ZETA**j for j in range(42)))
        points = {z+w for z in g3_points() for w in offsets}
    elif family in ('t6', 't6-spindle'):
        points = endpoint_points(pool)
        if family == 't6-spindle':
            points = set(g3_points(g2_points(points)))
    else:
        raise ValueError('Unknown family')
    if len(points) > max_vertices:
        raise ValueError(f'{len(points)} vertices exceed cap {max_vertices}; no truncation performed')
    return points, pool


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family', choices=['radius-two', 'g3-minkowski', 't6', 't6-spindle'], required=True)
    parser.add_argument('--powers', type=int, default=3)
    parser.add_argument('--spindle-layers', type=int, choices=[1, 2], default=1)
    parser.add_argument('--max-vertices', type=int, default=100000)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--diagnostic-five', action='store_true')
    args = parser.parse_args()
    if not 2 <= args.powers <= 4:
        parser.error('Supported powers count is 2 through 4')
    args.out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).read_bytes()
    (args.out/'program.py').write_bytes(source)
    started = time.monotonic()
    points, pool = construct(args.family, args.powers, args.spindle_layers, args.max_vertices)
    metadata = dict(vars(args), out=str(args.out), family='odd-power-'+args.family,
                    exact_direction_count=len(pool), point_generation_seconds=time.monotonic()-started,
                    direction_rule='u1^k*zeta42^j*rho^l; u1 is Haugland unit_directions()[1]',
                    program_sha256=hashlib.sha256(source).hexdigest(), success_certificate=False)
    (args.out/'configuration.json').write_text(json.dumps(metadata, indent=2)+'\n')
    export(points, args.out/'graph.json', metadata)
    result = solve_external(args.out/'graph.json', args.out/'six-color', args.timeout, configuration='sat')
    if args.diagnostic_five and result['status'] == 'SAT':
        solve_external(args.out/'graph.json', args.out/'five-color', args.timeout, configuration='sat', k=5)


if __name__ == '__main__':
    main()
