"""Exact L10,2 point sets from Voronov--Neopryatnaya--Dergachev.

Primary source: arXiv:2106.11824v4, sections 5--6, series 2. This module
constructs points only. Its base field has the previously established
five-color ring map, so no target SAT run is spent on this base alone.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path

from .multiquadratic import Point, SQRT2, SQRT3, ONE, ORIGIN
from .mq_construct import canonical_points

PHI0 = Point((SQRT2*SQRT3+SQRT2)/4, (SQRT2*SQRT3-SQRT2)/4)
PHI1 = Point(SQRT2*SQRT3/3, SQRT3/3)


@lru_cache(maxsize=1)
def m73_points():
    if not PHI0.is_unit() or not PHI1.is_unit() or PHI0**24 != Point(ONE):
        raise ArithmeticError('Invalid exact generators')
    directions = canonical_points(PHI0**j*PHI1**k for j in range(24) for k in (-1, 0, 1))
    if len(directions) != 72 or not all(p.is_unit() for p in directions):
        raise ArithmeticError('Expected 72 distinct exact unit directions')
    return (ORIGIN,)+directions


@lru_cache(maxsize=1)
def m2_points():
    candidates = set(a+b for a in m73_points() for b in m73_points())
    return canonical_points(p for p in candidates if p.norm_squared() <= ONE)


@lru_cache(maxsize=1)
def m32257_points():
    points = canonical_points(a+b for a in m2_points() for b in m73_points())
    if len(points) != 32257:
        raise ArithmeticError(f'Primary-source base count mismatch: {len(points)}')
    return points


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Output already exists; preserving previous artifact')
    points = m32257_points()
    record = {'schema': 'hn-multiquadratic-point-set-v1', 'radicands': [2, 3, 5, 11],
              'source': 'https://arxiv.org/html/2106.11824v4',
              'construction': 'M1 order24 exponents -1..1; M2=clip(M1+M1,1); M3=M2+M1',
              'stage_sizes': [len(m73_points()), len(m2_points()), len(points)],
              'program_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'vertices': [p.to_json() for p in points],
              'goal_satisfied': False, 'success_certificate': False}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, separators=(',', ':'))+'\n')
    print(json.dumps({'stage_sizes': record['stage_sizes'], 'output': str(args.out),
                      'sha256': hashlib.sha256(args.out.read_bytes()).hexdigest()}), flush=True)
