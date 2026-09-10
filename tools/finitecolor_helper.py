#!/usr/bin/env python3
"""Exact finite-field quotient coloring screen; never a Euclidean certificate.

Fp² is represented by a+b*w with w²=d, as a+p*b. The graph on Fp^4
joins differences whose two field-coordinate squares sum to one.
Rank-two Fp-linear quotient maps are enumerated once per row space, via
reduced row echelon form. A loop rejects the quotient; SAT models lift.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import threading
import time

import numpy as np
from pysat.solvers import Solver


def field_mul(a, b, p, d):
    return ((a % p)*(b % p)+d*(a//p)*(b//p)) % p + p*((a % p)*(b//p)+(a//p)*(b % p)) % (p*p)


def unit_directions(p, d):
    assert pow(d, (p-1)//2, p) == p-1
    q = p*p
    result = []
    for x, y in product(range(q), repeat=2):
        a, b = field_mul(x, x, p, d), field_mul(y, y, p, d)
        if ((a % p+b % p) % p, (a//p+b//p) % p) == (1, 0):
            result.append((x % p, x//p, y % p, y//p))
    assert len(result) == q-1
    return result


def rref_maps(p, rank):
    for pivots in combinations(range(4), rank):
        free = [(r, c) for r in range(rank) for c in range(pivots[r]+1, 4)
                if c not in pivots]
        for values in product(range(p), repeat=len(free)):
            rows = [[int(c == pivots[r]) for c in range(4)] for r in range(rank)]
            for (r, c), value in zip(free, values):
                rows[r][c] = value
            yield rows


def quotient_directions(rows, directions, p):
    return {tuple(sum(x*y for x, y in zip(row, s)) % p for row in rows)
            for s in directions}


def edges_from_directions(directions, p, rank=2):
    vertices = list(product(range(p), repeat=rank))
    index = {v: i for i, v in enumerate(vertices)}
    edges = set()
    for i, v in enumerate(vertices):
        for s in directions:
            j = index[tuple((x+y) % p for x, y in zip(v, s))]
            if i != j:
                edges.add((min(i, j), max(i, j)))
    return vertices, sorted(edges)


def color(n, edges, k, timeout):
    clauses = []
    for v in range(n):
        clauses.append([v*k+c+1 for c in range(k)])
        clauses.extend([-(v*k+c+1), -(v*k+d+1)]
                       for c in range(k) for d in range(c+1, k))
    clauses.extend([-(u*k+c+1), -(v*k+c+1)] for u, v in edges for c in range(k))
    clauses.append([1])
    with Solver(name='g4', bootstrap_with=clauses) as solver:
        timer = threading.Timer(timeout, solver.interrupt)
        timer.start()
        try:
            result = solver.solve_limited(expect_interrupt=True)
            if result is not True:
                return {'status': 'UNSAT_UNVERIFIED' if result is False else 'UNKNOWN_TIMEOUT'}
            positive = {v for v in solver.get_model() if v > 0}
            colors = [next(c for c in range(k) if v*k+c+1 in positive) for v in range(n)]
            assert all(colors[u] != colors[v] for u, v in edges)
            return {'status': 'SAT_VERIFIED', 'colors': colors}
        finally:
            timer.cancel()


def spectrum_counts(directions, p):
    """Exact coefficients of eigenvalues in Z[zeta_p], by additive characters."""
    multiplicities = Counter()
    representatives = {}
    for a in product(range(p), repeat=4):
        counts = [0]*p
        for s in directions:
            counts[sum(x*y for x, y in zip(a, s)) % p] += 1
        key = tuple(counts)
        multiplicities[key] += 1
        representatives.setdefault(key, a)
    result = []
    for counts, count in multiplicities.items():
        approx = sum(n*np.cos(2*np.pi*j/p) for j, n in enumerate(counts))
        result.append({'counts': counts, 'multiplicity': count,
                       'representative_character': representatives[counts],
                       'approximate_eigenvalue_only': float(approx)})
    return sorted(result, key=lambda x: x['approximate_eigenvalue_only'])


def run(p, d, out, timeout):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).read_bytes()
    (out/'program.py').write_bytes(source)
    started = time.monotonic()
    directions = unit_directions(p, d)
    spectrum = spectrum_counts(directions, p)
    (out/'spectrum.json').write_text(json.dumps(spectrum, indent=2)+'\n')
    rank1 = [rows for rows in rref_maps(p, 1)
             if (0,) not in quotient_directions(rows, directions, p)]
    quotient_types = {}
    map_counts = Counter()
    for rows in rref_maps(p, 2):
        qdirs = quotient_directions(rows, directions, p)
        if (0, 0) in qdirs:
            map_counts['looped'] += 1
            continue
        map_counts['loopless'] += 1
        key = tuple(sorted(qdirs))
        quotient_types.setdefault(key, rows)
    report = {'schema': 'hn-finite-field-quotient-screen-v1',
              'prime': p, 'extension_w_squared': d,
              'purpose': 'finite-field coloring obstruction screen; not a Euclidean realization',
              'program_sha256': hashlib.sha256(source).hexdigest(),
              'unit_directions': directions, 'rank1_loopless_maps': rank1,
              'rank2_map_counts': dict(map_counts), 'quotient_direction_types': len(quotient_types),
              'coloring_status': {}, 'success_certificate': False}
    best = []
    for qdirs, rows in sorted(quotient_types.items(), key=lambda item: len(item[0])):
        vertices, edges = edges_from_directions(qdirs, p)
        answer = color(p*p, edges, 6, timeout)
        record = {'rows': rows, 'degree': len(qdirs), **answer}
        best.append(record)
        report['coloring_status'][answer['status']] = report['coloring_status'].get(answer['status'], 0)+1
        if answer['status'] == 'SAT_VERIFIED':
            colors = answer['colors']
            table = {v: colors[j] for j, v in enumerate(vertices)}
            full = []
            for v in product(range(p), repeat=4):
                image = tuple(sum(x*y for x, y in zip(row, v)) % p for row in rows)
                full.append(table[image])
            full_index = {v: j for j, v in enumerate(product(range(p), repeat=4))}
            for v, i in full_index.items():
                for s in directions:
                    w = tuple((x+y) % p for x, y in zip(v, s))
                    assert full[i] != full[full_index[w]]
            report['verified_lifted_coloring'] = {'rows': rows, 'colors_lexicographic_Fp4': full}
            print(json.dumps({'prime': p, 'status': 'SAT_VERIFIED', 'quotient_degree': len(qdirs)}), flush=True)
            break
    report['quotient_results'] = best
    report['elapsed_seconds'] = time.monotonic()-started
    (out/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('unit_directions', 'quotient_results', 'verified_lifted_coloring')}), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, choices=(5, 7), required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=.25)
    args = parser.parse_args()
    run(args.prime, 2 if args.prime == 5 else 3, args.out, args.timeout)
