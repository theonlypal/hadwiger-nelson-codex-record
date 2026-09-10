"""Try an independent sixth-color class repairing additive modulo-five maps.

Variables s_v say which vertices receive color six. Each monochromatic
edge requires s_u or s_v; every graph edge requires not(s_u and s_v).
This is a 2-SAT ansatz, not the unrestricted six-coloring problem.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from tools import verify


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


@lru_cache(None)
def residue(value):
    q = Fraction(value)
    if q.denominator % 5 == 0:
        raise ValueError('A coefficient has no residue modulo five')
    return q.numerator*pow(q.denominator, -1, 5) % 5


def odd_cycle(n, bad):
    adjacency = [[] for _ in range(n)]
    for u,v in bad:
        adjacency[u].append(v); adjacency[v].append(u)
    parity = [-1]*n; parent = [-1]*n
    for root in range(n):
        if parity[root] >= 0 or not adjacency[root]:
            continue
        parity[root] = 0; queue = deque([root])
        while queue:
            u = queue.popleft()
            for v in adjacency[u]:
                if parity[v] < 0:
                    parent[v] = u; parity[v] = parity[u]^1; queue.append(v)
                elif parity[u] == parity[v]:
                    left = []; w = u
                    while w >= 0:
                        left.append(w); w = parent[w]
                    positions = {w:i for i,w in enumerate(left)}
                    right = []; w = v
                    while w not in positions:
                        right.append(w); w = parent[w]
                    cycle = left[:positions[w]+1] + list(reversed(right)) + [u]
                    return cycle, Counter(map(len, adjacency))
    return None, Counter(map(len, adjacency))


def implication_path(matrix, labels, source, target):
    queue = deque([source]); previous = {source:None}
    while queue:
        u = queue.popleft()
        if u == target:
            path = []
            while u is not None:
                path.append(u); u = previous[u]
            return list(reversed(path))
        for v in matrix.indices[matrix.indptr[u]:matrix.indptr[u+1]]:
            v = int(v)
            if labels[v] == labels[source] and v not in previous:
                previous[v] = u; queue.append(v)
    raise ArithmeticError('Missing path inside a strongly connected component')


def literal(index):
    return (index//2+1) * (-1 if index%2 else 1)


def repair(n, edges, colors):
    bad = [(u,v) for u,v in edges if colors[u] == colors[v]]
    base = np.asarray(edges, dtype=np.int32)
    extra = np.asarray(bad, dtype=np.int32).reshape((-1,2))
    rows = np.concatenate((2*base[:,0], 2*base[:,1], 2*extra[:,0]+1, 2*extra[:,1]+1))
    cols = np.concatenate((2*base[:,1]+1, 2*base[:,0]+1, 2*extra[:,1], 2*extra[:,0]))
    matrix = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows,cols)), shape=(2*n,2*n)).tocsr()
    component_count, labels = connected_components(matrix, directed=True, connection='strong')
    contradictory = np.flatnonzero(labels[0::2] == labels[1::2])
    cycle, degrees = odd_cycle(n, bad)
    report = {'monochromatic_edges':len(bad), 'monochromatic_degree_histogram':dict(degrees),
              'two_sat_variables':n, 'two_sat_clauses':len(edges)+len(bad),
              'strong_components':int(component_count),
              'contradictory_variables':len(contradictory),
              'scope':'A failed repair excludes only this additive map followed by one independent sixth-color class.'}
    if len(contradictory):
        report['status'] = 'REPAIR_ANSATZ_UNSAT'
        if cycle:
            bad_set = {tuple(sorted(edge)) for edge in bad}
            if (len(cycle)-1)%2 != 1 or len(set(cycle[:-1])) != len(cycle)-1:
                raise ArithmeticError('Not a simple odd-cycle obstruction')
            if any(tuple(sorted(edge)) not in bad_set for edge in zip(cycle,cycle[1:])):
                raise ArithmeticError('Odd-cycle obstruction uses an edge outside the bad graph')
            report.update(odd_cycle=cycle, odd_cycle_length=len(cycle)-1,
                          obstruction_checked=True,
                          obstruction_argument='Each cycle edge requires exactly one endpoint in S. Alternation around this odd cycle is impossible.')
        else:
            vertex = int(contradictory[0]); source = 2*vertex
            paths = [implication_path(matrix, labels, source, source+1),
                     implication_path(matrix, labels, source+1, source)]
            edge_set = {tuple(sorted(edge)) for edge in edges}
            bad_set = {tuple(sorted(edge)) for edge in bad}
            for path in paths:
                for u,v in zip(path,path[1:]):
                    pair = tuple(sorted((u//2,v//2)))
                    if u%2 == 0 and v%2 == 1:
                        if pair not in edge_set: raise ArithmeticError('Invalid independence implication')
                    elif u%2 == 1 and v%2 == 0:
                        if pair not in bad_set: raise ArithmeticError('Invalid cover implication')
                    else:
                        raise ArithmeticError('Invalid implication polarity')
            report.update(contradiction_vertex=vertex,
                          implication_paths=[[literal(i) for i in path] for path in paths],
                          obstruction_checked=True)
        return report, None, bad
    if cycle:
        raise ArithmeticError('A satisfiable repair cannot have a monochromatic odd cycle')
    # CSR SCC component numbering is normally reverse topological. Check
    # both orders explicitly rather than trusting an undocumented order.
    for selected in (labels[0::2] < labels[1::2], labels[0::2] > labels[1::2]):
        if (all(not(selected[u] and selected[v]) for u,v in edges)
                and all(selected[u] or selected[v] for u,v in bad)):
            repaired = [5 if selected[v] else colors[v] for v in range(n)]
            report.update(status='SIX_COLOR_REPAIR_FOUND', sixth_color_vertices=int(sum(selected)))
            return report, repaired, bad
    # A future SciPy implementation may number components differently.
    from pysat.solvers import Solver
    with Solver(name='g4', bootstrap_with=([-(u+1),-(v+1)] for u,v in edges)) as solver:
        for u,v in bad: solver.add_clause([u+1,v+1])
        if not solver.solve(): raise ArithmeticError('SAT solver and SCC result disagree')
        selected = {i-1 for i in solver.get_model() if i>0}
    repaired = [5 if v in selected else colors[v] for v in range(n)]
    report.update(status='SIX_COLOR_REPAIR_FOUND', sixth_color_vertices=len(selected))
    return report, repaired, bad


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph', type=Path, default=Path('runs/odd-power-t6-spindle-v1/graph.json'))
    parser.add_argument('--first-map', type=Path, default=Path('runs/mixed-haugland-g3-v31-linear5/result.json'))
    parser.add_argument('--maps', type=Path, default=Path('runs/mixed-five-forced-pair-screen/functional-sampling.json'))
    parser.add_argument('--out', type=Path, default=Path('runs/odd-power-t6-sixth-repair-v1'))
    parser.add_argument('--max-maps', type=int, default=20)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'program.py').write_bytes(Path(__file__).read_bytes())
    raw = args.graph.read_bytes(); graph = json.loads(raw)
    if graph['schema'] != 'hn-exact-v1' or graph['field']['cyclotomic_order'] != 42 or graph['field']['radicand'] != 5:
        raise ValueError('Expected the legacy zeta42/sqrt5 coefficient basis')
    n = len(graph['vertices']); edges = graph['edges']
    residues = np.asarray([[residue(c) for c in v['a']+v['b']] for v in graph['vertices']], dtype=np.int32)
    input_maps = [json.loads(args.first_map.read_text())['linear_functional']]
    input_maps.extend(json.loads(args.maps.read_text())['functionals'])
    maps = []; seen = set()
    for values in input_maps:
        functional = tuple(values[:24])
        if len(functional)!=24 or any(type(v) is not int or not 0<=v<5 for v in functional):
            raise ValueError('Invalid additive functional')
        if functional not in seen:
            maps.append(functional); seen.add(functional)
    maps = maps[:args.max_maps]
    config = {'graph':str(args.graph.resolve()), 'graph_sha256':hashlib.sha256(raw).hexdigest(),
              'first_map_file_sha256':verify._sha256(args.first_map),
              'map_file_sha256':verify._sha256(args.maps), 'maps_tested_limit':len(maps),
              'coordinate_basis':'a_j*zeta42^j + b_j*sqrt5*zeta42^j, j=0..11',
              'vertices':n, 'edges':len(edges), 'goal_satisfied':False,
              'source_program_sha256':verify._sha256(Path(__file__)),
              'formula':'For every edge uv: (!s_u OR !s_v). For each L-monochromatic edge uv: (s_u OR s_v).'}
    dump(args.out/'config.json', config); records = []
    for number, functional in enumerate(maps):
        start = time.monotonic()
        colors = list(map(int, residues @ np.asarray(functional, dtype=np.int32) % 5))
        report, repaired, bad = repair(n, edges, colors)
        cnf = (f'p cnf {n} {len(edges)+len(bad)}\n'+
               ''.join(f'-{u+1} -{v+1} 0\n' for u,v in edges)+
               ''.join(f'{u+1} {v+1} 0\n' for u,v in bad)).encode()
        compressed = gzip.compress(cnf, compresslevel=6, mtime=0)
        if gzip.decompress(compressed)!=cnf: raise ArithmeticError('CNF compression mismatch')
        (args.out/f'map-{number:02}.cnf.gz').write_bytes(compressed)
        report.update(map_index=number, linear_functional=list(functional),
                      tentative_color_distribution=dict(Counter(colors)),
                      cnf_sha256=hashlib.sha256(cnf).hexdigest(),
                      cnf_gzip_sha256=hashlib.sha256(compressed).hexdigest(),
                      elapsed_seconds=time.monotonic()-start, goal_satisfied=False)
        dump(args.out/f'map-{number:02}.json', report); records.append(report)
        print(json.dumps({k:v for k,v in report.items() if k not in ('linear_functional','implication_paths','odd_cycle','monochromatic_degree_histogram')}), flush=True)
        if repaired is not None:
            witness = {'k':6, 'colors':repaired, 'graph_sha256':config['graph_sha256']}
            verify.verify_coloring(graph, witness, 6)
            geometry = verify.verify_geometry(graph)
            if args.graph.read_bytes()!=raw: raise ValueError('Graph changed during repair')
            dump(args.out/'coloring.json', witness)
            dump(args.out/'verification.json', {'graph_sha256':config['graph_sha256'],
                  'geometry':geometry.__dict__, 'proper_six_coloring_verified':True,
                  'coloring_sha256':verify._sha256(args.out/'coloring.json'), 'goal_satisfied':False})
            dump(args.out/'result.json', dict(config, status='SIX_COLORING_INDEPENDENTLY_VERIFIED', attempts=records))
            return
        dump(args.out/'result.json', dict(config, status='NO_REPAIR_IN_TESTED_MAPS', attempts=records))


if __name__ == '__main__':
    main()
