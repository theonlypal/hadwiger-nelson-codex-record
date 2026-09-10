"""Upper-color filters using old-direction cosets modulo seven and a mod-five map.

A quotient coloring pulls back to the exact graph. Quotient loops, UNSAT,
and timeouts establish no lower bound for the original graph.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
from threading import Timer
import time

import numpy as np
from pysat.solvers import Solver

from hn.construct import unit_directions, spindle_rotation
from hn.solve import clauses, cnf_bytes
from tools import verify


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


@lru_cache(None)
def residue(value, prime, multiplier=1):
    q = multiplier*Fraction(value)
    if q.denominator % prime == 0:
        raise ValueError(f'Coordinate has a pole modulo {prime} after scaling')
    return q.numerator*pow(q.denominator, -1, prime) % prime


def nullspace(rows, prime):
    reduced = np.array(rows, dtype=np.int32); width = reduced.shape[1]
    rank = 0; pivots = []
    for column in range(width):
        choices = np.flatnonzero(reduced[rank:,column])
        if not len(choices): continue
        i = rank+choices[0]; reduced[[rank,i]] = reduced[[i,rank]]
        reduced[rank] = reduced[rank]*pow(int(reduced[rank,column]), -1, prime) % prime
        for i in range(len(reduced)):
            if i != rank and reduced[i,column]:
                reduced[i] = (reduced[i]-reduced[i,column]*reduced[rank]) % prime
        pivots.append(column); rank += 1
    free = [column for column in range(width) if column not in pivots]
    annihilator = np.zeros((len(free),width), dtype=np.int32)
    for i,column in enumerate(free):
        annihilator[i,column] = 1
        for row,pivot in enumerate(pivots):
            annihilator[i,pivot] = -reduced[row,column] % prime
    if np.any(np.array(rows)@annihilator.T % prime):
        raise ArithmeticError('The proposed annihilator does not kill the old directions')
    return reduced[:rank], annihilator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph', type=Path, default=Path('runs/odd-power-t6-spindle-v1/graph.json'))
    parser.add_argument('--first-map', type=Path, default=Path('runs/mixed-haugland-g3-v31-linear5/result.json'))
    parser.add_argument('--maps', type=Path, default=Path('runs/mixed-five-forced-pair-screen/functional-sampling.json'))
    parser.add_argument('--out', type=Path, default=Path('runs/odd-power-t6-coset-repair-v1'))
    parser.add_argument('--max-maps', type=int, default=20)
    parser.add_argument('--timeout-per-map', type=float, default=5)
    parser.add_argument('--split-loop-endpoints', action='store_true',
                        help='Give each endpoint of an initial quotient loop its own type')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'program.py').write_bytes(Path(__file__).read_bytes())
    raw = args.graph.read_bytes(); graph = json.loads(raw)
    if graph['schema'] != 'hn-exact-v1' or graph['field']['cyclotomic_order'] != 42 or graph['field']['radicand'] != 5:
        raise ValueError('Expected the zeta42/sqrt5 legacy basis')
    n = len(graph['vertices']); edges = graph['edges']
    points7 = np.asarray([[residue(c,7,7) for c in v['a']+v['b']] for v in graph['vertices']], dtype=np.int32)
    points5 = np.asarray([[residue(c,5) for c in v['a']+v['b']] for v in graph['vertices']], dtype=np.int32)
    rho = spindle_rotation(); old = unit_directions()
    directions = [d.to_json() for d in (*old, *(rho*d for d in old))]
    old7 = np.asarray([[residue(c,7,7) for c in d['a']+d['b']] for d in directions], dtype=np.int32)
    span, annihilator = nullspace(old7,7)
    cosets = [tuple(row) for row in (points7@annihilator.T % 7).tolist()]
    coset_list = sorted(set(cosets)); coset_index = {h:i for i,h in enumerate(coset_list)}
    vertex_cosets = [coset_index[h] for h in cosets]
    dump(args.out/'cosets.json', {'old_exact_directions':directions,
          'seven_times_old_directions_mod7':old7.tolist(), 'old_direction_span_rref':span.tolist(),
          'annihilator_rows':annihilator.tolist(), 'coset_values':coset_list,
          'vertex_coset_indices':vertex_cosets})
    input_maps = [json.loads(args.first_map.read_text())['linear_functional']]
    input_maps.extend(json.loads(args.maps.read_text())['functionals'])
    maps = []; seen = set()
    for values in input_maps:
        functional = tuple(values[:24])
        if functional not in seen:
            maps.append(functional); seen.add(functional)
    maps = maps[:args.max_maps]
    config = {'graph':str(args.graph.resolve()), 'graph_sha256':hashlib.sha256(raw).hexdigest(),
              'program_sha256':verify._sha256(Path(__file__)),
              'first_map_sha256':verify._sha256(args.first_map), 'maps_sha256':verify._sha256(args.maps),
              'old_direction_count':len(directions), 'old_direction_mod7_rank':len(span),
              'annihilator_dimension':len(annihilator), 'coset_count':len(coset_list),
              'coset_size_histogram':dict(Counter(Counter(vertex_cosets).values())),
              'graph_vertices':n, 'graph_edges':len(edges), 'maps_limit':len(maps),
              'timeout_per_map':args.timeout_per_map, 'goal_satisfied':False,
              'split_loop_endpoints':args.split_loop_endpoints,
              'scope':'A verified quotient coloring gives an upper bound; quotient failure is never an original-graph lower bound.'}
    dump(args.out/'config.json', config); records = []; quotient_hashes = {}
    for number,functional in enumerate(maps):
        start = time.monotonic()
        colors5 = list(map(int, points5@np.asarray(functional,dtype=np.int32)%5))
        base_types = list(zip(vertex_cosets,colors5))
        exceptional = {v for u,w in edges if base_types[u]==base_types[w] for v in (u,w)} if args.split_loop_endpoints else set()
        vertex_types = [(*t,v if v in exceptional else -1) for v,t in enumerate(base_types)]
        types = sorted(set(vertex_types)); type_index = {t:i for i,t in enumerate(types)}
        mapping = [type_index[t] for t in vertex_types]
        quotient_edges = {}; loops = []
        for edge_number,(u,v) in enumerate(edges):
            pair = tuple(sorted((mapping[u],mapping[v])))
            if pair[0] == pair[1]: loops.append(edge_number)
            quotient_edges.setdefault(pair,edge_number)
        report = {'map_index':number, 'functional':list(functional), 'quotient_vertices':len(types),
                  'quotient_edges':len(quotient_edges), 'quotient_loop_preimages':len(loops),
                  'exceptional_loop_endpoint_vertices':len(exceptional),
                  'goal_satisfied':False, 'scope':config['scope']}
        dump(args.out/f'map-{number:02}-types.json', {'types':types,'original_vertex_to_type':mapping})
        if loops:
            report.update(status='QUOTIENT_HAS_LOOPS', first_loop_original_edge_indices=loops[:32])
        else:
            quotient = {'schema':'hn-color-quotient-v1', 'vertices':types,
                        'edges':[list(pair) for pair in sorted(quotient_edges)],
                        'metadata':{'scope':'Coloring upper filter; no planar geometry asserted.'}}
            dump(args.out/f'map-{number:02}-quotient.json', quotient)
            formula = cnf_bytes(quotient,6); formula_hash = hashlib.sha256(formula).hexdigest()
            (args.out/f'map-{number:02}.cnf').write_bytes(formula)
            report['cnf_sha256'] = formula_hash
            if formula_hash in quotient_hashes:
                report.update(status='DUPLICATE_QUOTIENT', same_formula_as_map=quotient_hashes[formula_hash])
            else:
                quotient_hashes[formula_hash] = number
                with Solver(name='g4', with_proof=True, bootstrap_with=clauses(quotient,6)) as solver:
                    timer = Timer(args.timeout_per_map, solver.interrupt); timer.daemon = True; timer.start()
                    try:
                        answer = solver.solve_limited(expect_interrupt=True)
                    finally:
                        timer.cancel()
                    report['statistics'] = solver.accum_stats()
                    if answer:
                        positive = {v for v in solver.get_model() if v>0}
                        colors = [next(c for c in range(6) if 6*v+c+1 in positive) for v in range(len(types))]
                        verify.verify_coloring(quotient, {'k':6,'colors':colors}, 6)
                        lifted = {'k':6,'colors':[colors[v] for v in mapping], 'graph_sha256':config['graph_sha256']}
                        verify.verify_coloring(graph,lifted,6)
                        geometry = verify.verify_geometry(graph)
                        if args.graph.read_bytes()!=raw: raise ValueError('Graph changed during upper-bound search')
                        dump(args.out/'coloring.json', lifted)
                        dump(args.out/'verification.json', {'graph_sha256':config['graph_sha256'],
                             'geometry':geometry.__dict__, 'proper_six_coloring_verified':True,
                             'coloring_sha256':verify._sha256(args.out/'coloring.json'),
                             'quotient_colors':colors, 'successful_map':number, 'goal_satisfied':False})
                        report['status'] = 'SIX_COLORING_INDEPENDENTLY_VERIFIED'
                    elif answer is False:
                        proof = ('\n'.join(solver.get_proof())+'\n').encode()
                        (args.out/f'map-{number:02}.drat.gz').write_bytes(gzip.compress(proof,mtime=0))
                        report.update(status='QUOTIENT_UNSAT_UNVERIFIED', proof_sha256=hashlib.sha256(proof).hexdigest())
                    else:
                        report['status'] = 'QUOTIENT_UNKNOWN_TIMEOUT'
        report['elapsed_seconds'] = time.monotonic()-start
        records.append(report); dump(args.out/f'map-{number:02}.json',report)
        print(json.dumps({k:v for k,v in report.items() if k!='functional'}),flush=True)
        status = 'SIX_COLORING_INDEPENDENTLY_VERIFIED' if report['status']=='SIX_COLORING_INDEPENDENTLY_VERIFIED' else 'NO_COLORING_IN_TESTED_QUOTIENTS'
        dump(args.out/'result.json',dict(config,status=status,attempts=records))
        if status=='SIX_COLORING_INDEPENDENTLY_VERIFIED': return


if __name__ == '__main__':
    main()
