"""Find an upper-bound witness by repairing a saved coloring near new points.

Restricted UNSAT or timeout never establishes a graph chromatic lower bound.
The SAT model, if found, is checked against every edge of the full graph.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from threading import Timer
import time

from pysat.solvers import Solver

from .multiquadratic import Point
from .field import Element
from tools import verify as common
from tools import verify_mq


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph', type=Path, required=True)
    parser.add_argument('--base-graph', type=Path, required=True)
    parser.add_argument('--coloring', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-radius', type=int, default=3)
    parser.add_argument('--timeout-per-radius', type=float, default=20)
    parser.add_argument('--allow-improper-base', action='store_true',
                        help='For the identical graph, free every conflicting-edge endpoint initially')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'program.py').write_bytes(Path(__file__).read_bytes())
    raw = args.graph.read_bytes(); graph = json.loads(raw)
    base_raw = args.base_graph.read_bytes(); base = json.loads(base_raw)
    witness = json.loads(args.coloring.read_text()); k = witness['k']
    if witness.get('graph_sha256') != hashlib.sha256(base_raw).hexdigest():
        raise ValueError('Saved coloring does not bind the base graph')
    if (len(witness['colors']) != len(base['vertices'])
            or any(type(c) is not int or not 0<=c<k for c in witness['colors'])):
        raise ValueError('Malformed base color state')
    base_conflicts = [(u,v) for u,v in base['edges'] if witness['colors'][u]==witness['colors'][v]]
    if args.allow_improper_base:
        if raw != base_raw:
            raise ValueError('Improper-base repair requires the identical graph bytes')
    else:
        common.verify_coloring(base, witness, k)
    point_class = Element if graph['schema']=='hn-exact-v1' else Point
    points = [point_class.from_json(p) for p in graph['vertices']]
    point_index = {p:i for i,p in enumerate(points)}
    colors = [-1]*len(points)
    base_to_full = [point_index[point_class.from_json(p)] for p in base['vertices']]
    for v,c in zip(base_to_full, witness['colors']):
        colors[v] = c
    old_edges = {tuple(sorted((base_to_full[u], base_to_full[v])))
                 for u,v in base['edges']}
    if not old_edges.issubset(map(tuple, graph['edges'])):
        raise ValueError('The full graph omits a recorded base edge')
    adjacency = [set() for _ in points]
    for u,v in graph['edges']:
        adjacency[u].add(v); adjacency[v].add(u)
    free = {i for i,c in enumerate(colors) if c < 0}
    if args.allow_improper_base:
        free |= {base_to_full[v] for edge in base_conflicts for v in edge}
    config = {'graph_sha256':hashlib.sha256(raw).hexdigest(),
              'base_graph_sha256':hashlib.sha256(base_raw).hexdigest(),
              'coloring_sha256':common._sha256(args.coloring),
              'k':k, 'new_vertices':sum(c<0 for c in colors),
              'initial_free_vertices':len(free), 'max_radius':args.max_radius,
              'allow_improper_base':args.allow_improper_base,
              'initial_conflicting_edges':len(base_conflicts),
              'timeout_per_radius':args.timeout_per_radius,
              'goal_satisfied':False,
              'scope':'Only SAT proves a full coloring; restricted UNSAT and timeouts prove no lower bound.'}
    dump(args.out/'config.json', config)
    records = []
    for radius in range(args.max_radius+1):
        active = sorted(free); local = {v:i for i,v in enumerate(active)}
        clauses = []
        for v in range(len(active)):
            clauses.append([k*v+c+1 for c in range(k)])
            clauses.extend([-(k*v+c+1), -(k*v+d+1)] for c in range(k) for d in range(c+1,k))
        for u,v in graph['edges']:
            if u in local and v in local:
                clauses.extend([-(k*local[u]+c+1), -(k*local[v]+c+1)] for c in range(k))
            elif u in local:
                clauses.append([-(k*local[u]+colors[v]+1)])
            elif v in local:
                clauses.append([-(k*local[v]+colors[u]+1)])
            elif colors[u] == colors[v]:
                # A new edge between two previously nonadjacent old vertices
                # forces this restricted attempt to fail, without a global claim.
                clauses.append([])
        cnf = (f'p cnf {k*len(active)} {len(clauses)}\n'+
               ''.join(' '.join(map(str,c))+' 0\n' for c in clauses)).encode()
        (args.out/f'radius-{radius}.cnf').write_bytes(cnf)
        dump(args.out/f'radius-{radius}-free-indices.json', active)
        start = time.monotonic()
        with Solver(name='g4', bootstrap_with=clauses) as solver:
            timer = Timer(args.timeout_per_radius, solver.interrupt)
            timer.daemon = True; timer.start()
            try:
                answer = solver.solve_limited(expect_interrupt=True)
            finally:
                timer.cancel()
            record = {'radius':radius, 'free_vertices':len(active), 'clauses':len(clauses),
                      'cnf_sha256':hashlib.sha256(cnf).hexdigest(),
                      'status':'SAT' if answer else 'RESTRICTED_UNSAT' if answer is False else 'UNKNOWN_TIMEOUT',
                      'statistics':solver.accum_stats(), 'elapsed_seconds':time.monotonic()-start}
            records.append(record); print(json.dumps(record), flush=True)
            if answer:
                assignment = {v for v in solver.get_model() if v>0}
                repaired = colors[:]
                for v,i in local.items():
                    repaired[v] = next(c for c in range(k) if k*i+c+1 in assignment)
                result = {'k':k, 'colors':repaired, 'graph_sha256':config['graph_sha256']}
                geometry = (common.verify_geometry(graph) if graph['schema']=='hn-exact-v1'
                            else verify_mq.verify_geometry(graph))
                common.verify_coloring(graph, result, k)
                if args.graph.read_bytes()!=raw:
                    raise ValueError('Graph changed during search')
                dump(args.out/'coloring.json', result)
                dump(args.out/'verification.json', {'geometry':geometry.__dict__,
                      'graph_sha256':config['graph_sha256'], 'proper_coloring_verified':True,
                      'coloring_sha256':common._sha256(args.out/'coloring.json'), 'k':k,
                      'goal_satisfied':False})
                dump(args.out/'result.json', dict(config, status='SAT', attempts=records))
                return
        free |= {v for u in free for v in adjacency[u]}
        dump(args.out/'result.json', dict(config, status='NO_REPAIR_FOUND', attempts=records))


if __name__ == '__main__':
    main()
