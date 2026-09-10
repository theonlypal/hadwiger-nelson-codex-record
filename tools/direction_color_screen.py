#!/usr/bin/env python3
"""Search additive F5 colorings of an existing tower graph's directions.

This never constructs a new graph. A successful linear functional colors
every translate of every observed direction. Failure concerns this ansatz
only. Every supplied graph edge is checked against the final coloring.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import threading
import time

import numpy as np
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver


def run(graph_path, out, timeout=60):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).read_bytes()
    (out/'program.py').write_bytes(source)
    started = time.monotonic()
    graph_path = Path(graph_path)
    raw = graph_path.read_bytes()
    graph = json.loads(raw)
    assert graph['schema'] == 'hn-cyclotomic-tower-v1'
    assert graph['field']['cyclotomic_order'] == 42 and graph['field']['radicands'] == [5,33]
    @lru_cache(None)
    def residue(value):
        value = Fraction(value)
        if value.denominator % 5 == 0:
            raise ValueError('A coordinate denominator is divisible by5')
        return value.numerator*pow(value.denominator, -1, 5) % 5
    coords = np.array([[residue(v) for part in z['parts'] for v in part]
                       for z in graph['vertices']], dtype=np.int16)
    directions = set()
    for u, v in graph['edges']:
        delta = (coords[v]-coords[u]) % 5
        indices = np.flatnonzero(delta)
        if not len(indices):
            raise ValueError('An edge collapses modulo5, so no linear coloring can work')
        delta = delta*pow(int(delta[indices[0]]), -1, 5) % 5
        directions.add(tuple(map(int,delta)))
    directions = sorted(directions)
    width = coords.shape[1]
    assert tuple([1]+[0]*(width-1)) in directions
    pool = IDPool()
    cnf = CNF()
    def onehot(key):
        variables = [pool.id((key, c)) for c in range(5)]
        cnf.append(variables)
        for i in range(5):
            for j in range(i+1,5):
                cnf.append([-variables[i],-variables[j]])
        return variables
    variables = [onehot(('coefficient',j)) for j in range(width)]
    cnf.append([variables[0][1]])
    for d, direction in enumerate(directions):
        nonzero = [(j,a) for j,a in enumerate(direction) if a]
        previous = None
        for term,(j,a) in enumerate(nonzero):
            state = onehot(('sum',d,term))
            if previous is None:
                for value in range(5):
                    cnf.append([-variables[j][value],state[a*value%5]])
            else:
                for accum in range(5):
                    for value in range(5):
                        cnf.append([-previous[accum],-variables[j][value],state[(accum+a*value)%5]])
            previous = state
        cnf.append([-previous[0]])
    cnf.to_file(str(out/'linear.cnf'))
    print(json.dumps({'stage':'solve','projective_mod5_directions':len(directions),
                      'variables':pool.top,'clauses':len(cnf.clauses)}),flush=True)
    with Solver(name='g4',bootstrap_with=cnf.clauses) as solver:
        timer = threading.Timer(timeout,solver.interrupt)
        timer.start()
        try:
            answer=solver.solve_limited(expect_interrupt=True)
            status='SAT_VERIFIED' if answer is True else 'UNSAT_UNVERIFIED' if answer is False else 'UNKNOWN_TIMEOUT'
            result={'schema':'hn-direction-color-screen-v1','graph':str(graph_path.resolve()),
                    'graph_sha256':hashlib.sha256(raw).hexdigest(),'program_sha256':hashlib.sha256(source).hexdigest(),
                    'prime':5,'directions':directions,'status':status,'timeout_seconds':timeout,
                    'success_certificate':False,'solver_stats':solver.accum_stats()}
            if answer is True:
                model=set(solver.get_model())
                coefficients=[next(c for c in range(5) if variables[j][c] in model) for j in range(width)]
                assert all(sum(a*b for a,b in zip(coefficients,d))%5 for d in directions)
                colors=(coords@np.array(coefficients,dtype=np.int16))%5
                assert all(colors[u]!=colors[v] for u,v in graph['edges'])
                result.update({'linear_functional':coefficients,'vertices_checked':len(coords),
                               'edges_checked':len(graph['edges'])})
                coloring={'graph_sha256':result['graph_sha256'],'k':5,'colors':list(map(int,colors))}
                (out/'coloring.json').write_text(json.dumps(coloring,separators=(',',':'))+'\n')
            result['elapsed_seconds']=time.monotonic()-started
            (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps({k:v for k,v in result.items() if k!='directions'}),flush=True)
            return result
        finally:
            timer.cancel()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--timeout',type=float,default=60)
    args=parser.parse_args()
    run(args.graph,args.out,args.timeout)
