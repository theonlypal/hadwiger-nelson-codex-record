"""Proof-producing SAT runs; solver output alone is never a success certificate."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
from threading import Timer

from pysat.solvers import Solver


def clauses(graph, k):
    for v in range(len(graph['vertices'])):
        yield [v*k+c+1 for c in range(k)]
        for c in range(k):
            for d in range(c+1, k):
                yield [-(v*k+c+1), -(v*k+d+1)]
    for u, v in graph['edges']:
        for c in range(k):
            yield [-(u*k+c+1), -(v*k+c+1)]


def cnf_bytes(graph, k):
    n = len(graph['vertices'])
    m = n*(1+k*(k-1)//2) + len(graph['edges'])*k
    return (f'p cnf {n*k} {m}\n' + ''.join(' '.join(map(str,c))+' 0\n' for c in clauses(graph,k))).encode('ascii')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def solve(graph_path, k, out, timeout=60, solver_name='g4'):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    graph_path = Path(graph_path).resolve()
    graph = json.loads(graph_path.read_text())
    cnf_path = out / 'coloring.cnf'
    cnf_path.write_bytes(cnf_bytes(graph,k))
    result = {'schema':'hn-sat-run-v1', 'graph':str(graph_path),
              'graph_sha256':digest(graph_path), 'k':k,
              'vertices':len(graph['vertices']), 'edges':len(graph['edges']),
              'solver':solver_name, 'timeout_seconds':timeout,
              'python':sys.version, 'platform':platform.platform(),
              'cnf_sha256':digest(cnf_path), 'success_certificate':False}
    start = time.monotonic()
    with Solver(name=solver_name, with_proof=True, bootstrap_with=clauses(graph,k)) as sat:
        timer = Timer(timeout, sat.interrupt)
        timer.daemon = True
        timer.start()
        try:
            answer = sat.solve_limited(expect_interrupt=True)
        finally:
            timer.cancel()
        result['elapsed_seconds'] = time.monotonic()-start
        result['statistics'] = sat.accum_stats()
        if answer is True:
            positive = {i for i in sat.get_model() if i>0}
            colors = [next(c for c in range(k) if v*k+c+1 in positive) for v in range(len(graph['vertices']))]
            if any(colors[u] == colors[v] for u,v in graph['edges']):
                raise RuntimeError('SAT model does not give a proper coloring')
            (out/'coloring.json').write_text(json.dumps({'k':k,'colors':colors,'graph_sha256':result['graph_sha256']},indent=2)+'\n')
            result['status'] = 'SAT'
            result['proper_coloring_checked'] = True
        elif answer is False:
            proof = sat.get_proof()
            if proof is None:
                raise RuntimeError('UNSAT solver did not emit proof')
            (out/'proof.drat').write_text('\n'.join(proof)+'\n')
            result['status'] = 'UNSAT_UNVERIFIED'
            result['proof_sha256'] = digest(out/'proof.drat')
        else:
            result['status'] = 'UNKNOWN_TIMEOUT'
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('graph')
    p.add_argument('--k',type=int,required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--timeout',type=float,default=60)
    a = p.parse_args()
    if a.k < 1 or a.timeout <= 0:
        p.error('k and timeout must be positive')
    solve(a.graph,a.k,a.out,a.timeout)
