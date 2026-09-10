#!/usr/bin/env python3
"""Bounded exact unit-direction search and additive-F5 obstruction filter.

UNSAT here only excludes additive coefficient colorings modulo5. It is not
a certificate of graph non-five-colorability. No large point set is built.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import threading
import time

from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from hn.cyclotomic_tower import Element, ONE, ZETA, SQRT5, SQRT33
from hn.construct import unit_directions
from hn.mixed_cyclotomic import B, HALF_THETA3, THETA4


def residue(direction):
    if direction.denominator % 5 == 0:
        return None
    inverse=pow(direction.denominator,-1,5)
    values=tuple(v*inverse%5 for v in direction.coefficients)
    first=next((v for v in values if v),None)
    if first is None:
        raise ArithmeticError('An exact unit direction reduced to zero')
    scale=pow(first,-1,5)
    return tuple(v*scale%5 for v in values)


def encode(directions):
    pool=IDPool();cnf=CNF()
    def onehot(key):
        ids=[pool.id((key,c)) for c in range(5)]
        cnf.append(ids)
        for i in range(5):
            for j in range(i+1,5):cnf.append([-ids[i],-ids[j]])
        return ids
    variables=[onehot(('coefficient',j)) for j in range(48)]
    assert (1,)+(0,)*47 in directions
    cnf.append([variables[0][1]])
    for d,direction in enumerate(sorted(directions)):
        previous=None
        for term,(j,a) in enumerate((j,a) for j,a in enumerate(direction) if a):
            state=onehot(('sum',d,term))
            if previous is None:
                for v in range(5):cnf.append([-variables[j][v],state[a*v%5]])
            else:
                for r in range(5):
                    for v in range(5):cnf.append([-previous[r],-variables[j][v],state[(r+a*v)%5]])
            previous=state
        cnf.append([-previous[0]])
    return cnf,variables


def screen(directions,out,timeout):
    out.mkdir()
    cnf,variables=encode(directions)
    cnf.to_file(str(out/'linear.cnf'))
    started=time.monotonic()
    with Solver(name='g4',bootstrap_with=cnf.clauses,with_proof=True) as solver:
        timer=threading.Timer(timeout,solver.interrupt);timer.start()
        try:
            answer=solver.solve_limited(expect_interrupt=True)
            result={'status':'SAT_VERIFIED' if answer is True else 'UNSAT_UNVERIFIED' if answer is False else 'UNKNOWN_TIMEOUT',
                    'scope':'additive coefficient F5 colorings only; no graph lower bound',
                    'directions':len(directions),'cnf_variables':cnf.nv,'cnf_clauses':len(cnf.clauses),
                    'elapsed_seconds':time.monotonic()-started,'stats':solver.accum_stats()}
            if answer is True:
                model=set(solver.get_model())
                coeff=[next(c for c in range(5) if variables[j][c] in model) for j in range(48)]
                assert all(sum(a*b for a,b in zip(coeff,d))%5 for d in directions)
                result['linear_functional']=coeff
            elif answer is False:
                proof=solver.get_proof()
                (out/'proof.drat').write_text('\n'.join(proof)+'\n')
            (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
            print(json.dumps({'stage':out.name,**result}),flush=True)
            return result
        finally:timer.cancel()


def run(seed_result,out,timeout=60):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).read_bytes();(out/'program.py').write_bytes(source)
    seed_path=Path(seed_result);seed=json.loads(seed_path.read_text())
    constraints={tuple(d) for d in seed['directions']}
    old_functional=seed['linear_functional']
    origins={};rejected=[]
    haugland=[Element.from_legacy(d) for d in unit_directions()]
    def include(d,formula):
        if d in origins:return
        if not d.is_unit():raise ArithmeticError('Unit identity failed: '+formula)
        r=residue(d)
        if r is None:
            rejected.append({'formula':formula,'reason':'denominator divisible by5','point':d.to_json()})
            return
        origins[d]={'formula':formula,'residue':r,
                    'zero_under_seed_functional':sum(a*b for a,b in zip(old_functional,r))%5==0}
        constraints.add(r)
    results=[]
    stages=[('half-angle-products',range(-1,2),range(0,2)),
            ('multiple-rotation-products',range(-2,3),range(-2,3))]
    for name,ks,ms in stages:
        for k in ks:
            for m in ms:
                rho=HALF_THETA3**k*THETA4**m
                for j,u in enumerate(haugland):include(u*rho,f'u[{j}]*half_theta3^{k}*theta4^{m}')
        result=screen(constraints,out/name,timeout);results.append(result)
        if result['status']=='UNSAT_UNVERIFIED':break
    if results[-1]['status']!='UNSAT_UNVERIFIED':
        for radical_name,radical in [('1',ONE),('sqrt5',SQRT5),('sqrt33',SQRT33),('sqrt165',SQRT5*SQRT33)]:
            for q in range(1,6):
                r=q*radical
                rho=(ONE-3*r*r+2*B*r)/(ONE+3*r*r)
                for j in range(42):include(rho*ZETA**j,f'circle({q}*{radical_name})*zeta^{j}')
        result=screen(constraints,out/'rational-circle-directions',timeout);results.append(result)
    record={'schema':'hn-exact-direction-pool-v1','field':{'cyclotomic_order':42,'radicands':[5,33]},
            'seed':str(seed_path.resolve()),'seed_sha256':hashlib.sha256(seed_path.read_bytes()).hexdigest(),
            'program_sha256':hashlib.sha256(source).hexdigest(),
            'new_directions':[dict(v,point=d.to_json()) for d,v in origins.items()],
            'excluded_directions':rejected,'constraints':sorted(constraints),'results':results,
            'success_certificate':False}
    (out/'pool.json').write_text(json.dumps(record,separators=(',',':'))+'\n')
    print(json.dumps({'stage':'complete','exact_directions':len(origins),
                      'zero_under_seed_functional':sum(v['zero_under_seed_functional'] for v in origins.values()),
                      'excluded_5_denominators':len(rejected),'final_status':results[-1]['status']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed-result',type=Path,default=Path('runs/mixed-haugland-g3-v31-linear5/result.json'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--timeout',type=float,default=60)
    args=parser.parse_args();run(args.seed_result,args.out,args.timeout)
