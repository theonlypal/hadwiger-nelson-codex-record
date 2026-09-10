#!/usr/bin/env python3
"""Exact finite residue diagnostics, never a Euclidean lower-bound certificate.

The default report verifies the two characteristic-five projections and the
reduced characteristic-three projection of Haugland's direction families.
Optionally export/color the 729-vertex latter quotient, with finite-graph scope.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import sys

from flint import nmod_poly
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hn.construct import unit_directions
from hn.field import ONE, ZETA
from hn.solve import solve

PHI = [1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1]


def key(poly, degree=6):
    return tuple(int(poly[j]) for j in range(degree))


def reduce_element(z, prime, factor):
    inv = pow(z.denominator, -1, prime)
    assert not any(z.coefficients[12:])
    return nmod_poly([v*inv % prime for v in z.coefficients[:12]], prime) % factor


def residue_report():
    u = unit_directions()[1]
    directions = tuple(dict.fromkeys(unit_directions() + tuple(u*u*ZETA**j for j in range(42))))
    assert len(directions) == 126 and all(z.is_unit() for z in directions)
    q = [1]+[0]*10+[-4]
    def lattice_residue(z):
        assert 7 % z.denominator == 0
        return sum(a*b*(7//z.denominator) for a,b in zip(q,z.coefficients)) % 7
    generators = [ZETA**j for j in range(12)]+[u*ZETA**j for j in range(12)]
    assert all(lattice_residue(z) == 0 for z in generators)
    assert lattice_residue(u*u) == 3
    report = {'scope':'finite residue and lattice diagnostics; no Euclidean chromatic lower bound',
              'directions':126,'lattice_separation':{'covector':q,'modulus':7,
              'applied_to':'coefficient vector of 7z', 'generators_checked':24,
              'old_generator_residues':[lattice_residue(z) for z in generators],
              'u1_squared_residue':lattice_residue(u*u), 'u1_squared':(u*u).to_json()},
              'factors':[], 'success_certificate':False}
    for prime in (3,5):
        phi=nmod_poly(PHI,prime)
        for factor,multiplicity in phi.factor()[1]:
            t=reduce_element(ZETA,prime,factor)
            odd=reduce_element(u,prime,factor)
            one=nmod_poly([1],prime)
            t_order=next(k for k in range(1,127) if t.pow_mod(k,factor)==one)
            odd_order=next(k for k in range(1,127) if odd.pow_mod(k,factor)==one)
            conjugate_t=reduce_element(ZETA.conjugate(),prime,factor)
            assert t.pow_mod(prime**3,factor)==conjugate_t
            old={key(reduce_element(z,prime,factor)) for z in unit_directions()}
            new={key(reduce_element(z,prime,factor)) for z in directions}
            group_order=prime**3+1
            assert all(nmod_poly(list(v),prime).pow_mod(group_order,factor)==one for v in new)
            assert len(new)==group_order  # They are all roots of X^(p^3+1)-1.
            relations=[]
            for exponent in range(1,7):
                relations.append({'u1_power':exponent,'zeta_powers':[j for j in range(t_order)
                                  if odd.pow_mod(exponent,factor)==t.pow_mod(j,factor)]})
            report['factors'].append({'prime':prime,'factor_coefficients':list(key(factor,7)),
                  'factor_multiplicity':multiplicity,'irreducible_degree':factor.degree(),
                  'zeta_order':t_order,'u1_order':odd_order,'u1_residue':list(key(odd)),
                  'U84_distinct_residues':len(old),'U126_distinct_residues':len(new),
                  'unit_norm_group_order':group_order,'conjugation_is_frobenius_cubed':True,
                  'power_relations':relations,'directions':sorted(new)})
    return report


def finite_graph(report):
    record=next(x for x in report['factors'] if x['prime']==3)
    directions=record['directions']
    vertices=list(itertools.product(range(3),repeat=6))
    lookup={v:j for j,v in enumerate(vertices)}
    edges=set()
    for i,v in enumerate(vertices):
        for d in directions:
            w=tuple((a+b)%3 for a,b in zip(v,d));j=lookup[w]
            assert i!=j
            edges.add((min(i,j),max(i,j)))
    assert len(vertices)==729 and len(edges)==10206
    return {'schema':'hn-finite-cayley-diagnostic-v1',
            'warning':'Finite-field graph. These are not coordinates in the Euclidean plane.',
            'field':{'prime':3,'factor_coefficients':record['factor_coefficients']},
            'vertices':vertices,'edges':sorted(edges),'directions':directions,
            'success_certificate':False}


def quadratic_sign(a,b):
    """Exact sign of a+b*sqrt(5), for rational integers a,b."""
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if (a>0)==(b>0):return (a>0)-(a<0)
    comparison=(a*a>5*b*b)-(a*a<5*b*b)
    return comparison if a>0 else -comparison


def spectrum_report(record):
    """Every additive-character eigenvalue, from integer dot-product counts.

    Characters in Fp^6 are indexed by all coefficient vectors, regardless of
    which polynomial basis was used for the field. No eigenvalue tolerance
    enters either the ordering or the Hoffman lower-bound calculation.
    """
    p=record['prime']
    vectors=np.array(list(itertools.product(range(p),repeat=6)),dtype=np.int16)
    directions=np.array(record['directions'],dtype=np.int16)
    residues=(vectors@directions.T)%p
    counts=np.stack([(residues==j).sum(axis=1) for j in range(p)],axis=1)
    histogram=Counter(map(tuple,counts.tolist()))
    assert sum(histogram.values())==p**6
    rows=[]
    for c,multiplicity in sorted(histogram.items()):
        assert all(c[j]==c[-j] for j in range(1,p))
        # Store eigenvalues as (a+b*sqrt(5))/2, even in characteristic3.
        a,b=(2*(c[0]-c[1]),0) if p==3 else (2*c[0]-c[1]-c[2],c[1]-c[2])
        rows.append({'counts':c,'multiplicity':multiplicity,'twice_eigenvalue':[a,b]})
    minimum=min((r['twice_eigenvalue'] for r in rows),key=lambda v:v[0])
    for row in rows:
        a,b=row['twice_eigenvalue']
        if quadratic_sign(a-minimum[0],b-minimum[1])<0:minimum=[a,b]
    assert all(quadratic_sign(a-minimum[0],b-minimum[1])>=0
               for a,b in (r['twice_eigenvalue'] for r in rows))
    # Hoffman: chi >= 1+degree/(-lambda_min). Test whether this exceeds k.
    degree=len(directions);lower=1
    while quadratic_sign(2*degree+(lower-1)*minimum[0],(lower-1)*minimum[1])>0:
        lower+=1
    expected=[-16,0] if p==3 else [-18,-10]
    assert minimum==expected and lower==(5 if p==3 else 8)
    return {'scope':'finite residue Cayley graph only; not a Euclidean lower bound',
            'prime':p,'factor_coefficients':record['factor_coefficients'],
            'vertices':p**6,'degree':degree,'character_count_tuples':rows,
            'twice_minimum_eigenvalue':minimum,'hoffman_chromatic_lower_bound':lower,
            'success_certificate':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--color-reduced-three',action='store_true')
    parser.add_argument('--timeout',type=float,default=30)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).read_bytes();(args.out/'program.py').write_bytes(source)
    report=residue_report()
    report['program_sha256']=hashlib.sha256(source).hexdigest()
    (args.out/'residues.json').write_text(json.dumps(report,indent=2)+'\n')
    for index,record in enumerate(report['factors']):
        spectrum=spectrum_report(record)
        (args.out/f'spectrum-p{record["prime"]}-factor{index}.json').write_text(json.dumps(spectrum,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='factors'}),flush=True)
    if args.color_reduced_three:
        path=args.out/'finite-graph.json'
        path.write_text(json.dumps(finite_graph(report),separators=(',',':'))+'\n')
        for k in (4,5,6):
            result=solve(path,k,args.out/f'k{k}',timeout=args.timeout)
            if result['status']=='SAT':break


if __name__=='__main__':main()
