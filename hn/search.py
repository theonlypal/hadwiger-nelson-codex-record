"""Finite exact candidate families. A found coloring rejects a candidate."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
import numpy as np

from .field import Element, ZERO, ONE, ZETA, SQRT5, FIELD_METADATA
from .solve import solve


def all_candidate_edges(points, window=1e-8):
    """Numerical broad phase, exact acceptance. No completeness claim.

    The returned edge list defines the graph. Even if the broad phase missed
    a pair, this remains an exact unit-distance subgraph. Non-6-colorability
    of that subgraph would suffice. Norm equality alone admits each edge.
    """
    coords = np.array([[z.approximate().real,z.approximate().imag] for z in points])
    edges = []
    directions = {}
    for u in range(len(points)):
        delta = coords[u+1:] - coords[u]
        dsq = np.einsum('ij,ij->i',delta,delta)
        for offset in np.flatnonzero(np.abs(dsq-1)<window):
            v = u+1+int(offset)
            d = points[v]-points[u]
            if d not in directions:
                directions[d] = d.is_unit()
            if directions[d]:
                edges.append([u,v])
    return edges


def known_direction_edges(points,directions):
    directions = tuple(sorted(set(directions),key=lambda z:(z.denominator,z.coefficients)))
    if not all(d.is_unit() for d in directions):
        raise ValueError('A search direction is not exactly unit length')
    index = {z:i for i,z in enumerate(points)}
    edges = set()
    for u,z in enumerate(points):
        for d in directions:
            v = index.get(z+d)
            if v is not None and u!=v:
                edges.add((min(u,v),max(u,v)))
    return [list(e) for e in sorted(edges)]


def export(points,out,metadata,directions=None):
    start = time.monotonic()
    points = sorted(set(points),key=lambda z:(z.denominator,z.coefficients))
    print(json.dumps({'stage':'edges','vertices':len(points),'candidate':metadata}),flush=True)
    edges = all_candidate_edges(points) if directions is None else known_direction_edges(points,directions)
    graph = {'schema':'hn-exact-v1','field':FIELD_METADATA,
             'vertices':[z.to_json() for z in points],'edges':edges,
             'metadata':dict(metadata,edge_rule='Exact norm equality; candidate enumeration may omit unit pairs',
                             construction_seconds=time.monotonic()-start)}
    Path(out).write_text(json.dumps(graph,separators=(',',':'))+'\n')
    print(json.dumps({'stage':'graph_written','path':str(out),'vertices':len(points), 'edges':len(edges)}),flush=True)
    return graph


def spindle_rotation():
    return (7+(2*ZETA**7-ONE)*SQRT5)/8


def candidate_points(seed,kind,parameter):
    if kind == 'translated-star':
        offsets = [ZERO]+[ZETA**(j*(42//parameter)) for j in range(parameter)]
        return {z+w for z in seed for w in offsets}
    if kind == 'rotation-union':
        rho = spindle_rotation()
        return {z*rho**j for j in range(parameter) for z in seed}
    if kind == 'rational-rotation-union':
        r = parameter
        rho = ((1-3*r*r)+(2*ZETA**7-ONE)*(2*r))/(1+3*r*r)
        assert rho.is_unit()
        return set(seed)|{z*rho for z in seed}
    if kind == 'spindle-sum':
        rho = spindle_rotation()
        # A small exact unit hexagon in the rotated factor.
        offsets = [ZERO]+[rho*ZETA**(7*j) for j in range(6)]
        return {z+w for z in seed for w in offsets}
    raise ValueError(kind)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',required=True)
    p.add_argument('--family',choices=['translated-star','rotation-union','rational-rotation-union','spindle-sum'],required=True)
    p.add_argument('--parameter',type=int,default=2)
    p.add_argument('--out',required=True)
    p.add_argument('--timeout',type=float,default=60)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True,exist_ok=False)
    source = json.loads(Path(a.seed).read_text())
    if source['field']['cyclotomic_order']!=42 or source['field']['radicand']!=5:
        raise ValueError('Seed field does not match')
    seed = [Element.from_json(z) for z in source['vertices']]
    points = candidate_points(seed,a.family,a.parameter)
    export(points,out/'graph.json',{'family':a.family,'parameter':a.parameter,'seed':str(Path(a.seed).resolve())})
    result = solve(out/'graph.json',6,out/'six-color',a.timeout)
    if result['status']=='SAT':
        solve(out/'graph.json',5,out/'five-color',a.timeout)


if __name__=='__main__':
    main()
