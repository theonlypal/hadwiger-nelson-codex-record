"""Glue exact sqrt(2)-rotated Heule copies using actual unit bridge pairs."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np

from .mq_construct import heule_union_points, SQRT2_RATIONAL_ROTATION, export, canonical_points
from .multiquadratic import MQ, Point
from .solve import solve


def anchor_pairs(points):
    target=MQ(Fraction(3,4))
    coords=np.array([[p.x.approximate(),p.y.approximate()] for p in points])
    pairs=[]
    for i,p in enumerate(points):
        delta=coords[i+1:]-coords[i]
        dsq=np.einsum('ij,ij->i',delta,delta)
        for off in np.flatnonzero(np.abs(dsq-.75)<1e-8):
            j=i+1+int(off)
            if (points[j]-p).norm_squared()==target:
                pairs.append((i,j))
    return pairs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',required=True)
    p.add_argument('--anchors',type=int,default=3)
    p.add_argument('--timeout',type=float,default=180)
    a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    (out/'program.py').write_bytes(Path(__file__).read_bytes())
    points=heule_union_points()
    pairs=anchor_pairs(points)
    freq=Counter(i for pair in pairs for i in pair)
    # Prefer anchors with multiple exact bridge endpoints. Break ties by stable index.
    selected=sorted(freq,key=lambda i:(-freq[i],i))[:a.anchors]
    config={'family':'anchored-sqrt2-rotations','anchors':selected,'anchor_pairs':pairs,
            'exact_squared_anchor_distance':'3/4','rotation':'1/3+i*2*sqrt(2)/3',
            'rotation_norm_squared':'1','bridge_identity':'(3/4)*|rho-1|^2 = (3/4)*(4/3) = 1',
            'timeout_seconds':a.timeout,'program_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    rho=SQRT2_RATIONAL_ROTATION
    union=set(points)
    bridges=[]
    for index in selected:
        center=points[index]
        union.update(center+rho*(z-center) for z in points)
        for i,j in pairs:
            if index in (i,j):
                q=points[j if i==index else i]
                image=center+rho*(q-center)
                if not (image-q).is_unit():
                    raise ArithmeticError('A promised bridge is not unit')
                bridges.append((q,image))
    graph=export(union,out/'graph.json',config)
    index={p:i for i,p in enumerate(canonical_points(union))}
    edges={tuple(e) for e in graph['edges']}
    bridge_indices=[sorted((index[u],index[v])) for u,v in bridges]
    if not all(tuple(e) in edges for e in bridge_indices):
        raise ArithmeticError('Broadphase omitted a required exact bridge')
    (out/'bridge-provenance.json').write_text(json.dumps({'guaranteed_unit_bridges':bridge_indices},indent=2)+'\n')
    result=solve(out/'graph.json',6,out/'six-color',a.timeout)
    if result['status']=='SAT':
        solve(out/'graph.json',5,out/'five-color',60)


if __name__=='__main__':
    main()
