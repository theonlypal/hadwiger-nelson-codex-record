"""Exactly couple a full tower seed with its Haugland-u1 rotation about zero.

Every seed edge and its rotated image are retained. Additional cross-copy
edges are proposed numerically, then admitted only by exact norm equality.
No induced-graph or non-colorability conclusion follows from construction.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .construct import unit_directions
from .cyclotomic_tower import Element, FIELD_METADATA
from .mixed_cyclotomic import canonical_points
from .solve_external import solve_external


def dump(path,value):
    Path(path).write_text(json.dumps(value,separators=(',',':'))+'\n')


def construct(seed_path,out,window=1e-8):
    seed_path=Path(seed_path).resolve();out=Path(out)
    out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).read_bytes();(out/'program.py').write_bytes(source)
    raw=seed_path.read_bytes();seed=json.loads(raw)
    if seed.get('schema')!='hn-cyclotomic-tower-v1':raise ValueError('Expected exact tower seed')
    if seed['field']['cyclotomic_order']!=42 or seed['field']['radicands']!=[5,33]:
        raise ValueError('Seed field does not match the supported tower')
    started=time.monotonic()
    original=tuple(Element.from_json(v) for v in seed['vertices'])
    if len(set(original))!=len(original):raise ValueError('Seed has duplicate coordinates')
    rho=Element.from_legacy(unit_directions()[1])
    if not rho.is_unit():raise ArithmeticError('Rotation is not exactly unit')
    rotated=tuple(rho*z for z in original)
    points=canonical_points(original+rotated);index={z:i for i,z in enumerate(points)}
    original_ids=[index[z] for z in original];rotated_ids=[index[z] for z in rotated]
    masks=[0]*len(points)
    for i in original_ids:masks[i]|=1
    for i in rotated_ids:masks[i]|=2
    counts=Counter(masks)
    print(json.dumps({'stage':'points','seed_vertices':len(original),
                      'union_vertices':len(points),'intersection_vertices':counts[3],
                      'exclusive_vertices_per_copy':[counts[1],counts[2]],
                      'seconds':time.monotonic()-started}),flush=True)
    norm_cache={}
    edges=set();inherited=set()
    for a,b in seed['edges']:
        d=original[b]-original[a]
        if d not in norm_cache:norm_cache[d]=d.is_unit()
        if not norm_cache[d]:raise ArithmeticError('Seed edge is not exactly unit')
        for ids in (original_ids,rotated_ids):
            edge=tuple(sorted((ids[a],ids[b])))
            edges.add(edge);inherited.add(edge)
    inherited_count=len(inherited)
    coordinates=np.array([[z.approximate().real,z.approximate().imag] for z in rotated])
    candidates=0;rejected=0
    for a,z in enumerate(original):
        approximation=z.approximate()
        delta=coordinates-np.array([approximation.real,approximation.imag])
        distance=np.einsum('ij,ij->i',delta,delta)
        for b in np.flatnonzero(np.abs(distance-1)<window):
            edge=tuple(sorted((original_ids[a],rotated_ids[int(b)])))
            if edge in edges:continue
            candidates+=1
            d=rotated[int(b)]-z
            if d not in norm_cache:norm_cache[d]=d.is_unit()
            if norm_cache[d]:edges.add(edge)
            else:rejected+=1
        if a and a%10000==0:
            print(json.dumps({'stage':'cross_copy_edges','original_vertices_checked':a,
                              'edges_so_far':len(edges),'seconds':time.monotonic()-started}),flush=True)
    sorted_edges=sorted(edges)
    categories=Counter();strict=[];strict_vertices=set()
    for a,b in sorted_edges:
        common=masks[a]&masks[b]
        if not common:
            categories['strict_cross_copy']+=1;strict.append([a,b]);strict_vertices.update((a,b))
        elif common==3:categories['within_both_copies']+=1
        elif common==1:categories['within_original_only']+=1
        elif common==2:categories['within_rotated_only']+=1
    inherited_preserved=all(tuple(sorted((ids[a],ids[b]))) in edges
                            for ids in (original_ids,rotated_ids) for a,b in seed['edges'])
    assert inherited_preserved
    report={'scope':'exact Euclidean candidate; no coloring conclusion',
            'seed':str(seed_path),'seed_sha256':hashlib.sha256(raw).hexdigest(),
            'seed_vertices':len(original),'seed_edges':len(seed['edges']),
            'vertices':len(points),'edges':len(sorted_edges),'intersection_vertices':counts[3],
            'exclusive_vertices_per_copy':[counts[1],counts[2]],
            'edge_categories':dict(categories),'strict_cross_copy_incident_vertices':len(strict_vertices),
            'inherited_union_edges':inherited_count,'inherited_seed_edges_all_preserved':inherited_preserved,
            'new_edges_beyond_inherited':len(edges)-inherited_count,
            'numerical_proposals_checked_exactly':candidates,'rejected_nonunit_proposals':rejected,
            'unique_directions_norm_checked':len(norm_cache),'rotation':rho.to_json(),
            'construction_seconds':time.monotonic()-started,'success_certificate':False}
    metadata={'family':'mixed-haugland-g3-v31-union-odd-rotation','rotation_formula':'u1 from Haugland direction orbit; rotate about0',
              'edge_rule':'Retain both exact seed copies; cross-copy numerical candidates require exact norm equality',
              'induced_graph_claimed':False,'seed_metadata':seed.get('metadata',{}),
              'published_sources':['https://arxiv.org/html/2608.04542v1','https://arxiv.org/html/1805.12181v1'],
              'construction_program_sha256':hashlib.sha256(source).hexdigest(),
              'arithmetic_program_sha256':hashlib.sha256(Path(__file__).with_name('cyclotomic_tower.py').read_bytes()).hexdigest(),
              **report}
    graph={'schema':'hn-cyclotomic-tower-v1','field':dict(FIELD_METADATA),
           'vertices':[z.to_json() for z in points],'edges':sorted_edges,'metadata':metadata}
    dump(out/'graph.json',graph);dump(out/'coupling.json',report)
    dump(out/'copy-memberships.json',{'original_vertex_ids':original_ids,'rotated_vertex_ids':rotated_ids,
                                     'vertex_masks':masks,'strict_cross_copy_edges':strict})
    print(json.dumps({'stage':'graph_written',**report}),flush=True)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=Path,default=Path('artifacts/mixed_haugland_g3_v31_tower.json'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--timeout',type=float,default=300)
    parser.add_argument('--no-solve',action='store_true')
    args=parser.parse_args()
    if not 0<args.timeout<=300:parser.error('Solver timeout must lie in (0,300]')
    report=construct(args.seed,args.out)
    if not args.no_solve and report['edge_categories'].get('strict_cross_copy',0):
        solve_external(args.out/'graph.json',args.out/'six-color',timeout=args.timeout,configuration='sat',k=6)


if __name__=='__main__':main()
