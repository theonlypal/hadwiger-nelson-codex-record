"""Coloring-guided exact vertex augmentation; six-color SAT is the target test.

New centers z+u are scored by the colors of their current unit neighbors.
Approximate grouping only proposes centers. Exact equality decides all
neighbors and all accepted edges. Every graph has its own unrestricted CNF.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .field import Element
from .search import export
from .solve import solve


def candidate_pool(points,colors,directions,max_proposals=2000):
    coords=np.array([[z.approximate().real,z.approximate().imag] for z in points])
    offsets=np.array([[z.approximate().real,z.approximate().imag] for z in directions])
    # The maximum live table has n*m integer coordinates and one origin index.
    table=np.rint((coords[:,None,:]+offsets[None,:,:]).reshape(-1,2)*1e7).astype(np.int64)
    keys=np.ascontiguousarray(table).view(np.dtype((np.void,table.dtype.itemsize*2))).ravel()
    unique,first,inverse,counts=np.unique(keys,return_index=True,return_inverse=True,return_counts=True)
    masks=np.zeros(len(unique),dtype=np.uint8)
    bits=np.repeat(np.left_shift(np.uint8(1),np.array(colors,dtype=np.uint8)),len(directions))
    np.bitwise_or.at(masks,inverse,bits)
    popcount=np.array([int(i).bit_count() for i in range(256)],dtype=np.uint8)
    scores=popcount[masks].astype(np.int64)*100000+counts
    qualified=np.flatnonzero((counts>=6)&(popcount[masks]>=4))
    ranked=qualified[np.argsort(scores[qualified])[::-1]][:max_proposals]
    index={z:i for i,z in enumerate(points)}
    proposals=[]
    seen=set()
    for group in ranked:
        origin,direction=divmod(int(first[group]),len(directions))
        point=points[origin]+directions[direction]
        if point in index or point in seen:
            continue
        seen.add(point)
        neighbors={index[point-u] for u in directions if point-u in index}
        palette={colors[v] for v in neighbors}
        if len(neighbors)>=6:
            proposals.append((len(palette),len(neighbors),point,sorted(neighbors)))
    proposals.sort(key=lambda item:(item[0],item[1]),reverse=True)
    return proposals,{'approximate_groups':len(unique),'eligible_groups':len(qualified),
                      'exact_proposals':len(proposals),'palette_histogram':dict(Counter(p[0] for p in proposals))}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',required=True)
    p.add_argument('--coloring',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--rounds',type=int,default=10)
    p.add_argument('--batch',type=int,default=100)
    p.add_argument('--timeout',type=float,default=120)
    p.add_argument('--max-directions',type=int,default=400)
    p.add_argument('--max-proposals',type=int,default=2000)
    p.add_argument('--native',action='store_true',help='Use the proof-producing independently verified native solver')
    p.add_argument('--linear-five-rule',type=Path,help='JSON result containing a known tower-field linear_functional')
    p.add_argument('--screen-linear-five',action='store_true',help='Try a bounded additive five-coloring search on each tower export')
    a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=False)
    (out/'program.py').write_bytes(Path(__file__).read_bytes())
    config=dict(vars(a),algorithm='rank exact centers by neighbor palette count then degree',
                program_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    config['linear_five_rule']=str(a.linear_five_rule) if a.linear_five_rule is not None else None
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    graph=json.loads(Path(a.seed).read_text())
    coloring=json.loads(Path(a.coloring).read_text())
    if graph['schema']=='hn-exact-v1' and graph['field']['cyclotomic_order']==42 and graph['field']['radicand']==5:
        point_type=Element
        export_graph=export
    elif graph['schema']=='hn-multiquadratic-v1' and graph['radicands']==[2,3,5,11]:
        from .multiquadratic import Point
        from .mq_construct import export as export_mq
        point_type=Point
        export_graph=export_mq
    elif graph['schema']=='hn-cyclotomic-tower-v1' and graph['field'].get('cyclotomic_order')==42 and graph['field'].get('radicands')==[5,33]:
        from .cyclotomic_tower import Element as TowerElement
        from .mixed_cyclotomic import export as export_tower
        point_type=TowerElement
        export_graph=export_tower
    else:
        raise ValueError('Unsupported exact field')
    points=[point_type.from_json(v) for v in graph['vertices']]
    colors=coloring['colors']
    if coloring['k']!=6 or len(colors)!=len(points) or any(colors[u]==colors[v] for u,v in graph['edges']):
        raise ValueError('Seed coloring is not a proper6coloring')
    frequencies=Counter()
    for u,v in graph['edges']:
        d=points[v]-points[u]
        frequencies[d]+=1;frequencies[-d]+=1
    directions=[d for d,n in frequencies.most_common(a.max_directions)]
    if not all(d.is_unit() for d in directions):
        raise ValueError('Seed has nonunit edge')
    for iteration in range(a.rounds):
        start=time.monotonic()
        candidates,stats=candidate_pool(points,colors,directions,a.max_proposals)
        destination=out/f'round-{iteration:03d}'
        destination.mkdir()
        if not candidates:
            (destination/'screen.json').write_text(json.dumps(dict(stats,status='NO_CANDIDATE',goal_satisfied=False),indent=2)+'\n')
            print(json.dumps({'stage':'no_candidate',**stats}),flush=True)
            break
        chosen=candidates[:a.batch]
        (destination/'screen.json').write_text(json.dumps(dict(stats,selected=len(chosen),
            blocks_current_coloring=any(c[0]==6 for c in chosen),
            selected_neighbor_palettes=[c[0] for c in chosen],
            selected_degrees=[c[1] for c in chosen],
            screening_seconds=time.monotonic()-start),indent=2)+'\n')
        print(json.dumps({'stage':'augment','iteration':iteration,'selected':len(chosen),**stats}),flush=True)
        graph=export_graph(points+[c[2] for c in chosen],destination/'graph.json',
                     {'family':'coloring-guided','iteration':iteration,'seed':a.seed,'selected':len(chosen)})
        if a.linear_five_rule is not None or a.screen_linear_five:
            from .color_filter import check_functional, search_functional
            certificate = None
            if a.linear_five_rule is not None:
                functional = json.loads(a.linear_five_rule.read_text())['linear_functional']
                certificate = check_functional(destination/'graph.json', functional, destination/'known-five-rule')
            if a.screen_linear_five and (certificate is None or certificate['status'] != 'FIVE_COLORING_VERIFIED'):
                certificate = search_functional(destination/'graph.json', destination/'five-rule-search')
            if certificate['status'] == 'FIVE_COLORING_VERIFIED':
                print(json.dumps({'stage':'candidate_excluded_by_verified_five_coloring',
                                  'path':str(destination),'goal_satisfied':False}),flush=True)
                break
        if a.native:
            from .solve_external import solve_external
            result=solve_external(destination/'graph.json',destination/'six-color',a.timeout,configuration='sat')
        else:
            result=solve(destination/'graph.json',6,destination/'six-color',a.timeout)
        if result['status']!='SAT':
            print(json.dumps({'stage':'requires_followup','status':result['status'],'path':str(destination)}),flush=True)
            break
        points=[point_type.from_json(v) for v in graph['vertices']]
        colors=json.loads((destination/'six-color'/'coloring.json').read_text())['colors']


if __name__=='__main__':
    main()
