"""One exact pole-5 extension that defeats a specified five-color model."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from .adaptive import candidate_pool
from .multiquadratic import Point
from .mq_construct import THETA4,export,canonical_points
from .mq_intersections import core_indices
from .mq_pole5_search import pole5
from .solve_external import solve_external


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph',type=Path,default=Path('runs/mq-pole5-v1/graph.json'))
    parser.add_argument('--coloring',type=Path,default=Path('runs/mq-pole5-five-core-v1/full-coloring.json'))
    parser.add_argument('--out',type=Path,default=Path('runs/mq-pole5-adaptive-v1'))
    parser.add_argument('--max-new',type=int,default=1000)
    parser.add_argument('--screen-cap',type=int,default=6000)
    parser.add_argument('--timeout',type=float,default=300)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).read_bytes();(args.out/'program.py').write_bytes(source)
    graph=json.loads(args.graph.read_text());witness=json.loads(args.coloring.read_text())
    full=[Point.from_json(v) for v in graph['vertices']];full_set=set(full)
    colors=witness['colors']
    if witness['k']!=5 or len(colors)!=len(full) or any(colors[u]==colors[v] for u,v in graph['edges']):
        raise ValueError('Expected a proper full five-coloring')
    if witness['graph_sha256']!=hashlib.sha256(args.graph.read_bytes()).hexdigest():
        raise ValueError('Coloring graph hash differs')
    core=sorted(core_indices(graph,6));points=[full[v] for v in core];core_colors=[colors[v] for v in core]
    frequencies=Counter()
    for u,v in graph['edges']:
        d=full[v]-full[u];frequencies[d]+=1;frequencies[-d]+=1
    directions=tuple(d for d,_ in frequencies.most_common())
    assert all(d.is_unit() for d in directions)
    direction_set=set(directions)
    theta4_companions=[d for d in directions if d*THETA4 in direction_set]
    heule_report={'directed_unit_vectors':len(directions),'theta4_companion_count':len(theta4_companions),
                  'no_isometric_published_Heule_realization':not theta4_companions,
                  'argument':'The published Heule realization contains unit rays1 andtheta4 from its origin. Every rotated or reflected copy therefore requires two graph-directed unit vectors differing bytheta4. None occur in the complete stored direction set.',
                  'limitation':'This excludes isometric copies of the published coordinate realization, not all abstract Heule-subgraph embeddings or all five-chromatic subgraphs.',
                  'unit_vectors':[d.to_json() for d in directions]}
    (args.out/'heule-pattern-audit.json').write_text(json.dumps(heule_report,separators=(',',':'))+'\n')
    config={'source_graph':str(args.graph),'source_graph_sha256':witness['graph_sha256'],
            'source_coloring':str(args.coloring),'source_coloring_sha256':hashlib.sha256(args.coloring.read_bytes()).hexdigest(),
            'program_sha256':hashlib.sha256(source).hexdigest(),'max_new':args.max_new,
            'screen_cap':args.screen_cap,'directions':len(directions),'screen_vertices':len(core),
            'selection':'New exact centers with at least6 old-six-core neighbors covering all5 saved colors, including both pole5 and5-integral neighbors. Rank balanced cross-field support then degree.',
            'native_six_timeout_seconds':args.timeout,'goal_satisfied':False}
    (args.out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    print(json.dumps({'stage':'screen_started',**config}),flush=True)
    start=time.monotonic()
    pool,stats=candidate_pool(points,core_colors,directions,max_proposals=args.screen_cap)
    qualified=[]
    for palette_size,degree,p,neighbors in pool:
        if palette_size!=5 or p in full_set:continue
        actual_palette={core_colors[v] for v in neighbors}
        assert actual_palette==set(range(5))
        poles=sum(pole5(points[v]) for v in neighbors)
        integral=len(neighbors)-poles
        if not poles or not integral:continue
        qualified.append({'point':p,'neighbors':[core[v] for v in neighbors],
                          'neighbor_colors':[core_colors[v] for v in neighbors],
                          'pole_neighbors':poles,'integral_neighbors':integral})
    qualified.sort(key=lambda p:(min(p['pole_neighbors'],p['integral_neighbors']),len(p['neighbors'])),reverse=True)
    chosen=qualified[:args.max_new]
    report=dict(stats,qualified_mixed_five_palette_centers=len(qualified),selected=len(chosen),
                screening_seconds=time.monotonic()-start,
                chosen_degrees=[len(p['neighbors']) for p in chosen],
                minimum_known_six_core_degree=min((len(p['neighbors']) for p in chosen),default=0),
                every_new_point_blocks_saved_five_coloring=True,
                completeness_claimed=False,
                selection_provenance=[dict(p,point=p['point'].to_json(),
                                           generator_origin=p['neighbors'][0],
                                           generator_unit_vector=(p['point']-full[p['neighbors'][0]]).to_json()) for p in chosen])
    (args.out/'screen.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('selection_provenance','chosen_degrees')}),flush=True)
    if not chosen:return
    path=args.out/'graph.json'
    expanded=export(full+[p['point'] for p in chosen],path,
                    dict(config,family='pole5-five-model-blocking',selected=len(chosen)))
    canonical=canonical_points(full+[p['point'] for p in chosen]);index={p:i for i,p in enumerate(canonical)}
    edge_set={tuple(e) for e in expanded['edges']}
    for item in chosen:
        v=index[item['point']]
        for neighbor in item['neighbors']:
            u=index[full[neighbor]]
            assert (min(u,v),max(u,v)) in edge_set
    enlarged_core=core_indices(expanded,6)
    new_ids={index[p['point']] for p in chosen}
    assert new_ids<=enlarged_core
    summary={'vertices':len(canonical),'edges':len(expanded['edges']),
             'new_vertices':len(chosen),'six_core_vertices':len(enlarged_core),
             'six_core_edges':sum(u in enlarged_core and v in enlarged_core for u,v in expanded['edges']),
             'new_vertices_in_six_core':len(new_ids&enlarged_core),
             'every_new_vertex_blocks_saved_five_coloring':True,'goal_satisfied':False}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'stage':'core_verified',**summary}),flush=True)
    result=solve_external(path,args.out/'native-six-color',timeout=args.timeout,configuration='sat')
    summary['native_six_status']=result['status'];summary['goal_satisfied']=result.get('goal_satisfied',False)
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
