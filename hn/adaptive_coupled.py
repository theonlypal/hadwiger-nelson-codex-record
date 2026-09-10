"""Coloring-guided centers with neighbors on both sides of a field cut."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from .field import Element
from .construct import spindle_rotation
from .search import export
from .solve import solve


def extension_label(point):
    return int(any(point.coefficients[len(point.coefficients)//2:]))


def cross_stats(graph, points, excluded_hub, selected=None):
    labels = [extension_label(p) for p in points]
    population = range(len(points)) if selected is None else selected
    crossing = [(u,v) for u,v in graph['edges'] if labels[u] != labels[v]]
    counts = Counter(v for edge in crossing for v in edge)
    without_hub = [(u,v) for u,v in crossing if points[u] != excluded_hub and points[v] != excluded_hub]
    return {'partition_sizes':[sum(labels[v]==i for v in population) for i in (0,1)],
            'cross_edges':len(crossing),'cross_endpoints':len(counts),
            'largest_cross_degree':max(counts.values(),default=0),
            'cross_edges_excluding_seed_hub':len(without_hub),
            'cross_endpoints_excluding_seed_hub':len({v for edge in without_hub for v in edge})}


def candidate_pool(points, colors, directions, excluded_hub, max_proposals=10000,
                   minimum_existing_degree=6):
    """Approximate grouping proposes; exact lookup decides qualification."""
    coordinate = np.array([[p.approximate().real,p.approximate().imag] for p in points])
    offsets = np.array([[p.approximate().real,p.approximate().imag] for p in directions])
    table = np.rint((coordinate[:,None,:]+offsets[None,:,:]).reshape(-1,2)*1e7).astype(np.int64)
    keys = np.ascontiguousarray(table).view(np.dtype((np.void,16))).ravel()
    unique,first,inverse,counts = np.unique(keys,return_index=True,return_inverse=True,return_counts=True)
    color_masks = np.zeros(len(unique),dtype=np.uint8)
    color_bits = np.repeat(np.left_shift(np.uint8(1),np.array(colors,dtype=np.uint8)),len(directions))
    np.bitwise_or.at(color_masks,inverse,color_bits)
    sides = np.array([0 if p==excluded_hub else 1+extension_label(p) for p in points],dtype=np.uint8)
    repeated_sides = np.repeat(sides,len(directions))
    side_a = np.bincount(inverse,weights=(repeated_sides==1),minlength=len(unique)).astype(np.int32)
    side_b = np.bincount(inverse,weights=(repeated_sides==2),minlength=len(unique)).astype(np.int32)
    popcounts = np.array([j.bit_count() for j in range(256)],dtype=np.uint8)[color_masks]
    mixed = np.flatnonzero((side_a>0)&(side_b>0))
    qualified = np.flatnonzero((counts>=6)&(popcounts>=4)&(side_a>0)&(side_b>0))
    eligible = mixed
    scores = ((popcounts[eligible]==6).astype(np.int64)*10**12
              +np.minimum(side_a[eligible],side_b[eligible]).astype(np.int64)*10**8
              +popcounts[eligible].astype(np.int64)*10**6+counts[eligible])
    ranked = eligible[np.argsort(scores)[::-1]][:max_proposals]
    index = {p:i for i,p in enumerate(points)}
    proposals, seen = [], set()
    existing = 0
    exact_degrees = Counter()
    exact_palette_sizes = Counter()
    for group in ranked:
        origin,direction = divmod(int(first[group]),len(directions))
        point = points[origin]+directions[direction]
        if point in index:
            existing += 1
            continue
        if point in seen:
            continue
        seen.add(point)
        neighbors = set()
        for unit in directions:
            neighbor = index.get(point-unit)
            if neighbor is not None:
                neighbors.add(neighbor)
        a = [v for v in neighbors if points[v]!=excluded_hub and extension_label(points[v])==0]
        b = [v for v in neighbors if points[v]!=excluded_hub and extension_label(points[v])==1]
        palette = {colors[v] for v in neighbors}
        if a and b:
            exact_degrees[len(neighbors)] += 1
            exact_palette_sizes[len(palette)] += 1
        if len(neighbors)>=minimum_existing_degree and a and b:
            proposals.append({'point':point,'neighbors':sorted(neighbors),'palette':sorted(palette),
                              'base_neighbors_excluding_hub':len(a),
                              'extension_neighbors_excluding_hub':len(b),
                              'origin':origin,'direction':direction})
    proposals.sort(key=lambda p:(len(p['palette'])==6,
                                min(p['base_neighbors_excluding_hub'],p['extension_neighbors_excluding_hub']),
                                len(p['palette']),len(p['neighbors'])),reverse=True)
    return proposals, {'approximate_groups':len(unique),'approximate_eligible_groups':len(qualified),
                       'approximate_mixed_groups':len(mixed),
                       'mixed_groups_examined':len(ranked),
                       'mixed_screen_truncated':len(mixed)>max_proposals,
                       'existing_points_excluded':existing,
                       'new_exact_mixed_degree_histogram':dict(sorted(exact_degrees.items())),
                       'new_exact_mixed_palette_histogram':dict(sorted(exact_palette_sizes.items())),
                       'exact_proposals':len(proposals),
                       'palette_histogram':dict(Counter(len(p['palette']) for p in proposals))}


def joint_experiment(graph, points, colors, directions, excluded_hub, out, timeout):
    """Retain the complete screened mixed-center batch, then peel jointly."""
    from collections import deque
    from tools.core_audit import core_numbers
    proposals,stats = candidate_pool(points,colors,directions,excluded_hub,
                                    minimum_existing_degree=2)
    if stats['mixed_screen_truncated']:
        raise RuntimeError('Joint screen is truncated; refusing an ALL-centers claim')
    old_colors = dict(zip(points,colors))
    old_points = set(points)
    (out/'centers.json').write_text(json.dumps({
        'stats':stats,
        'centers':[dict(p,point=p['point'].to_json(),
                        exact_direction=directions[p['direction']].to_json()) for p in proposals]
    },separators=(',',':'))+'\n')
    print(json.dumps({'stage':'joint_centers','new_centers':len(proposals),**stats}),flush=True)
    graph_path = out/'graph.json'
    expanded = export(points+[p['point'] for p in proposals],graph_path,
                      {'family':'joint-cross-field-centers','seed_vertices':len(points),
                       'new_centers':len(proposals),'all_screened_mixed_centers':True,
                       'all_old_degree_two_centers_retained_before_peeling':True})
    new_points = [Element.from_json(p) for p in expanded['vertices']]
    adjacency = [set() for _ in new_points]
    for u,v in expanded['edges']:
        adjacency[u].add(v);adjacency[v].add(u)
    cores = core_numbers(adjacency)
    core = {v for v,c in enumerate(cores) if c>=6}
    old_ids = {v for v,p in enumerate(new_points) if p in old_points}
    core_graph = dict(expanded,edges=[e for e in expanded['edges'] if e[0] in core and e[1] in core])
    result = {'vertices':len(new_points),'edges':len(expanded['edges']),
              'degeneracy':max(cores,default=0),'six_core_vertices':len(core),
              'six_core_edges':len(core_graph['edges']),
              'new_vertices_surviving_six_core':len(core-old_ids),
              'new_vertices_removed':len(new_points)-len(old_points)-len(core-old_ids),
              'old_vertices_retained':len(core&old_ids),
              'new_vertices_original_degree_histogram':dict(Counter(len(adjacency[v]) for v in range(len(new_points)) if v not in old_ids)),
              'coupling_full':cross_stats(expanded,new_points,excluded_hub),
              'coupling_core':cross_stats(core_graph,new_points,excluded_hub,core),
              'core_collapses_to_seed':core==old_ids,'goal_satisfied':False}
    (out/'core.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'stage':'joint_core',**result}),flush=True)
    color_path = out/'six-color'/'coloring.json'
    if core==old_ids:
        degree = [len(a) for a in adjacency]
        alive = set(range(len(new_points)))
        queue = deque(v for v,d in enumerate(degree) if d<6)
        removed = []
        while queue:
            v = queue.popleft()
            if v not in alive:
                continue
            alive.remove(v);removed.append(v)
            for w in adjacency[v]:
                if w in alive:
                    degree[w]-=1
                    if degree[w]<6:queue.append(w)
        assert alive==old_ids
        extended_colors = [old_colors.get(p,-1) for p in new_points]
        for v in reversed(removed):
            forbidden = {extended_colors[w] for w in adjacency[v] if extended_colors[w]>=0}
            extended_colors[v] = next(c for c in range(6) if c not in forbidden)
        assert all(extended_colors[u]!=extended_colors[v] for u,v in expanded['edges'])
        color_path.parent.mkdir()
        sha = hashlib.sha256(graph_path.read_bytes()).hexdigest()
        color_path.write_text(json.dumps({'k':6,'colors':extended_colors,'graph_sha256':sha})+'\n')
        (color_path.parent/'result.json').write_text(json.dumps({
            'status':'SAT','method':'reverse-six-core-peeling','solver':None,
            'graph_sha256':sha,'k':6,'vertices':len(new_points),'edges':len(expanded['edges']),
            'proper_coloring_checked':True,'success_certificate':False},indent=2)+'\n')
        (out/'elimination.json').write_text(json.dumps({'removed_in_order':removed,
                                                      'remaining_seed_vertex_ids':sorted(alive)})+'\n')
    else:
        answer = solve(graph_path,6,color_path.parent,timeout)
        if answer['status']!='SAT':
            return
    verify = subprocess.run([sys.executable,'tools/verify.py','--graph',str(graph_path),
                             '--k','6','--coloring',str(color_path)],capture_output=True,text=True)
    (out/'independent-verification.json').write_text(verify.stdout)
    if verify.returncode:
        (out/'independent-verification.stderr').write_text(verify.stderr)
        raise RuntimeError('Independent verification failed')
    print(json.dumps({'stage':'independently_verified_proper_six_coloring',
                      'solver_used':not result['core_collapses_to_seed']}),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',default='artifacts/haugland_g3.json')
    parser.add_argument('--coloring',default='runs/baseline-g3-k6/coloring.json')
    parser.add_argument('--out',default='runs/adaptive-g3-coupled-v1')
    parser.add_argument('--rounds',type=int,default=10)
    parser.add_argument('--batch',type=int,default=50)
    parser.add_argument('--timeout',type=float,default=60)
    parser.add_argument('--joint',action='store_true',help='Add all mixed centers jointly before six-core peeling')
    args = parser.parse_args()
    out = Path(args.out);out.mkdir(parents=True,exist_ok=False)
    (out/'program.py').write_bytes(Path(__file__).read_bytes())
    graph = json.loads(Path(args.seed).read_text())
    if graph.get('schema')!='hn-exact-v1' or graph['field']['cyclotomic_order']!=42 or graph['field']['radicand']!=5:
        raise ValueError('This experiment requires the specified Haugland field')
    points = [Element.from_json(v) for v in graph['vertices']]
    witness = json.loads(Path(args.coloring).read_text());colors = witness['colors']
    if witness['k']!=6 or len(colors)!=len(points) or any(colors[u]==colors[v] for u,v in graph['edges']):
        raise ValueError('Initial six-coloring is not proper')
    crosses = Counter(v for edge in graph['edges'] if extension_label(points[edge[0]])!=extension_label(points[edge[1]]) for v in edge)
    excluded_hub = points[crosses.most_common(1)[0][0]]
    frequency = Counter()
    for u,v in graph['edges']:
        step = points[v]-points[u]
        frequency[step]+=1;frequency[-step]+=1
    original_directions = tuple(p for p,_ in frequency.most_common())
    if not all(p.is_unit() for p in original_directions):
        raise ValueError('Seed contains a nonunit edge')
    config = dict(vars(args),program_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  seed_sha256=hashlib.sha256(Path(args.seed).read_bytes()).hexdigest(),
                  initial_directions=len(original_directions),excluded_hub=excluded_hub.to_json(),
                  selection=('Retain all screened mixed centers jointly before six-core peeling' if args.joint else
                             'Exact degree>=6, neighbors on both field sides excluding seed hub; rank six-palette blocking, minimum side support, palette count, degree'),
                  direction_expansion='D union rho*D union rho^2*D if original D yields no qualifying center',
                  initial_coupling=cross_stats(graph,points,excluded_hub),goal_satisfied=False)
    (out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    if args.joint:
        joint_experiment(graph,points,colors,original_directions,excluded_hub,out,args.timeout)
        return
    directions = original_directions
    expanded = False
    previous_graph = str(Path(args.seed).resolve())
    for iteration in range(args.rounds):
        destination = out/f'round-{iteration:03d}';destination.mkdir()
        start = time.monotonic()
        proposals,stats = candidate_pool(points,colors,directions,excluded_hub)
        if not proposals and not expanded:
            (destination/'original-directions-screen.json').write_text(json.dumps(stats,indent=2)+'\n')
            rho = spindle_rotation()
            directions = tuple(sorted({rho**j*d for j in range(3) for d in original_directions},
                                      key=lambda d:(d.denominator,d.coefficients)))
            assert all(d.is_unit() for d in directions)
            expanded = True
            proposals,stats = candidate_pool(points,colors,directions,excluded_hub)
        chosen = proposals[:args.batch]
        screen = dict(stats,iteration=iteration,directions=len(directions),expanded=expanded,
                      selected=len(chosen),previous_graph=previous_graph,
                      coupling_before=cross_stats(graph,points,excluded_hub),
                      selected_provenance=[dict(p,point=p['point'].to_json(),
                                                exact_direction=directions[p['direction']].to_json()) for p in chosen],
                      screening_seconds=time.monotonic()-start,goal_satisfied=False)
        (destination/'screen.json').write_text(json.dumps(screen,indent=2)+'\n')
        print(json.dumps({k:v for k,v in screen.items() if k!='selected_provenance'}),flush=True)
        if not chosen:
            print(json.dumps({'status':'NO_QUALIFYING_CENTER','goal_satisfied':False}),flush=True)
            break
        graph_path = destination/'graph.json'
        graph = export(points+[p['point'] for p in chosen],graph_path,
                       {'family':'cross-field-coloring-guided','iteration':iteration,
                        'seed':args.seed,'selected':len(chosen),'expanded_directions':expanded})
        points = [Element.from_json(v) for v in graph['vertices']]
        after = cross_stats(graph,points,excluded_hub)
        (destination/'coupling.json').write_text(json.dumps(after,indent=2)+'\n')
        print(json.dumps({'stage':'coupling_after',**after}),flush=True)
        result = solve(graph_path,6,destination/'six-color',args.timeout)
        if result['status']!='SAT':
            break
        colors = json.loads((destination/'six-color'/'coloring.json').read_text())['colors']
        previous_graph = str(graph_path.resolve())


if __name__=='__main__':
    main()
