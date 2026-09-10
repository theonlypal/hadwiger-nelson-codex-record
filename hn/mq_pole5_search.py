"""Exact circle intersections that can introduce rational poles at five.

All center pairs are screened before orbit selection. Each accepted pair has
one of five exact squared separations; the proposed common unit neighbors
are then checked algebraically. Selection keeps complete 60-degree orbits.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from .multiquadratic import MQ, Point, SQRT5, SQRT3, SQRT11
from .mq_construct import HEX_ROTATION, canonical_points, export
from .mq_intersections import core_indices
from .solve_external import solve_external

# One half of sqrt((4-s)/s), the multiplier of i*(b-a).
OFFSETS = {Fraction(5,4):SQRT5*SQRT11/10,
           Fraction(5,2):SQRT3*SQRT5/10,
           Fraction(10,3):SQRT5/10,
           Fraction(20,9):SQRT5/5,
           Fraction(25,9):SQRT11/10}


def pole5(point):
    return point.x.denominator%5==0 or point.y.denominator%5==0


def scan(points, candidate_cap):
    index = {p:i for i,p in enumerate(points)}
    coordinates = np.array([[p.x.approximate(),p.y.approximate()] for p in points])
    candidates = {}
    pair_counts = Counter()
    rejected = 0
    for a,p in enumerate(points):
        differences = coordinates[a+1:]-coordinates[a]
        squared = np.einsum('ij,ij->i',differences,differences)
        for separation,multiplier in OFFSETS.items():
            for offset in np.flatnonzero(np.abs(squared-float(separation))<1e-8):
                b = a+1+int(offset)
                delta = points[b]-p
                if delta.norm_squared()!=MQ(separation):
                    rejected+=1
                    continue
                pair_counts[str(separation)]+=1
                midpoint = (p+points[b])/2
                displacement = Point(-delta.y,delta.x)*multiplier
                for sign in (1,-1):
                    c = midpoint+sign*displacement
                    assert (c-p).is_unit() and (c-points[b]).is_unit()
                    if c in index:
                        continue
                    if c not in candidates:
                        if len(candidates)>=candidate_cap:
                            return candidates,{'candidate_cap_reached':True,'center_pairs':dict(pair_counts),
                                               'last_center_index':a,'approximate_pairs_rejected':rejected}
                        candidates[c]={'parents':[a,b],'separation_squared':str(separation),
                                       'sign':sign,'known_old_neighbors':set()}
                    candidates[c]['known_old_neighbors'].update((a,b))
        if a and a%2000==0:
            print(json.dumps({'stage':'pair_scan','old_centers_processed':a,
                              'new_centers':len(candidates),'center_pairs':dict(pair_counts)}),flush=True)
    return candidates,{'candidate_cap_reached':False,'center_pairs':{str(s):pair_counts[str(s)] for s in OFFSETS},
                       'approximate_pairs_rejected':rejected,'old_centers_processed':len(points)}


def select_orbits(points,candidates,vertex_cap):
    index = {p:i for i,p in enumerate(points)}
    rotations = [HEX_ROTATION**j for j in range(6)]
    if any(HEX_ROTATION*p not in index for p in points):
        raise ValueError('Seed does not have exact 60-degree rotational symmetry')
    unassigned = set(candidates)
    orbits = []
    while unassigned:
        representative = next(iter(unassigned))
        members = canonical_points(r*representative for r in rotations)
        source = candidates[representative]
        for j,r in enumerate(rotations):
            member = r*representative
            if member not in candidates:
                # Rotating both verified parent circles certifies the missing
                # screened representative without any completeness assumption.
                parents = [index[r*points[v]] for v in source['parents']]
                candidates[member]=dict(source,parents=parents,known_old_neighbors=set(parents),
                                        rotation_derived=True)
            unassigned.discard(member)
        degrees = [len(candidates[p]['known_old_neighbors']) for p in members]
        radius = members[0].norm_squared().approximate()
        orbit_key = tuple((p.x.denominator,p.x.coefficients,p.y.denominator,p.y.coefficients) for p in members)
        orbits.append({'members':members,'has_pole5':any(map(pole5,members)),
                       'minimum_known_old_degree':min(degrees),'known_old_degree_sum':sum(degrees),
                       'squared_radius_display':radius,'sort_key':orbit_key})
    orbits.sort(key=lambda o:(-int(o['has_pole5']),-o['minimum_known_old_degree'],
                              -o['known_old_degree_sum'],-len(o['members']),
                              o['squared_radius_display'],o['sort_key']))
    selected = []
    current = len(points)
    for orbit in orbits:
        if current+len(orbit['members'])<=vertex_cap:
            selected.append(orbit)
            current+=len(orbit['members'])
    return selected,{'candidate_orbits':len(orbits),'selected_orbits':len(selected),
                     'candidate_pole5_orbits':sum(o['has_pole5'] for o in orbits),
                     'selected_pole5_orbits':sum(o['has_pole5'] for o in selected),
                     'selection_truncated':len(selected)<len(orbits),
                     'selection_rule':'Complete sixfold orbits; prioritize actual pole5, minimum known old degree, total known old degree, then radius and canonical coordinates'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=Path,default=Path('runs/mq-intersections-v1939-r1/graph.json'))
    parser.add_argument('--out',type=Path,default=Path('runs/mq-pole5-v1'))
    parser.add_argument('--max-vertices',type=int,default=40000)
    parser.add_argument('--candidate-cap',type=int,default=300000)
    parser.add_argument('--timeout',type=float,default=300)
    args = parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    source = Path(__file__).read_bytes();(args.out/'program.py').write_bytes(source)
    old_graph = json.loads(args.seed.read_text())
    points = [Point.from_json(v) for v in old_graph['vertices']]
    if args.max_vertices<=len(points):
        raise ValueError('Vertex cap must leave space beyond the seed')
    config = {'seed':str(args.seed),'seed_sha256':hashlib.sha256(args.seed.read_bytes()).hexdigest(),
              'program_sha256':hashlib.sha256(source).hexdigest(),
              'separations_squared':list(map(str,OFFSETS)),
              'max_vertices':args.max_vertices,'candidate_cap':args.candidate_cap,
              'native_solver_timeout_seconds':args.timeout,'goal_satisfied':False,
              'candidate_pair_filter':'float squared distance within1e-8; exact equality and exact two unit edges for acceptance; no completeness claim'}
    (args.out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    started = time.monotonic()
    candidates,scan_stats = scan(points,args.candidate_cap)
    scan_stats.update({'new_centers':len(candidates),'new_centers_with_pole5':sum(map(pole5,candidates)),
                       'known_old_degree_histogram':dict(sorted(Counter(len(v['known_old_neighbors']) for v in candidates.values()).items())),
                       'elapsed_seconds':time.monotonic()-started})
    (args.out/'scan.json').write_text(json.dumps(scan_stats,indent=2)+'\n')
    print(json.dumps({'stage':'scan_complete',**scan_stats}),flush=True)
    if scan_stats['candidate_cap_reached']:
        print(json.dumps({'status':'CANDIDATE_CAP_REACHED','graph_not_exported':True,'goal_satisfied':False}),flush=True)
        return
    if not scan_stats['new_centers_with_pole5']:
        print(json.dumps({'status':'NO_GENUINE_POLE5','goal_satisfied':False}),flush=True)
        return
    selected,selection_stats = select_orbits(points,candidates,args.max_vertices)
    additions = [p for orbit in selected for p in orbit['members']]
    print(json.dumps({'stage':'orbit_selection',**selection_stats,'new_vertices':len(additions)}),flush=True)
    graph_path = args.out/'graph.json'
    graph = export(points+additions,graph_path,dict(config,scan=scan_stats,selection=selection_stats,
                                                    family='denominator5-common-unit-neighbor-closure'))
    canonical = canonical_points(points+additions)
    graph_index = {p:i for i,p in enumerate(canonical)}
    provenance = {'seed_to_graph_index':[graph_index[p] for p in points],
                  'graph_sha256':hashlib.sha256(graph_path.read_bytes()).hexdigest(),
                  'orbits':[]}
    for orbit in selected:
        provenance['orbits'].append({k:v for k,v in orbit.items() if k not in ('members','sort_key')})
        provenance['orbits'][-1]['members']=[dict(graph_index=graph_index[p],
                                                **{k:sorted(v) if isinstance(v,set) else v
                                                   for k,v in candidates[p].items()}) for p in orbit['members']]
    (args.out/'provenance.json').write_text(json.dumps(provenance,separators=(',',':'))+'\n')
    core = core_indices(graph,6)
    old_ids = {graph_index[p] for p in points}
    poles = {i for i,p in enumerate(canonical) if pole5(p)}
    crossing = [(u,v) for u,v in graph['edges'] if (u in poles)!=(v in poles)]
    summary = {'vertices':len(canonical),'edges':len(graph['edges']),
               'new_vertices':len(additions),'new_vertices_in_six_core':len(core-old_ids),
               'six_core_vertices':len(core),'six_core_edges':sum(u in core and v in core for u,v in graph['edges']),
               'pole5_vertices':len(poles),'pole5_vertices_in_six_core':len(poles&core),
               'pole5_cross_edges':len(crossing),
               'pole5_cross_edges_in_six_core':sum(u in core and v in core for u,v in crossing),
               'selection':selection_stats,'scan':scan_stats,'goal_satisfied':False}
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'stage':'six_core',**summary}),flush=True)
    result = solve_external(graph_path,args.out/'native-six-color',timeout=args.timeout)
    summary['native_status']=result['status']
    summary['goal_satisfied']=result.get('goal_satisfied',False)
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    main()
