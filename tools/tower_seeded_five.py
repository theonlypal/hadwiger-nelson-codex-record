#!/usr/bin/env python3
"""Unrestricted five-color SAT with phases from certified rotated seed copies.

The overlap/permutation analysis selects preferences only. It never adds
clauses or assumptions and does not freeze any vertex's color.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from hn.solve_external import solve_external, DEFAULT_SOLVER
from tools.verify_tower import verify_coloring


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')


def run(candidate,seed_coloring,out,timeout):
    candidate=Path(candidate).resolve();out=Path(out).resolve()
    out.mkdir(parents=True,exist_ok=False)
    (out/'program.py').write_bytes(Path(__file__).read_bytes())
    coupling=json.loads((candidate/'coupling.json').read_text())
    seed=Path(coupling['seed']);seed_graph=json.loads(seed.read_text())
    seed_coloring=Path(seed_coloring).resolve();coloring=json.loads(seed_coloring.read_text())
    if digest(seed)!=coupling['seed_sha256'] or coloring['graph_sha256']!=digest(seed):
        raise ValueError('Seed or coloring graph identity mismatch')
    verify_coloring(seed_graph,coloring,5)
    memberships=json.loads((candidate/'copy-memberships.json').read_text())
    original=memberships['original_vertex_ids'];rotated=memberships['rotated_vertex_ids']
    assert len(original)==len(rotated)==len(coloring['colors'])
    assert len(set(original))==len(original) and len(set(rotated))==len(rotated)
    original_index={v:i for i,v in enumerate(original)}
    rotated_index={v:i for i,v in enumerate(rotated)}
    shared=sorted(set(original_index)&set(rotated_index))
    overlap=[{'vertex':v,'original_seed_vertex':original_index[v],
              'rotated_seed_vertex':rotated_index[v],
              'original_color':coloring['colors'][original_index[v]],
              'rotated_color_before_permutation':coloring['colors'][rotated_index[v]]}
             for v in shared]
    counts=[[0]*5 for _ in range(5)]
    for row in overlap:counts[row['rotated_color_before_permutation']][row['original_color']]+=1
    scored=[(sum(counts[c][permutation[c]] for c in range(5)),permutation)
            for permutation in itertools.permutations(range(5))]
    best_matches,best_permutation=max(scored)
    compatible=[permutation for matches,permutation in scored if matches==len(shared)]
    hinted={v:coloring['colors'][i] for i,v in enumerate(original)}
    selection='original complete copy only; shared colors admit no global permutation'
    chosen=None;cross_conflicts=None
    if compatible:
        choices=[]
        for permutation in compatible:
            full=dict(hinted)
            for i,v in enumerate(rotated):
                c=permutation[coloring['colors'][i]]
                assert v not in full or full[v]==c
                full[v]=c
            conflicts=sum(full[a]==full[b] for a,b in memberships['strict_cross_copy_edges'])
            choices.append((conflicts,permutation,full))
        cross_conflicts,chosen,hinted=min(choices,key=lambda row:(row[0],row[1]))
        selection='both complete copies; choose overlap-compatible permutation with fewest strict cross-copy conflicts'
    phases=[]
    for vertex,color in sorted(hinted.items()):
        for c in range(5):
            variable=5*vertex+c+1;phases.append(variable if c==color else -variable)
    phase_path=out/'phases.txt'
    phase_path.write_text(' '.join(map(str,phases))+' 0\n')
    record={'scope':'search phase preferences only; no restrictions added to CNF',
            'candidate_graph':str(candidate/'graph.json'),'candidate_graph_sha256':digest(candidate/'graph.json'),
            'seed_graph':str(seed),'seed_graph_sha256':digest(seed),
            'seed_coloring':str(seed_coloring),'seed_coloring_sha256':digest(seed_coloring),
            'copy_memberships_sha256':digest(candidate/'copy-memberships.json'),
            'shared_vertices':len(shared),'shared_vertices_details':overlap,
            'overlap_color_counts_rotated_by_original':counts,
            'permutations_examined':len(scored),'compatible_permutations':compatible,
            'maximum_shared_matches':best_matches,'best_permutation':best_permutation,
            'chosen_permutation':chosen,'selection':selection,
            'hinted_vertices':len(hinted),'hinted_variables':len(phases),
            'strict_cross_copy_conflicts_in_complete_hint':cross_conflicts,
            'phase_file_sha256':digest(phase_path),'success_certificate':False}
    dump(out/'overlap-and-phases.json',record)
    print(json.dumps({k:v for k,v in record.items() if k!='shared_vertices_details'}),flush=True)
    return solve_external(candidate/'graph.json',out/'native',timeout=timeout,
                          solver=DEFAULT_SOLVER.with_name('cadical-phase'),
                          configuration='sat',k=5,phase_file=phase_path)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,default=Path('runs/mixed-g3-v31-odd-rotation-v1'))
    parser.add_argument('--seed-coloring',type=Path,default=Path('runs/mixed-haugland-g3-v31-linear5/coloring.json'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--timeout',type=float,default=300)
    args=parser.parse_args()
    if not 0<args.timeout<=300:parser.error('timeout must lie in (0,300]')
    run(args.candidate,args.seed_coloring,args.out,args.timeout)
