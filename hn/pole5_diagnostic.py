"""Native five-color diagnostic on an exact five-core.

Even an independently verified non-five-colorability result leaves the
requested non-six-colorability target false. The public certificate gate is
not used for this diagnostic.
"""
from __future__ import annotations

import argparse
from collections import deque
import hashlib
import json
from pathlib import Path
import signal
import subprocess
import time

from tools import verify as common
from tools import verify_mq
from .solve_external import DEFAULT_SOLVER,DEFAULT_CHECKER,parse_model,verify_full_cnf_model,compress_completed_run


def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2)+'\n')


def five_core(graph):
    adjacency=[set() for _ in graph['vertices']]
    for u,v in graph['edges']:
        adjacency[u].add(v);adjacency[v].add(u)
    degrees=[len(a) for a in adjacency]
    alive=set(range(len(adjacency)))
    queue=deque(v for v,d in enumerate(degrees) if d<5)
    removed=[]
    removal_degrees=[]
    while queue:
        v=queue.popleft()
        if v not in alive:continue
        alive.remove(v);removed.append(v);removal_degrees.append(degrees[v])
        for w in adjacency[v]:
            if w in alive:
                degrees[w]-=1
                if degrees[w]<5:queue.append(w)
    retained=sorted(alive)
    renumber={v:i for i,v in enumerate(retained)}
    core={k:v for k,v in graph.items() if k not in ('vertices','edges','metadata')}
    core.update(vertices=[graph['vertices'][v] for v in retained],
                edges=[[renumber[u],renumber[v]] for u,v in graph['edges'] if u in alive and v in alive],
                metadata={'family':'exact-five-core-diagnostic','color_test':5,'target_non6_satisfied':False})
    return core,adjacency,retained,removed,removal_degrees


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph',type=Path,default=Path('runs/mq-pole5-v1/graph.json'))
    parser.add_argument('--out',type=Path,default=Path('runs/mq-pole5-five-core-v1'))
    parser.add_argument('--timeout',type=float,default=600)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).read_bytes();(args.out/'program.py').write_bytes(source)
    raw=args.graph.read_bytes();original=json.loads(raw)
    original_hash=hashlib.sha256(raw).hexdigest()
    core,adjacency,retained,removed,removal_degrees=five_core(original)
    graph_path=args.out/'graph.json'
    graph_path.write_text(json.dumps(core,separators=(',',':'))+'\n')
    core_hash=common._sha256(graph_path)
    dump(args.out/'peeling.json',{'full_graph_sha256':original_hash,'core_graph_sha256':core_hash,
                                'core_to_full_vertex_index':retained,'removed_in_order':removed,
                                'degree_at_removal':removal_degrees,
                                'color_equivalence':'Proper five-colorings extend by reversing this elimination order.'})
    geometry=verify_mq.verify_geometry(core)
    cnf=args.out/'coloring.cnf';canonical=common.coloring_cnf(core,5);cnf.write_bytes(canonical)
    proof=args.out/'proof.drat';solution=args.out/'solution.out'
    command=[str(DEFAULT_SOLVER),'--no-colors','--sat','--seed=0','-w',str(solution),str(cnf),str(proof)]
    config={'schema':'hn-native-five-color-diagnostic-v1','k':5,'goal_satisfied':False,'success_certificate':False,
            'full_graph':str(args.graph.resolve()),'full_graph_sha256':original_hash,
            'graph':str(graph_path.resolve()),'graph_sha256':core_hash,
            'vertices':len(core['vertices']),'edges':len(core['edges']),
            'full_vertices':len(original['vertices']),'full_edges':len(original['edges']),
            'removed_vertices':len(removed),'geometry':geometry.__dict__,
            'cnf_sha256':hashlib.sha256(canonical).hexdigest(),'program_sha256':hashlib.sha256(source).hexdigest(),
            'solver_sha256':common._sha256(DEFAULT_SOLVER),'checker_sha256':common._sha256(DEFAULT_CHECKER),
            'command':command,'timeout_seconds':args.timeout,
            'interpretation':'A non-five-colorability proof would establish a lower bound of six only, never the requested seven.'}
    dump(args.out/'config.json',config)
    print(json.dumps({'stage':'native_five_color_started',**config}),flush=True)
    started=time.monotonic();timed_out=False
    with (args.out/'solver.stdout.log').open('wb') as stdout,(args.out/'solver.stderr.log').open('wb') as stderr:
        process=subprocess.Popen(command,stdout=stdout,stderr=stderr,start_new_session=True)
        try:code=process.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out=True;process.send_signal(signal.SIGTERM)
            try:code=process.wait(timeout=5)
            except subprocess.TimeoutExpired:process.kill();code=process.wait()
    result=dict(config,native_exit_code=code,timed_out=timed_out,elapsed_seconds=time.monotonic()-started,
                proof_sha256=common._sha256(proof),proof_bytes=proof.stat().st_size,
                stdout_sha256=common._sha256(args.out/'solver.stdout.log'),
                stderr_sha256=common._sha256(args.out/'solver.stderr.log'))
    if args.graph.read_bytes()!=raw or common._sha256(graph_path)!=core_hash or cnf.read_bytes()!=canonical:
        raise ValueError('A bound graph or formula changed during the diagnostic')
    if timed_out:result['status']='UNKNOWN_TIMEOUT'
    elif code==10:
        assignment=parse_model(solution,5*len(core['vertices']))
        result['cnf_clauses_checked']=verify_full_cnf_model(cnf,assignment)
        colors=[next(c for c in range(5) if assignment[5*v+c+1]) for v in range(len(core['vertices']))]
        witness={'k':5,'colors':colors,'graph_sha256':core_hash}
        common.verify_coloring(core,witness,5);dump(args.out/'coloring.json',witness)
        full_colors=[-1]*len(original['vertices'])
        for v,c in zip(retained,colors):full_colors[v]=c
        for v in reversed(removed):
            forbidden={full_colors[w] for w in adjacency[v] if full_colors[w]>=0}
            full_colors[v]=next(c for c in range(5) if c not in forbidden)
        full_witness={'k':5,'colors':full_colors,'graph_sha256':original_hash}
        full_geometry=verify_mq.verify_geometry(original)
        common.verify_coloring(original,full_witness,5)
        dump(args.out/'full-coloring.json',full_witness)
        dump(args.out/'verification.json',{'goal_satisfied':False,'k':5,'core_geometry':geometry.__dict__,
                                          'full_geometry':full_geometry.__dict__,
                                          'core_graph_sha256':core_hash,'full_graph_sha256':original_hash,
                                          'core_and_full_proper_colorings_verified':True,
                                          'full_coloring_sha256':common._sha256(args.out/'full-coloring.json')})
        result.update(status='SAT',proper_coloring_checked=True,full_coloring_extended_and_verified=True)
    elif code==20:
        result['status']='UNSAT_UNVERIFIED';dump(args.out/'result.json',result)
        checked=common._verify_coloring_proof(core,cnf,proof,DEFAULT_CHECKER,k=5)
        dump(args.out/'verification.json',{'goal_satisfied':False,'success_certificate':False,
                                          'chromatic_lower_bound':6,'geometry':geometry.__dict__,**checked})
        result.update(status='NON5_INDEPENDENTLY_VERIFIED',chromatic_lower_bound=6,
                      goal_satisfied=False,success_certificate=False)
    elif code==0:result['status']='UNKNOWN_NATIVE'
    else:result['status']='ERROR_NATIVE_EXIT'
    dump(args.out/'result.json',result)
    if result['status'] in ('SAT','UNKNOWN_TIMEOUT','UNKNOWN_NATIVE'):
        result=compress_completed_run(args.out)
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
