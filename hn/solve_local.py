"""Bounded local coloring search; only independently checked proper models count."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

from tools import verify


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph',type=Path,default=Path('runs/odd-power-t6-spindle-v1/graph.json'))
    parser.add_argument('--old-graph',type=Path,default=Path('runs/t6-spindle/graph.json'))
    parser.add_argument('--source-coloring',type=Path,default=Path('runs/t6-spindle/six-color/coloring.json'))
    parser.add_argument('--overlap',type=Path,default=Path('runs/odd-power-t6-phased-v1/overlap.json'))
    parser.add_argument('--out',type=Path,default=Path('runs/odd-power-t6-tabu-seed0-v1'))
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--seconds',type=float,default=120)
    parser.add_argument('--max-moves',type=int,default=2_000_000)
    parser.add_argument('--sample',type=int,default=64)
    args = parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    raw=args.graph.read_bytes();graph=json.loads(raw);old=json.loads(args.old_graph.read_text())
    witness=json.loads(args.source_coloring.read_text());overlap=json.loads(args.overlap.read_text())
    graph_hash=hashlib.sha256(raw).hexdigest()
    if overlap['new_graph_sha256']!=graph_hash or overlap['old_graph_sha256']!=verify._sha256(args.old_graph):
        raise ValueError('Overlap graph hashes differ')
    if overlap['source_color_sha256']!=verify._sha256(args.source_coloring):
        raise ValueError('Overlap coloring hash differs')
    verify.verify_coloring(old,witness,6)
    preferences=[-1]*len(graph['vertices']);seen_old=set();seen_new=set()
    for u,v in overlap['mapping']:
        if u in seen_old or v in seen_new or old['vertices'][u]!=graph['vertices'][v]:
            raise ValueError('Overlap is not an exact injective coordinate match')
        seen_old.add(u);seen_new.add(v);preferences[v]=witness['colors'][u]
    source=Path('tools/tabu_color.cpp');snapshot=args.out/'tabu_color.cpp';snapshot.write_bytes(source.read_bytes())
    (args.out/'program.py').write_bytes(Path(__file__).read_bytes())
    compiler=Path(shutil.which('clang++') or shutil.which('g++') or '').resolve(strict=True)
    binary=(args.out/'tabu-color').resolve()
    compile_command=[str(compiler),'-std=c++17','-O3','-DNDEBUG',str(snapshot),'-o',str(binary)]
    compiled=subprocess.run(compile_command,capture_output=True)
    (args.out/'compiler.stdout.log').write_bytes(compiled.stdout);(args.out/'compiler.stderr.log').write_bytes(compiled.stderr)
    if compiled.returncode:raise RuntimeError('Candidate generator compilation failed')
    compiler_version=subprocess.run([str(compiler),'--version'],capture_output=True,check=True).stdout.decode()
    input_path=args.out/'graph.edges';init=args.out/'initial-preferences.txt'
    input_path.write_text(f"{len(graph['vertices'])} {len(graph['edges'])} 6\n"+
                          ''.join(f'{u} {v}\n' for u,v in graph['edges']))
    init.write_text(''.join(f'{c}\n' for c in preferences))
    best=args.out/'best-model.txt';stats=args.out/'native-result.json'
    command=[str(binary),str(input_path),str(init),str(best),str(stats),str(args.seed),
             str(args.seconds),str(args.max_moves),str(args.sample)]
    config={'graph_sha256':graph_hash,'graph':str(args.graph.resolve()),'vertices':len(graph['vertices']),
            'edges':len(graph['edges']),'k':6,'seed':args.seed,'seconds':args.seconds,'max_moves':args.max_moves,
            'sample_conflicted_vertices':args.sample,'source_preferences':len(seen_new),
            'initialization':'Exact old-color overlap as preferences; remaining vertices greedy by degree. All vertices free afterward.',
            'scope':'Untrusted heuristic candidate generator. Failure establishes no lower bound.',
            'overlap_sha256':verify._sha256(args.overlap),'source_coloring_sha256':verify._sha256(args.source_coloring),
            'program_sha256':verify._sha256(Path(__file__)),'cpp_sha256':verify._sha256(snapshot),
            'binary_sha256':verify._sha256(binary),'compiler':str(compiler),'compiler_version':compiler_version,
            'compile_command':compile_command,'command':command,'graph_input_sha256':verify._sha256(input_path),
            'initial_preferences_sha256':verify._sha256(init),'goal_satisfied':False,'success_certificate':False}
    dump(args.out/'config.json',config);print(json.dumps({'stage':'tabu_started',**config}),flush=True)
    start=time.monotonic()
    with (args.out/'stdout.log').open('wb') as stdout,(args.out/'stderr.log').open('wb') as stderr:
        try:
            result=subprocess.run(command,stdout=stdout,stderr=stderr,timeout=args.seconds+20)
            code=result.returncode;timed_out=False
        except subprocess.TimeoutExpired:
            code=None;timed_out=True
    report=dict(config,native_exit_code=code,outer_timeout=timed_out,elapsed_seconds=time.monotonic()-start,
                stdout_sha256=verify._sha256(args.out/'stdout.log'),stderr_sha256=verify._sha256(args.out/'stderr.log'))
    if args.graph.read_bytes()!=raw:raise ValueError('Exact graph changed during local search')
    if best.exists():
        colors=[int(c)for c in best.read_text().split()]
        if len(colors)!=len(graph['vertices']) or any(not 0<=c<6 for c in colors):
            raise ValueError('Malformed candidate coloring')
        bad=[i for i,(u,v) in enumerate(graph['edges']) if colors[u]==colors[v]]
        candidate={'k':6,'colors':colors,'graph_sha256':graph_hash,'proper_coloring':not bad,
                   'conflicting_edge_count':len(bad),'scope':'Search state only unless independently verified proper.'}
        dump(args.out/'best-coloring.json',candidate)
        dump(args.out/'conflicting-edge-indices.json',bad)
        report.update(best_conflicts=len(bad),best_coloring_sha256=verify._sha256(args.out/'best-coloring.json'))
        if stats.exists():
            native=json.loads(stats.read_text());report['native']=native
            if native['best_conflicts']!=len(bad):raise ArithmeticError('Native and exact edge conflict counts disagree')
        if not bad:
            geometry=verify.verify_geometry(graph);verify.verify_coloring(graph,candidate,6)
            dump(args.out/'coloring.json',candidate)
            dump(args.out/'verification.json',{'graph_sha256':graph_hash,'geometry':geometry.__dict__,
                  'proper_six_coloring_verified':True,'coloring_sha256':verify._sha256(args.out/'coloring.json'),
                  'goal_satisfied':False})
            report['status']='SIX_COLORING_INDEPENDENTLY_VERIFIED'
        else:
            report['status']='UNKNOWN_IMPROPER_BEST'
    else:
        report['status']='UNKNOWN_NO_MODEL'
    dump(args.out/'result.json',report);print(json.dumps(report),flush=True)


if __name__=='__main__':main()
