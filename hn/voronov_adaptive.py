"""Add exact centers that block a saved proper six-color assignment.

Blocking one assignment is only a proposal rule. Each augmented graph gets
a new unrestricted coloring search and an independent geometry check.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from .adaptive import candidate_pool
from .voronov_field import Point
from .solve_external import solve_external
from tools import verify_voronov as verifier


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--graph', type=Path, required=True)
    p.add_argument('--coloring', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--rounds', type=int, default=3)
    p.add_argument('--batch', type=int, default=1000)
    p.add_argument('--directions', type=int, default=216)
    p.add_argument('--max-proposals', type=int, default=3000)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'program.py').write_bytes(Path(__file__).read_bytes())
    graph_path, color_path = args.graph, args.coloring
    for iteration in range(args.rounds):
        start = time.monotonic()
        graph_raw = graph_path.read_bytes(); graph = json.loads(graph_raw)
        witness = json.loads(color_path.read_text())
        digest = hashlib.sha256(graph_raw).hexdigest()
        if witness.get('graph_sha256') != digest:
            raise ValueError('Source model hash does not bind this graph')
        verifier.verify_geometry(graph); verifier.verify_coloring(graph, witness, 6)
        points = [Point.from_json(v) for v in graph['vertices']]
        frequencies = Counter()
        for u,v in graph['edges']:
            d = points[v]-points[u]
            frequencies[d] += 1; frequencies[-d] += 1
        directions = [d for d, _ in frequencies.most_common(args.directions)]
        if not all(d.is_unit() for d in directions):
            raise ArithmeticError('Non-unit candidate direction')
        dest = args.out/f'round-{iteration:03d}'; dest.mkdir()
        print(json.dumps(dict(stage='candidate_screen', iteration=iteration,
                              vertices=len(points), directions=len(directions))), flush=True)
        proposals, screen = candidate_pool(points, witness['colors'], directions, args.max_proposals)
        chosen = [entry for entry in proposals if entry[0] == 6][:args.batch]
        (dest/'screen.json').write_text(json.dumps(dict(screen, selected=len(chosen),
            current_coloring_blocked=bool(chosen), graph_sha256=digest,
            direction_count=len(directions), goal_satisfied=False), indent=2)+'\n')
        if not chosen:
            print(json.dumps(dict(stage='no_blocking_center', **screen)), flush=True)
            break
        old_n = len(points)
        points.extend(entry[2] for entry in chosen)
        # Reuse the already verified old edges; propose and exactly accept
        # every new incident edge. Missing edges cannot help a false lower bound.
        edges = list(graph['edges']); cache = {}
        xy = np.array([[z.approximate().real, z.approximate().imag] for z in points])
        for v in range(old_n, len(points)):
            delta = xy[:v]-xy[v]
            dsq = np.einsum('ij,ij->i', delta, delta)
            for u in np.flatnonzero(np.abs(dsq-1) < 1e-8):
                d = points[v]-points[int(u)]
                if d not in cache: cache[d] = d.is_unit()
                if cache[d]: edges.append([int(u),v])
        graph = dict(graph, vertices=[point.to_json() for point in points], edges=edges,
            metadata=dict(family='voronov-current-color-blocking-centers',
                source_graph_sha256=digest, source_coloring_sha256=hashlib.sha256(color_path.read_bytes()).hexdigest(),
                iteration=iteration, added_vertices=len(chosen),
                added_edges=len(edges)-len(json.loads(graph_raw)['edges']),
                construction_seconds=time.monotonic()-start,
                induced_graph_claimed=False, success_certificate=False))
        graph_path = dest/'graph.json'
        graph_path.write_text(json.dumps(graph, separators=(',', ':'))+'\n')
        print(json.dumps(dict(stage='augmented_graph', vertices=len(points), edges=len(edges),
                              metadata=graph['metadata'])), flush=True)
        command = [sys.executable, '-m', 'tools.hex_seed_repair', '--graph', str(graph_path),
                   '--out', str(dest/'hex-six'), '--seconds', '60', '--trials', '96']
        with (dest/'hex-six.log').open('wb') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        result = json.loads((dest/'hex-six/result.json').read_text())
        if result['status'] == 'SIX_COLORING_INDEPENDENTLY_VERIFIED':
            color_path = dest/'hex-six/coloring.json'
            print(json.dumps(dict(stage='proper_six_coloring_verified', iteration=iteration)), flush=True)
            continue
        result = solve_external(graph_path, dest/'native-six', 300, configuration='sat')
        if result['status'] != 'SAT':
            print(json.dumps(dict(stage='requires_followup', result=result)), flush=True)
            break
        color_path = dest/'native-six/coloring.json'


if __name__ == '__main__':
    main()
