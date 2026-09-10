"""Complete a pole-5 candidate with Heule's pinned exact 510-point realization.

The source embedding is checked edge by edge. The source's published
chromatic claim is not treated as an independently verified lower bound.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from .heule_import import load_vtx
from .multiquadratic import Point
from .mq_construct import canonical_points, export, graph_from_points
from .mq_intersections import core_indices
from .solve_external import solve_external


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph', type=Path, default=Path('runs/mq-pole5-adaptive-v1/graph.json'))
    parser.add_argument('--source', type=Path, default=Path('sources/heule-cnp-sat/510.vtx'))
    parser.add_argument('--out', type=Path, default=Path('runs/mq-pole5-heule510-v1'))
    parser.add_argument('--timeout', type=float, default=300)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'program.py').write_bytes(Path(__file__).read_bytes())
    original = json.loads(args.graph.read_text())
    points = [Point.from_json(v) for v in original['vertices']]
    old_set = set(points)
    source_points = load_vtx(args.source)
    source_graph = graph_from_points(source_points)
    source_canonical = canonical_points(source_points)
    missing = [p for p in source_canonical if p not in old_set]
    if len(missing) > 471:
        raise ValueError('The completion would exceed the remaining bounded addition budget')
    config = {'source_graph':str(args.graph), 'source_graph_sha256':digest(args.graph),
              'source_coordinates':str(args.source), 'source_coordinate_sha256':digest(args.source),
              'source_repository_commit':'bb414955a6ef5f49f7df2b245b1e778aa67c068a',
              'source_vertices':len(source_points), 'source_edges':len(source_graph['edges']),
              'new_vertices':len(missing), 'total_additions_since_original_pole_graph':529+len(missing),
              'program_sha256':digest(__file__), 'native_six_timeout_seconds':args.timeout,
              'goal_satisfied':False,
              'scope':'Exact isometric inclusion of the published source realization; no source chromatic assertion is accepted as a fresh proof.'}
    dump(args.out/'config.json', config)
    print(json.dumps({'stage':'completion_started', **config}), flush=True)
    graph = export(points+missing, args.out/'graph.json',
                   dict(config, family='pole5-five-model-extension-plus-heule510'))
    expanded = canonical_points(points+missing)
    index = {p:i for i,p in enumerate(expanded)}
    old_to_new = [index[p] for p in points]
    source_to_full = [index[p] for p in source_canonical]
    new_ids = {index[p] for p in missing}
    edges = {tuple(e) for e in graph['edges']}
    for u,v in original['edges']:
        a,b = sorted((old_to_new[u], old_to_new[v]))
        if (a,b) not in edges:
            raise ArithmeticError('Completion export omitted an original edge')
    for u,v in source_graph['edges']:
        a,b = sorted((source_to_full[u], source_to_full[v]))
        if (a,b) not in edges:
            raise ArithmeticError('Completion export omitted a source edge')
    dump(args.out/'embedding.json', {'full_graph_sha256':digest(args.out/'graph.json'),
                                    'old_to_full_vertex_index':old_to_new,
                                    'source_canonical_coordinates':[p.to_json() for p in source_canonical],
                                    'source_to_full_vertex_index':source_to_full,
                                    'source_edges':source_graph['edges'],
                                    'new_full_vertex_indices':sorted(new_ids),
                                    'all_original_and_source_edges_present':True,
                                    'fresh_source_non4_certificate_obtained':False})
    core = core_indices(graph, 6)
    kinds = Counter('old-old' if u not in new_ids and v not in new_ids else
                    'new-new' if u in new_ids and v in new_ids else 'old-new'
                    for u,v in graph['edges'])
    summary = {'vertices':len(expanded), 'edges':len(edges), 'new_vertices':len(missing),
               'edge_classes':dict(kinds), 'six_core_vertices':len(core),
               'six_core_edges':sum(u in core and v in core for u,v in edges),
               'new_vertices_in_six_core':len(new_ids&core),
               'source_embedding_edges_checked':len(source_graph['edges']),
               'goal_satisfied':False}
    dump(args.out/'summary.json', summary)
    print(json.dumps({'stage':'completion_exported', **summary}), flush=True)
    result = solve_external(args.out/'graph.json', args.out/'native-six-color',
                            timeout=args.timeout, configuration='sat')
    summary.update(native_six_status=result['status'], goal_satisfied=result.get('goal_satisfied',False))
    dump(args.out/'summary.json', summary)


if __name__ == '__main__':
    main()
