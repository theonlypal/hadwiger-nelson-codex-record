"""Check explicit additive five-colorings before a target SAT search.

The only pruning result is an independently verified proper coloring of
the entire exported exact graph. Failure of this ansatz says nothing about
its chromatic number. No result from this module satisfies the target.
"""
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path


def check_functional(graph_path, functional, out):
    from tools import verify_tower
    path = Path(graph_path)
    raw = path.read_bytes()
    graph = json.loads(raw)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    report = {'graph_sha256': hashlib.sha256(raw).hexdigest(),
              'goal_satisfied': False, 'success_certificate': False,
              'status': 'NOT_APPLICABLE'}
    if (graph.get('schema') != 'hn-cyclotomic-tower-v1'
            or graph.get('field', {}).get('cyclotomic_order') != 42
            or graph.get('field', {}).get('radicands') != [5, 33]):
        report['reason'] = 'Unsupported coefficient basis'
    elif (len(functional) != 48
          or any(type(c) is not int or not 0 <= c < 5 for c in functional)):
        raise ValueError('Expected 48 coefficients in {0,1,2,3,4}')
    else:
        @lru_cache(None)
        def residue(value):
            value = Fraction(value)
            if value.denominator % 5 == 0:
                raise ZeroDivisionError('A coordinate has a pole in this coefficient reduction')
            return value.numerator * pow(value.denominator, -1, 5) % 5
        try:
            colors = [sum(a*residue(b) for a, b in zip(functional,
                      (c for part in vertex['parts'] for c in part))) % 5
                      for vertex in graph['vertices']]
        except ZeroDivisionError as exc:
            report['reason'] = str(exc)
        else:
            bad = [edge for edge in graph['edges'] if colors[edge[0]] == colors[edge[1]]]
            report.update({'linear_functional': functional, 'monochromatic_edge_count': len(bad),
                           'first_monochromatic_edges': bad[:32]})
            if bad:
                report['status'] = 'FUNCTIONAL_REJECTED'
            else:
                coloring = {'k': 5, 'colors': colors, 'graph_sha256': report['graph_sha256']}
                # These independent checks validate every edge and every exact
                # coordinate, rather than trusting the candidate direction set.
                verify_tower.verify_coloring(graph, coloring, 5)
                geometry = verify_tower.verify_geometry(graph)
                target = out/'coloring.json'
                target.write_text(json.dumps(coloring, separators=(',', ':'))+'\n')
                report.update({'status': 'FIVE_COLORING_VERIFIED',
                               'geometry': geometry.__dict__,
                               'coloring_sha256': hashlib.sha256(target.read_bytes()).hexdigest()})
    (out/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def search_functional(graph_path, out, timeout=60):
    """An optional bounded ansatz search; only a checked model prunes."""
    from tools.direction_color_screen import run
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    try:
        candidate = run(graph_path, out/'search', timeout)
    except (ValueError, AssertionError) as exc:
        report = {'status': 'NO_LINEAR_CERTIFICATE', 'reason': str(exc),
                  'goal_satisfied': False, 'success_certificate': False}
    else:
        if candidate['status'] == 'SAT_VERIFIED':
            report = check_functional(graph_path, candidate['linear_functional'], out/'verification')
        else:
            report = {'status': 'NO_LINEAR_CERTIFICATE', 'ansatz_status': candidate['status'],
                      'goal_satisfied': False, 'success_certificate': False}
    (out/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    return report
