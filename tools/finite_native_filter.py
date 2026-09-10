"""Native coloring search on an abstract finite quotient, for upper bounds.

This module cannot emit a Euclidean success certificate. A positive model
may later be pulled back through a separately verified edge-preserving map.
An UNSAT outcome only describes this finite graph (and is retained as
unverified here). Triangle color normalization is explicit and reversible
by a global permutation of color names.
"""
import argparse
import hashlib
import json
from pathlib import Path
import signal
import subprocess
import time

from hn.solve import clauses
from hn.solve_external import DEFAULT_SOLVER, parse_model, verify_full_cnf_model


def run(graph_path, out, k=6, timeout=180):
    graph_path = Path(graph_path).resolve()
    raw = graph_path.read_bytes(); graph = json.loads(raw)
    if graph.get('schema') != 'hn-finite-cayley-diagnostic-v1':
        raise ValueError('This filter only accepts the explicitly abstract finite schema')
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    (out/'program.py').write_bytes(Path(__file__).read_bytes())
    edge_set = {tuple(sorted(edge)) for edge in graph['edges']}
    neighbor = next(v for u, v in graph['edges'] if u == 0)
    third = next(v for v in range(len(graph['vertices']))
                 if (0, v) in edge_set and tuple(sorted((neighbor, v))) in edge_set)
    triangle = [0, neighbor, third]
    formula = list(clauses(graph, k)) + [[k*v+c+1] for c, v in enumerate(triangle)]
    cnf = out/'coloring.cnf'
    cnf.write_text(f'p cnf {k*len(graph["vertices"])} {len(formula)}\n'+
                   ''.join(' '.join(map(str, clause))+' 0\n' for clause in formula))
    proof = out/'proof.drat'; solution = out/'solution.out'
    command = [str(DEFAULT_SOLVER), '--no-colors', '--sat', '-w', str(solution), str(cnf), str(proof)]
    result = {'schema': 'hn-finite-native-filter-v1', 'scope': 'abstract finite quotient only',
              'graph_sha256': hashlib.sha256(raw).hexdigest(), 'graph': str(graph_path),
              'k': k, 'vertices': len(graph['vertices']), 'edges': len(graph['edges']),
              'triangle_normalization': triangle, 'command': command,
              'goal_satisfied': False, 'success_certificate': False,
              'solver_sha256': hashlib.sha256(DEFAULT_SOLVER.read_bytes()).hexdigest(),
              'cnf_sha256': hashlib.sha256(cnf.read_bytes()).hexdigest()}
    (out/'configuration.json').write_text(json.dumps(result, indent=2)+'\n')
    started = time.monotonic(); timed_out = False
    with (out/'stdout.log').open('wb') as stdout, (out/'stderr.log').open('wb') as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True; process.send_signal(signal.SIGTERM)
            try: returncode = process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill(); returncode = process.wait()
    result.update(elapsed_seconds=time.monotonic()-started, native_exit_code=returncode)
    if timed_out:
        result['status'] = 'UNKNOWN_TIMEOUT'
    elif returncode == 10:
        model = parse_model(solution, k*len(graph['vertices']))
        result['cnf_clauses_checked'] = verify_full_cnf_model(cnf, model)
        colors = [next(c for c in range(k) if model[k*v+c+1]) for v in range(len(graph['vertices']))]
        if any(colors[u] == colors[v] for u, v in graph['edges']):
            raise ValueError('Invalid quotient coloring')
        result['status'] = 'FINITE_COLORING_VERIFIED'
        (out/'coloring.json').write_text(json.dumps({'k': k, 'colors': colors,
             'graph_sha256': result['graph_sha256']})+'\n')
    else:
        result['status'] = 'FINITE_UNSAT_UNVERIFIED' if returncode == 20 else 'NATIVE_ERROR'
    if proof.exists():
        digest = hashlib.sha256()
        with proof.open('rb') as stream:
            for chunk in iter(lambda: stream.read(1024*1024), b''): digest.update(chunk)
        result.update(proof_sha256=digest.hexdigest(), proof_bytes=proof.stat().st_size)
    (out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--k', type=int, choices=[5, 6], default=6)
    parser.add_argument('--timeout', type=float, default=180)
    args = parser.parse_args()
    run(args.graph, args.out, args.k, args.timeout)
