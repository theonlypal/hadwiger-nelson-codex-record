"""Generate checked seven-color starts, then search for a proper six-coloring.

Floating coordinates choose candidates only. Every positive coloring is
checked against every exported edge and the independent exact verifier.
An improper state or exhausted search establishes no chromatic lower bound.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from hn.solve_external import common, verify_mq, verify_tower, verify_voronov


def load_points(graph):
    schema = graph['schema']
    if schema == 'hn-exact-v1':
        from hn.field import Element as Point
        verifier = common
    elif schema == 'hn-multiquadratic-v1':
        from hn.multiquadratic import Point
        verifier = verify_mq
    elif schema == 'hn-cyclotomic-tower-v1':
        from hn.cyclotomic_tower import Element as Point
        verifier = verify_tower
    elif schema == 'hn-voronov-field-v1':
        from hn.voronov_field import Point
        verifier = verify_voronov
    else:
        raise ValueError('Unsupported exact field')
    z = [Point.from_json(v).approximate() for v in graph['vertices']]
    return np.array(z, dtype=complex), verifier


def hex_colors(z, angle, side, shift):
    z = z * np.exp(1j * angle) + shift
    # Pointy regular hexagons: Cartesian to fractional axial coordinates.
    q = z.real / (np.sqrt(3) * side) - z.imag / (3 * side)
    r = 2 * z.imag / (3 * side)
    cube = np.stack((q, -q-r, r), axis=1)
    rounded = np.rint(cube).astype(np.int64)
    largest = np.argmax(np.abs(rounded-cube), axis=1)
    rows = np.arange(len(z))
    rounded[rows, largest] -= rounded.sum(axis=1)
    return (rounded[:, 0] + 3 * rounded[:, 2]) % 7


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--trials', type=int, default=96)
    parser.add_argument('--seconds', type=float, default=120)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    raw = args.graph.read_bytes(); graph = json.loads(raw)
    digest = hashlib.sha256(raw).hexdigest()
    z, verifier = load_points(graph)
    edges = np.array(graph['edges'], dtype=np.int64)
    u, v = edges.T
    rng = np.random.default_rng(args.seed)
    best = None; best_count = len(edges)+1; history = []
    for trial in range(args.trials):
        angle = float(rng.uniform(0, np.pi/3))
        side = float(rng.uniform(.4, .499))
        shift = complex(*rng.uniform(-1, 1, size=2))
        colors7 = hex_colors(z, angle, side, shift)
        if np.any(colors7[u] == colors7[v]):
            raise ArithmeticError('Hex candidate failed actual graph edge check')
        counts = np.bincount(np.concatenate((u*7+colors7[v], v*7+colors7[u])),
                             minlength=len(z)*7).reshape(len(z), 7)
        for removed in range(7):
            keep = [c for c in range(7) if c != removed]
            mask = colors7 == removed
            choices = counts[mask][:, keep]
            residual = int(choices.min(axis=1).sum())
            if residual < best_count:
                colors = colors7.copy()
                colors[mask] = np.array(keep)[choices.argmin(axis=1)]
                colors -= colors > removed
                if int(np.sum(colors[u] == colors[v])) != residual:
                    raise ArithmeticError('Color removal conflict accounting failed')
                best = colors; best_count = residual
                chosen = dict(trial=trial, angle=angle, side=side,
                              shift=[shift.real, shift.imag], removed=removed)
                seven = colors7.tolist()
        history.append(dict(trial=trial, best_conflicts=best_count))
    (args.out/'hex-search.json').write_text(json.dumps(dict(history=history, selected=chosen,
        best_conflicts=best_count, graph_sha256=digest, seed=args.seed))+'\n')
    witness7 = dict(k=7, colors=seven, graph_sha256=digest)
    verifier.verify_coloring(graph, witness7, 7)
    (args.out/'seven-coloring.json').write_text(json.dumps(witness7)+'\n')
    (args.out/'program.py').write_bytes(Path(__file__).read_bytes())
    source = Path('tools/tabu_color.cpp')
    (args.out/'tabu_color.cpp').write_bytes(source.read_bytes())
    binary = (args.out/'tabu-color').resolve()
    subprocess.run(['c++', '-O3', '-std=c++17', str(args.out/'tabu_color.cpp'), '-o', str(binary)], check=True)
    inp = args.out/'graph.edges'; init = args.out/'initial.txt'
    inp.write_text(f'{len(z)} {len(edges)} 6\n'+''.join(f'{a} {b}\n' for a,b in edges))
    init.write_text('\n'.join(map(str, best))+'\n')
    model = args.out/'best.txt'; stats = args.out/'native-result.json'
    cmd = [str(binary), str(inp), str(init), str(model), str(stats), str(args.seed),
           str(args.seconds), '1000000000', '64']
    print(json.dumps(dict(stage='hex_repair_started', best_conflicts=best_count,
                         selected=chosen, command=cmd)), flush=True)
    started = time.monotonic()
    with (args.out/'stdout.log').open('wb') as out, (args.out/'stderr.log').open('wb') as err:
        result = subprocess.run(cmd, stdout=out, stderr=err, timeout=args.seconds+20)
    colors = np.array([int(c) for c in model.read_text().split()])
    if len(colors) != len(z) or np.any((colors < 0) | (colors >= 6)):
        raise ValueError('Malformed local model')
    bad = np.flatnonzero(colors[u] == colors[v]).tolist()
    report = dict(graph_sha256=digest, elapsed_seconds=time.monotonic()-started,
                  native_exit_code=result.returncode, best_conflicts=len(bad),
                  native=json.loads(stats.read_text()), goal_satisfied=False,
                  success_certificate=False, status='UNKNOWN_IMPROPER_BEST')
    if args.graph.read_bytes() != raw:
        raise ValueError('Graph changed during search')
    if report['native']['best_conflicts'] != len(bad):
        raise ArithmeticError('Native model conflict count mismatch')
    witness = dict(k=6, colors=colors.tolist(), graph_sha256=digest,
                   proper_coloring=not bad, conflicting_edge_count=len(bad))
    (args.out/'best-coloring.json').write_text(json.dumps(witness)+'\n')
    (args.out/'conflicting-edge-indices.json').write_text(json.dumps(bad)+'\n')
    if not bad:
        report['geometry'] = verifier.verify_geometry(graph).__dict__
        verifier.verify_coloring(graph, witness, 6)
        report['status'] = 'SIX_COLORING_INDEPENDENTLY_VERIFIED'
        (args.out/'coloring.json').write_text(json.dumps(witness)+'\n')
    (args.out/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
