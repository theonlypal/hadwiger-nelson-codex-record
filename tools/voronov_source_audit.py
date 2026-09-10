#!/usr/bin/env python3
"""Safely parse the author's base coordinate lists and compare exact points.

Only integer constants, rational arithmetic, and square roots of positive
rationals with squarefree part dividing6 are accepted. This is a data parser;
it neither evaluates Python source nor executes the downloaded notebook.
"""
from __future__ import annotations

import ast
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from hn.multiquadratic import MQ,Point,ONE,SQRT2,SQRT3
from hn.voronov_base import m73_points,m2_points


def root_rational(value):
    if any(value.coefficients[1:]) or value.coefficients[0]<=0:
        raise ValueError('Expected a positive rational radicand')
    n=value.coefficients[0];d=value.denominator;rest=n*d;squarefree=1
    for prime in (2,3):
        exponent=0
        while rest%prime==0:rest//=prime;exponent+=1
        if exponent%2:squarefree*=prime
    # Remaining factors must be a rational square.
    scale=isqrt(n*d//squarefree)
    if scale*scale*squarefree!=n*d:raise ValueError('Unexpected squarefree factor')
    return (ONE*scale/d)*{1:ONE,2:SQRT2,3:SQRT3,6:SQRT2*SQRT3}[squarefree]


def evaluate(node):
    if isinstance(node,ast.Constant) and type(node.value) is int:return ONE*node.value
    if isinstance(node,ast.UnaryOp):
        value=evaluate(node.operand)
        if isinstance(node.op,ast.USub):return -value
        if isinstance(node.op,ast.UAdd):return value
    if isinstance(node,ast.BinOp):
        a,b=evaluate(node.left),evaluate(node.right)
        if isinstance(node.op,ast.Add):return a+b
        if isinstance(node.op,ast.Sub):return a-b
        if isinstance(node.op,ast.Mult):return a*b
        if isinstance(node.op,ast.Div):return a/b
    if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=='Sqrt':
        if len(node.args)==1 and not node.keywords:return root_rational(evaluate(node.args[0]))
    raise ValueError('Unrecognized coordinate data syntax')


@lru_cache(maxsize=None)
def coordinate(text):
    return evaluate(ast.parse(text.strip().replace('[','(').replace(']',')'),mode='eval').body)


def parse_points(path):
    points=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip():continue
        if not line.startswith('{') or not line.endswith('}'):raise ValueError('Expected one coordinate pair per line')
        pair=line[1:-1].split(',')
        if len(pair)!=2:raise ValueError('Unexpected coordinate pair')
        points.append(Point(coordinate(pair[0]),coordinate(pair[1])))
    return points


def main():
    started=time.monotonic();base=Path('sources/voronov-dist-graphs')
    manifest=json.loads((base/'manifest.json').read_text())
    hashes={row['path']:row['sha256'] for row in manifest['files']}
    report={'scope':'exact primary-coordinate provenance; no chromatic conclusion',
            'repository':manifest['repository'],'revision':manifest['revision'],
            'program_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'stages':[],
            'success_certificate':False}
    for stage,expected in ((1,set(m73_points())),(2,set(m2_points())),(3,None)):
        relative=f'plane/series 2/vtx/s2_M{stage}.vtx';path=base/relative
        assert hashlib.sha256(path.read_bytes()).hexdigest()==hashes[relative]
        parsed=parse_points(path);actual=set(parsed)
        if stage==3:
            artifact=Path('artifacts/voronov_m32257_points.json')
            expected={Point.from_json(p) for p in json.loads(artifact.read_text())['vertices']}
        assert len(actual)==len(parsed) and actual==expected
        row={'stage':stage,'primary_file':str(path),'primary_sha256':hashes[relative],
             'vertices':len(actual),'all_coordinates_match_exactly':True}
        report['stages'].append(row);print(json.dumps(row),flush=True)
    report['elapsed_seconds']=time.monotonic()-started
    Path('research/voronov_source_audit.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
