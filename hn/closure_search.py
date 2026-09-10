"""Complete both equilateral triangles on each existing unit edge, exactly."""
import argparse
import json
from pathlib import Path

from .field import Element, ZETA
from .search import export
from .solve import solve


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--timeout',type=float,default=60)
    a=p.parse_args()
    out=Path(a.out)
    out.mkdir(parents=True,exist_ok=False)
    graph=json.loads(Path(a.seed).read_text())
    if graph['field']['cyclotomic_order']!=42 or graph['field']['radicand']!=5:
        raise ValueError('Unexpected seed field')
    old=[Element.from_json(v) for v in graph['vertices']]
    points=set(old)
    plus,minus=ZETA**7,ZETA**(-7)
    for u,v in graph['edges']:
        d=old[v]-old[u]
        if not d.is_unit():
            raise ValueError('Nonunit seed edge')
        points.add(old[u]+d*plus)
        points.add(old[u]+d*minus)
    export(points,out/'graph.json',{'family':'equilateral-completion','seed':str(Path(a.seed).resolve())})
    result=solve(out/'graph.json',6,out/'six-color',a.timeout)
    if result['status']=='SAT':
        solve(out/'graph.json',5,out/'five-color',a.timeout)


if __name__=='__main__':
    main()
