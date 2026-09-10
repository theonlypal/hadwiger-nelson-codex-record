"""Inventory completed evidence without promoting solver output to a proof."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def main():
    rows=[]
    for path in sorted((ROOT/'runs').rglob('result.json')):
        result=json.loads(path.read_text())
        rows.append(dict(path=str(path.relative_to(ROOT)),**result))
    target_rows=[row for row in rows if row.get('k')==6]
    all_rejected=bool(target_rows) and all(row.get('status')=='SAT' for row in target_rows)
    summary={
        'target':'finite exact Euclidean unit-distance graph not 6-colorable',
        'goal_satisfied':False,
        'certificate_present':False,
        'target_solver_runs':sum(row.get('k')==6 for row in rows),
        'target_SAT_runs':sum(row.get('k')==6 and row.get('status')=='SAT' for row in rows),
        'runs':rows,
        'interpretation':('All recorded target candidates have proper 6-colorings; no non-6-colorability certificate has been obtained.'
                          if all_rejected else 'This inventory does not establish non-6-colorability; inspect each run and independently verify any claimed proof.')
    }
    (ROOT/'runs'/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    paths=[]
    for top in ['hn','tools','tests','research','artifacts','runs']:
        paths.extend(p for p in (ROOT/top).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    paths.extend(p for p in (ROOT/'sources').glob('*') if p.is_file())
    paths.extend(ROOT/p for p in ['README.md','requirements.lock','.gitignore'])
    lines=[]
    for path in sorted(paths):
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f'{digest}  {path.relative_to(ROOT)}\n')
    (ROOT/'MANIFEST.sha256').write_text(''.join(lines))
    print(json.dumps({k:v for k,v in summary.items() if k!='runs'}))


if __name__=='__main__':
    main()
