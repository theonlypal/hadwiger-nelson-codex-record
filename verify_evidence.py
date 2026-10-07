#!/usr/bin/env python3
"""Offline integrity and recorded-outcome checks. Does not run a solver."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda: f.read(1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def check_manifest(root):
    count = 0
    for line in (root / 'SHA256SUMS.txt').read_text().splitlines():
        expected, name = line.split('  ', 1)
        p = (root / name).resolve()
        if not p.is_relative_to(root.resolve()) or digest(p) != expected:
            raise ValueError('Hash mismatch or unsafe path: ' + name)
        count += 1
    return count


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive-root', type=Path)
    args = ap.parse_args()
    root = Path(__file__).resolve().parent
    n = check_manifest(root)
    index = json.loads((root / 'transcript/extraction.json').read_text())
    lines = (root / 'transcript/records.jsonl').read_bytes().splitlines(keepends=True)
    if len(lines) != len(index['records']):
        raise ValueError('Transcript length mismatch')
    for raw, item in zip(lines, index['records']):
        if hashlib.sha256(raw).hexdigest() != item['sha256_original_line_including_newline']:
            raise ValueError('Transcript record hash mismatch')
        if json.loads(raw)['timestamp'] != item['timestamp']:
            raise ValueError('Transcript timestamp mismatch')
    run = root / 'runs/mixed-g3-v31-odd-rotation-v1'
    six = json.loads((run / 'six-color/result.json').read_text())
    five = json.loads((run / 'five-color-phases-v1/native/result.json').read_text())
    if (six['status'], six['k'], six['proper_coloring_checked']) != ('SAT', 6, True):
        raise ValueError('Unexpected six-color record')
    if (five['status'], five['k'], five['partial_proof_retained']) != ('UNKNOWN_TIMEOUT', 5, True):
        raise ValueError('Unexpected five-color record')
    for record in (six, five):
        if record['goal_satisfied'] or record['success_certificate']:
            raise ValueError('Unexpected success claim')
        if (record['vertices'], record['edges']) != (125929, 991423):
            raise ValueError('Graph size mismatch')
    print(f'PASS: {n} repository file hashes; {len(lines)} transcript records; recorded outcomes')
    if args.archive_root:
        archive = args.archive_root.resolve()
        n = check_manifest(archive)
        for item in json.loads((root/'provenance/source-files.json').read_text()):
            p = archive / item['path']
            if digest(p) != item['sha256'] or p.stat().st_size != item['bytes']:
                raise ValueError('Source file mismatch: '+item['path'])
        print(f'PASS: {n} archive file hashes; archived run files match source manifest')


if __name__ == '__main__':
    main()
