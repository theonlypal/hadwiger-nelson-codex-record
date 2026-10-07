# Hadwiger-Nelson Codex record

Rayan Pal's local research record from September 10, 2026. Published October 6, 2026 (America/Los_Angeles).

The primary target was a unit-distance graph requiring seven colors. The secondary target was a graph requiring at least six colors, establishing `χ(R²) ≥ 6`.

At `2026-09-10T18:00:24.737Z`, the recorded user instruction to Codex said:

> Secondary success condition: construct a finite exact unit-distance graph that is not 5-colorable.

## Evidence

- [Timestamped original records](transcript/records.jsonl): eight selected records, unchanged.
- [Readable transcript](transcript/readable.txt): initial target, progress reports, execution output, and explicit secondary target.
- [Source line numbers, hashes, and omissions](transcript/extraction.json).
- [Codex Desktop session metadata](transcript/session-metadata.json).
- [Original commit IDs and dates](provenance/commits.json).
- [Browse the original code and Git history](https://github.com/theonlypal/hadwiger-nelson-codex-record/tree/original-history).
- [Source file timestamps and hashes](provenance/source-files.json).

## Recorded outcome

One constructed graph has **125,929 vertices and 991,423 edges**.

| Test | Recorded result |
| --- | --- |
| Six-color search | SAT in 101.986 seconds; valid six-coloring checked |
| Five-color search | UNKNOWN_TIMEOUT after 300.242 seconds |

[Six-color result](runs/mixed-g3-v31-odd-rotation-v1/six-color/result.json) · [Verification](runs/mixed-g3-v31-odd-rotation-v1/six-color/verification.json) · [Five-color result](runs/mixed-g3-v31-odd-rotation-v1/five-color-phases-v1/native/result.json)

No new lower bound was proved by these runs. The five-color trace is incomplete. The six-color trace came from a SAT run. Neither is a completed UNSAT certificate.

## Downloads and checks

[Release v1.0.0](https://github.com/theonlypal/hadwiger-nelson-codex-record/releases/tag/v1.0.0) contains:

- `original-history.bundle`: original Git history, preserving all three commit IDs.
- `hadwiger-nelson-evidence.tar.gz`: graph, CNFs, coloring, logs, code at `966b629`, and solver traces.
- `SHA256SUMS.txt`: release-asset checksums.

Run `python3 verify_evidence.py` to check the repository records. After extracting the evidence archive, run `python3 verify_evidence.py --archive-root /path/to/extracted/evidence`. These are integrity checks, not a new solver run or mathematical proof verification.

September 10 dates are local transcript timestamps, filesystem timestamps, and Git metadata. This publication does not independently attest to those historical dates. Original files were left untouched. Excerpts omit unrelated runtime context; retained source records were not rewritten.
