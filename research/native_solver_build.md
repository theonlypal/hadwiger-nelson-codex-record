# Native SAT portfolio and independent verification

The native runner is `hn/solve_external.py`. It uses CaDiCaL, giving a different
solver implementation and preprocessing/search configuration from the earlier
PySAT Glucose runs. A timeout remains UNKNOWN.

## Pinned build

The source is the official [CaDiCaL repository](https://github.com/arminbiere/cadical),
pinned at commit `c60730422e758ef1cebe7aeddf2dda31c996bf04`. Its reported version
is 3.0.1. The source lives in `sources/cadical`.

Commands actually used:

```sh
git clone --depth 1 https://github.com/arminbiere/cadical.git sources/cadical
git -C sources/cadical checkout --detach c60730422e758ef1cebe7aeddf2dda31c996bf04
cd sources/cadical
./configure
make -j2 cadical
```

The resulting macOS arm64 executable is `sources/cadical/build/cadical`, with
SHA-256 `39143d3955b59965d9e8ec2b7df71206fa4380d6f17c3c8e3254dce653897fd0`.
The build used the host compiler selected by `configure`, with `-O3 -DNDEBUG`.
Build configuration is preserved in the source checkout's generated makefiles,
`configure.log`, and `build/build.hpp`. This binary hash is platform-specific.

The official [README at the pinned commit](https://github.com/arminbiere/cadical/blob/c60730422e758ef1cebe7aeddf2dda31c996bf04/README.md)
documents the configure/make build. The [pinned command-line implementation](https://github.com/arminbiere/cadical/blob/c60730422e758ef1cebe7aeddf2dda31c996bf04/src/cadical.cpp)
documents the solution-file option, predefined SAT/UNSAT configurations, and
binary DRAT output. These options were also checked against the built binary's
`--help` output.

## Runner contract

Example from the repository root:

```sh
./.venv/bin/python -m hn.solve_external artifacts/heule_sqrt2_star_mq.json \
  --out runs/native-cadical-sqrt2-star-k6 --timeout 300 --configuration sat
```

Every invocation creates a new output directory. It validates the exact
geometry independently, regenerates the unrestricted canonical 6-coloring
CNF, and stores the native command, source revision, hashes, and environment.
No fixed colors or symmetry clauses are added.

The directory retains `coloring.cnf`, `solver.stdout.log`, `solver.stderr.log`,
`solution.out`, and the streamed binary `proof.drat`, even when incomplete.
The raw solver outputs are not rewritten. A hard timeout sends SIGTERM, then
SIGKILL only if the process fails to exit within five seconds.

After a completed SAT/UNKNOWN run, the proof is losslessly compressed to
`proof.drat.gz`. The runner hashes the decompressed bytes and rechecks the
original before removing the uncompressed file. `result.json` records both
hashes, both byte counts, the compressed path, and the successful decompression
check. Active proofs and UNSAT proofs undergoing verification are not touched.

Exit 10 is accepted as SAT only after checking a complete assignment to every
original CNF variable, every original clause, exactly one color per vertex,
and every graph edge. The proper coloring and verification report are saved.

Exit 20 triggers the field-specific independent certificate checker. Success
requires exact geometry, exact regeneration of the proof's unrestricted
6-coloring CNF, and the external DRAT check. The source solver's status alone
does not satisfy the goal. A timeout preserves its partial proof and returns
UNKNOWN without invoking an UNSAT success gate.

`--k 5` is an optional diagnostic; the default remains six. A five-color UNSAT
result is checked against the regenerated five-color CNF and exact geometry,
then labeled `NON5_VERIFIED` with `goal_satisfied: false`. It cannot pass either
public non-six success gate. The planned expensive MQ five-color diagnostic
was cancelled after explicit algebraic five-colorings eliminated that ring.

## Pipeline checks

`runs/native-cadical-selftest-sat` contains a native SAT run on an exact triangle
at six colors, including the checked complete model.

`runs/native-cadical-selftest-unsat` contains a native UNSAT run on the triangle
at **two** colors. The emitted binary DRAT proof was independently accepted by
the existing `drat-trim` build. Five malformed model files were rejected by
the full-model parser. These are pipeline checks and do not satisfy the
requested non-6-colorability condition.
