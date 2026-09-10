# Hadwiger–Nelson exact construction search

**Search active; the requested result has not been obtained.** There is no
finite graph here with a verified proof of non-6-colorability. This repository
records ongoing research, exact constructions, independently checked
colorings, and unresolved solver runs. The current continuation state is in
[RESEARCH_STATE.md](RESEARCH_STATE.md).

The required certificate has two parts: distinct exact points of the
Euclidean plane with every listed edge of length exactly one, and an
independently checked refutation of the unrestricted 6-coloring CNF for
that same graph. The success gate accepts only that conjunction.

## First construction batch

| Construction | Vertices | Exact unit edges | Six-color result |
|---|---:|---:|---|
| Published Haugland G3, reconstructed | 2,131 | 12,530 | Proper coloring verified |
| G3 translated over a unit hexagon and its center | 11,848 | 89,898 | Proper coloring verified |
| Union of seven spindle-angle rotations of G3 | 14,911 | 90,404 | Proper coloring verified |
| Both equilateral completions of every G3 edge | 13,861 | 96,917 | Proper coloring verified |
| Radius-two ball from two rotated direction sets | 14,113 | 63,294 | Proper coloring verified |
| Radius-two ball from three rotated direction sets | 31,753 | 137,340 | Proper coloring verified |
| Spindle assembly of the full exact T6 path set | 32,241 | 229,108 | Proper coloring verified |

The baseline and the two mixed-direction balls also have independently
checked 5-colorings. Other secondary 5-color runs use a 60-second interruption
budget; their saved result files distinguish timeouts from SAT. A baseline
4-color run timed out. Neither a timeout nor literature's known baseline
lower bound is counted as a fresh UNSAT certificate.

See `runs/summary.json` for the complete machine-readable outcomes,
`runs/**/result.json` for solver statistics, and the accompanying coloring,
CNF, and verification files for evidence. Every target candidate in the first
batch is rejected by a coloring witness. Later batches and their algebraic
obstruction checks are recorded in the continuation state and research notes.

## Exact geometry and independent verification

Points are complex numbers `a(t) + sqrt(5)*b(t)`, with
`t = exp(pi*i/21)`, and rational coefficient polynomials of degree below 12.
The fixed cyclotomic polynomial is

```
Phi_42(t) = t^12 + t^11 - t^9 - t^8 + t^6 - t^4 - t^3 + t + 1.
```

The selected real embedding takes `sqrt(5)>0`; plane coordinates are
`(Re(z), Im(z))`. Equality and norms are reduced exactly in this number
field. Complex conjugation sends `t` to `t^-1` and fixes `sqrt(5)`.
An edge is valid precisely when `(z-w)*conjugate(z-w)=1`.

`hn/field.py` uses integer coefficient vectors and common denominators for
construction. `tools/verify.py` imports none of the construction code and
uses FLINT rational polynomials to recheck all vertices and all edges.
It also independently regenerates the CNF and checks coloring assignments.
The quadratic extension is nontrivial because prime 5 ramifies in
`Q(sqrt(5))` but not in `Q(zeta_42)`, which makes the coefficient-pair
representation unique.

Candidate-edge enumeration sometimes uses floating-point filtering to
reduce work. **Every included edge is then checked exactly**, and is checked
again independently. No induced-graph completeness claim is made. Missing
a unit pair weakens a graph; it cannot create a false lower-bound witness.

The CNF has exactly one color per vertex and unequal colors on every edge.
There are no virtual edges, fixed colors, or symmetry-breaking assumptions.
`hn/solve.py` can emit DRAT proofs, but an UNSAT status alone never sets a
success flag. The verifier requires a proof for the regenerated **six**-color
formula and invokes the separately built DRAT-trim checker. Its executable
and proof inputs are hashed. `runs/verifier-proof-selftest` contains a real
non-2-colorability test and a forged-proof rejection; neither is a target
certificate.

## Reproduction

From this repository:

```sh
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.lock
git clone https://github.com/marijnheule/drat-trim.git sources/drat-trim
git -C sources/drat-trim checkout 2e3b2dc0ecf938addbd779d42877b6ed69d9a985
make -C sources/drat-trim drat-trim
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python tools/verify.py --graph runs/t6-spindle/graph.json --k 6 --coloring runs/t6-spindle/six-color/coloring.json
```

Reconstruct or run a fresh candidate into a new directory:

```sh
.venv/bin/python -m hn.construct g3 --output artifacts/rebuilt-g3.json
.venv/bin/python -m hn.ball_search --layers 3 --out runs/new-mixed-ball3 --timeout 60
.venv/bin/python -m hn.path_search --out runs/new-t6-spindle --timeout 60
```

A future claimed target must pass this command with an actual proof:

```sh
.venv/bin/python tools/verify.py --graph candidate.json --k 6 --cnf coloring.cnf --proof proof.drat --checker sources/drat-trim/drat-trim
```

## Literature and provenance

The [research notes](research/literature.md) link primary papers and describe
which coordinate fields already have five-colorability obstructions.
The inspected literature supports bounds 5 through 7, including
[Haugland's August 2026 revision](https://arxiv.org/html/2608.04542v4).
The actual Appendix A input retained for reconstruction is explicitly
pinned to [v1](https://arxiv.org/html/2608.04542v1), licensed CC BY 4.0 and
attributed to Jan Kristian Haugland. Its SHA-256 appears in baseline graph
metadata. Source and dependency details are in `sources/README.md` and
`requirements.lock`; retained files are inventoried in `MANIFEST.sha256`.

Work was performed in `/Users/rayanpal/hadwiger-nelson`, initialized as a
local Git repository. Both available GitHub access paths reported the
supplied remote `rayanpal/hadwiger-nelson` as not found.
