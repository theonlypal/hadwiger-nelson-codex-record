#!/usr/bin/env python3
"""Independent finite-certificate checker for hn-exact-v1 graphs.

No search/construction code is imported.  Coordinates mean the complex number
    a(exp(2*pi*i/N)) + sqrt(d)*b(exp(2*pi*i/N)),
where N is 42 or 84, d is 5 or 15, and all coefficients are rational.
All calculations use FLINT rational polynomials modulo the indicated exact
cyclotomic polynomial.  Complex conjugation sends t to t**(-1) and fixes
sqrt(d).  Q(sqrt(5)) and Q(sqrt(15)) have discriminants 5 and 60, respectively,
so the prime 5 ramifies in each.  Only primes dividing N can ramify in
Q(zeta_N); ramification in a subfield implies ramification in the extension.
Consequently sqrt(d) is not in Q(zeta_N).  The two coefficient polynomials are
therefore a unique representation, and conjugate multiplication is exact.

Non-6-colorability is accepted ONLY after an external DRAT checker verifies a
proof against the regenerated, unrestricted 6-coloring CNF.  The checker is a
trusted executable: its SHA-256 is reported so an independent reviewer can
identify/rebuild it.  A SAT solver's status is never accepted as a proof.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from flint import fmpq, fmpq_poly


class VerificationError(ValueError):
    """The supplied object does not establish the requested claim."""


# Coefficients in ascending powers.  Phi_84(t) = Phi_42(t**2).
PHI42 = [1, 1, 0, -1, -1, 0, 1, 0, -1, -1, 0, 1, 1]
PHI84 = [PHI42[i // 2] if i % 2 == 0 else 0 for i in range(25)]
RATIONAL = re.compile(r"[+-]?[0-9]+(?:/[+]?[0-9]+)?\Z")


def _integer(value: Any, label: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise VerificationError(f"{label} must be an integer >= {minimum}")
    return value


def _structure(graph: Any) -> tuple[int, list[tuple[int, int]]]:
    if not isinstance(graph, dict) or graph.get("schema") not in (
        "hn-exact-v1", "hn-multiquadratic-v1", "hn-cyclotomic-tower-v1", "hn-voronov-field-v1"
    ):
        raise VerificationError("expected a supported exact graph schema")
    vertices = graph.get("vertices")
    edges = graph.get("edges")
    if not isinstance(vertices, list) or not isinstance(edges, list):
        raise VerificationError("vertices and edges must be arrays")
    n = len(vertices)
    canonical_edges: set[tuple[int, int]] = set()
    for index, edge in enumerate(edges):
        if not isinstance(edge, list) or len(edge) != 2:
            raise VerificationError(f"edge {index} must have two endpoints")
        u = _integer(edge[0], f"edge {index} endpoint")
        v = _integer(edge[1], f"edge {index} endpoint")
        if u >= n or v >= n or u == v:
            raise VerificationError(f"edge {index} has an invalid endpoint or loop")
        pair = (min(u, v), max(u, v))
        if pair in canonical_edges:
            raise VerificationError(f"edge {index} is a duplicate")
        canonical_edges.add(pair)
    return n, sorted(canonical_edges)


class ExactField:
    def __init__(self, specification: Any):
        if not isinstance(specification, dict):
            raise VerificationError("field must be an object")
        self.order = specification.get("cyclotomic_order")
        self.radicand = specification.get("radicand")
        if type(self.order) is not int or self.order not in (42, 84):
            raise VerificationError("supported cyclotomic orders are 42 and 84")
        if type(self.radicand) is not int or self.radicand not in (5, 15):
            raise VerificationError("supported radicands are 5 and 15")
        self.modulus = fmpq_poly(PHI42 if self.order == 42 else PHI84)
        self.degree = self.modulus.degree()
        t = fmpq_poly([0, 1])
        # N is even, and t**N = 1 modulo Phi_N, hence t**(N-1) = 1/t.
        inverse_t = t ** (self.order - 1) % self.modulus
        self.conjugate_powers = [fmpq_poly([1])]
        for _ in range(1, self.degree):
            self.conjugate_powers.append(
                self.conjugate_powers[-1] * inverse_t % self.modulus
            )

    def parse(self, vertex: Any, index: int) -> tuple[fmpq_poly, fmpq_poly]:
        if not isinstance(vertex, dict):
            raise VerificationError(f"vertex {index} must be an object")
        result = []
        for part in ("a", "b"):
            values = vertex.get(part)
            if not isinstance(values, list) or len(values) != self.degree:
                raise VerificationError(
                    f"vertex {index}.{part} must contain {self.degree} rational strings"
                )
            coefficients = []
            for value in values:
                if not isinstance(value, str) or RATIONAL.fullmatch(value) is None:
                    raise VerificationError(f"vertex {index}.{part} has a non-rational string")
                try:
                    coefficients.append(fmpq(value))
                except (ValueError, ZeroDivisionError) as exc:
                    raise VerificationError(f"vertex {index}.{part} has an invalid rational") from exc
            result.append(fmpq_poly(coefficients))
        return result[0], result[1]

    def conjugate(self, polynomial: fmpq_poly) -> fmpq_poly:
        result = fmpq_poly([])
        for power, coefficient in enumerate(polynomial.coeffs()):
            if coefficient:
                result += coefficient * self.conjugate_powers[power]
        return result

    def unit_difference(
        self,
        left: tuple[fmpq_poly, fmpq_poly],
        right: tuple[fmpq_poly, fmpq_poly],
    ) -> bool:
        a, b = left[0] - right[0], left[1] - right[1]
        ac, bc = self.conjugate(a), self.conjugate(b)
        rational_part = (a * ac + self.radicand * b * bc) % self.modulus
        radical_part = (a * bc + b * ac) % self.modulus
        return rational_part == 1 and radical_part == 0


@dataclass(frozen=True)
class GeometryResult:
    vertices: int
    edges: int
    cyclotomic_order: int
    radicand: int


def verify_geometry(graph: Any) -> GeometryResult:
    if not isinstance(graph, dict) or graph.get("schema") != "hn-exact-v1":
        raise VerificationError("cyclotomic geometry requires schema hn-exact-v1")
    n, edges = _structure(graph)
    field = ExactField(graph.get("field"))
    points = [field.parse(vertex, i) for i, vertex in enumerate(graph["vertices"])]
    seen: dict[tuple[tuple[fmpq, ...], tuple[fmpq, ...]], int] = {}
    for index, point in enumerate(points):
        key = (tuple(point[0].coeffs()), tuple(point[1].coeffs()))
        if key in seen:
            raise VerificationError(f"vertices {seen[key]} and {index} are the same exact point")
        seen[key] = index
    for u, v in edges:
        if not field.unit_difference(points[u], points[v]):
            raise VerificationError(f"edge ({u}, {v}) does not have exact squared length 1")
    return GeometryResult(n, len(edges), field.order, field.radicand)


def coloring_cnf(graph: Any, k: int = 6) -> bytes:
    """Canonical DIMACS: x_(v,c) = v*k+c+1; colors c are zero-based.

    Clauses: vertex at-least-one, then all its pairwise at-most-one clauses;
    then, for sorted undirected edges, the k edge inequalities.  There are no
    symmetry-breaking constraints, preprocessing assumptions, or fixed colors.
    """
    _integer(k, "k", 1)
    n, edges = _structure(graph)
    clauses = n * (1 + k * (k - 1) // 2) + k * len(edges)
    lines = [f"p cnf {n * k} {clauses}\n"]
    for vertex in range(n):
        variables = [vertex * k + color + 1 for color in range(k)]
        lines.append(" ".join(map(str, variables)) + " 0\n")
        lines.extend(f"{-a} {-b} 0\n" for a, b in itertools.combinations(variables, 2))
    for u, v in edges:
        lines.extend(f"{-(u*k+c+1)} {-(v*k+c+1)} 0\n" for c in range(k))
    return "".join(lines).encode("ascii")


def verify_coloring(graph: Any, coloring: Any, k: int | None = None) -> dict[str, int]:
    n, edges = _structure(graph)
    if not isinstance(coloring, dict):
        raise VerificationError("coloring must be an object with k and colors")
    actual_k = _integer(coloring.get("k"), "coloring.k", 1)
    if k is not None and actual_k != k:
        raise VerificationError("coloring.k does not match the requested k")
    colors = coloring.get("colors")
    if not isinstance(colors, list) or len(colors) != n:
        raise VerificationError("coloring must assign a color to every vertex")
    for index, color in enumerate(colors):
        _integer(color, f"color of vertex {index}")
        if color >= actual_k:
            raise VerificationError(f"color of vertex {index} is outside [0,k)")
    for u, v in edges:
        if colors[u] == colors[v]:
            raise VerificationError(f"edge ({u}, {v}) is monochromatic")
    return {"k": actual_k, "vertices": n, "edges": len(edges)}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_coloring_proof(
    graph: Any,
    cnf_path: str | Path,
    proof_path: str | Path,
    checker_path: str | Path,
    k: int = 6,
    timeout: float | None = None,
) -> dict[str, Any]:
    """Bind/check non-k-colorability; does not establish geometry or success.

    The public field-specific certificate wrappers must first verify geometry.
    """
    _integer(k, "proof color count", 1)
    expected = coloring_cnf(graph, k)
    cnf = Path(cnf_path).resolve(strict=True)
    proof = Path(proof_path).resolve(strict=True)
    checker = Path(checker_path).resolve(strict=True)
    if cnf.read_bytes() != expected:
        raise VerificationError(f"proof CNF bytes do not match regenerated unrestricted {k}-coloring CNF")
    checker_sha256 = _sha256(checker)
    with tempfile.TemporaryDirectory(prefix="hn-independent-check-") as temporary:
        directory = Path(temporary)
        frozen_cnf = directory / "coloring.cnf"
        frozen_proof = directory / "proof.drat"
        frozen_cnf.write_bytes(expected)
        shutil.copyfile(proof, frozen_proof)
        proof_sha256 = _sha256(frozen_proof)
        try:
            completed = subprocess.run(
                [str(checker), str(frozen_cnf), str(frozen_proof)],
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise VerificationError(f"external proof checker did not complete: {exc}") from exc
        if completed.returncode != 0 or re.search(
            r"^s VERIFIED\s*$", completed.stdout, flags=re.MULTILINE
        ) is None:
            raise VerificationError(
                f"external DRAT proof check failed (exit {completed.returncode}): "
                + (completed.stdout + completed.stderr)[-2000:]
            )
        if frozen_cnf.read_bytes() != expected or _sha256(frozen_proof) != proof_sha256:
            raise VerificationError("checker inputs changed during verification")
    if _sha256(checker) != checker_sha256:
        raise VerificationError("proof checker executable changed during verification")
    return {
        "verified_color_count": k,
        "cnf_sha256": hashlib.sha256(expected).hexdigest(),
        "proof_sha256": proof_sha256,
        "checker_sha256": checker_sha256,
        "checker_path": str(checker),
        "checker_output": completed.stdout,
    }


def _verify_non6_proof(
    graph: Any,
    cnf_path: str | Path,
    proof_path: str | Path,
    checker_path: str | Path,
    k: int = 6,
    timeout: float | None = None,
) -> dict[str, Any]:
    if type(k) is not int or k != 6:
        raise VerificationError("the target certificate must establish non-6-colorability")
    return _verify_coloring_proof(graph, cnf_path, proof_path, checker_path, k, timeout)


def verify_non6_certificate(
    graph: Any,
    cnf_path: str | Path,
    proof_path: str | Path,
    checker_path: str | Path,
    k: int = 6,
    timeout: float | None = None,
) -> dict[str, Any]:
    """Verify cyclotomic geometry AND independently check its non-6 proof."""
    if type(k) is not int or k != 6:
        raise VerificationError("the target certificate must establish non-6-colorability")
    geometry = verify_geometry(graph)
    proof = _verify_non6_proof(graph, cnf_path, proof_path, checker_path, k, timeout)
    return {
        "goal_satisfied": True,
        "certificate_status": "exact-unit-distance-and-non-6-colorable",
        "geometry": geometry.__dict__,
        **proof,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph", required=True, type=Path)
    parser.add_argument("--k", type=int, default=6)
    parser.add_argument("--coloring", type=Path)
    parser.add_argument("--write-cnf", type=Path)
    parser.add_argument("--cnf", type=Path)
    parser.add_argument("--proof", type=Path)
    parser.add_argument("--checker", type=Path)
    parser.add_argument("--checker-timeout", type=float)
    args = parser.parse_args(argv)
    try:
        graph = json.loads(args.graph.read_text())
        proof_arguments = (args.cnf, args.proof, args.checker)
        if any(proof_arguments) and not all(proof_arguments):
            raise VerificationError("a non-6 certificate requires --cnf, --proof, and --checker together")
        if all(proof_arguments):
            report = verify_non6_certificate(
                graph, args.cnf, args.proof, args.checker, args.k, args.checker_timeout
            )
        else:
            report = {
                "goal_satisfied": False,
                "certificate_status": "geometry-only",
                "geometry": verify_geometry(graph).__dict__,
            }
        if args.coloring:
            report["coloring"] = verify_coloring(graph, json.loads(args.coloring.read_text()), args.k)
            if report["goal_satisfied"]:
                raise VerificationError("a valid 6-coloring contradicts the claimed non-6 certificate")
            report["certificate_status"] = "geometry-and-valid-coloring"
        if args.write_cnf:
            canonical = coloring_cnf(graph, args.k)
            args.write_cnf.write_bytes(canonical)
            report["cnf_sha256"] = hashlib.sha256(canonical).hexdigest()
        report["graph_sha256"] = _sha256(args.graph)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    except (VerificationError, OSError, ValueError) as exc:
        print(json.dumps({"goal_satisfied": False, "error": str(exc)}, sort_keys=True), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
