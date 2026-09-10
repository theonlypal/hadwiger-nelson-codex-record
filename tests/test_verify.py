"""Verifier checks intentionally independent of hn construction arithmetic."""

import copy
import importlib.util
import itertools
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


MODULE = Path(__file__).resolve().parents[1] / "tools" / "verify.py"
SPEC = importlib.util.spec_from_file_location("independent_verify", MODULE)
verify = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = verify
SPEC.loader.exec_module(verify)


def point(a=None, b=None, degree=12):
    result = {"a": ["0"] * degree, "b": ["0"] * degree}
    for name, coefficients in (("a", a or {}), ("b", b or {})):
        for power, value in coefficients.items():
            result[name][power] = str(value)
    return result


def triangle():
    # 0, 1, exp(i*pi/3); zeta_42**7 = exp(i*pi/3).
    return {
        "schema": "hn-exact-v1",
        "field": {"cyclotomic_order": 42, "radicand": 5},
        "vertices": [point(), point({0: 1}), point({7: 1})],
        "edges": [[0, 1], [0, 2], [1, 2]],
        "metadata": {},
    }


class GeometryTests(unittest.TestCase):
    def test_exact_triangle(self):
        self.assertEqual(verify.verify_geometry(triangle()).edges, 3)

    def test_radical_rotation_is_unit(self):
        # 7/8 + i*sqrt(15)/8 has squared norm (49+15)/64 = 1.
        # Here i*sqrt(3) = 2*zeta_42**7 - 1, so the sqrt(5)
        # coefficient is (2*t^7-1)/8.
        graph = triangle()
        graph["vertices"] = [point(), point({0: "7/8"}, {0: "-1/8", 7: "1/4"})]
        graph["edges"] = [[0, 1]]
        self.assertEqual(verify.verify_geometry(graph).edges, 1)

    def test_order84(self):
        graph = triangle()
        graph["field"] = {"cyclotomic_order": 84, "radicand": 15}
        graph["vertices"] = [point(degree=24), point({0: 1}, degree=24), point({14: 1}, degree=24)]
        self.assertEqual(verify.verify_geometry(graph).edges, 3)

    def test_bad_edge_rejected(self):
        graph = triangle()
        graph["vertices"][1] = point({0: 2})
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_geometry(graph)

    def test_numerically_near_unit_rejected(self):
        graph = triangle()
        graph["vertices"][1] = point({0: "1000000000000000000000001/1000000000000000000000000"})
        with self.assertRaises(verify.VerificationError):
            verify.verify_geometry(graph)

    def test_duplicate_rational_point_rejected(self):
        graph = triangle()
        graph["vertices"].append(point({0: "2/2"}))
        with self.assertRaisesRegex(verify.VerificationError, "same exact point"):
            verify.verify_geometry(graph)

    def test_distinct_radical_points(self):
        graph = triangle()
        graph["vertices"] = [point({0: 1}), point(b={0: 1})]
        graph["edges"] = []
        self.assertEqual(verify.verify_geometry(graph).vertices, 2)

    def test_bad_inputs(self):
        for replacement in ([0, 0], [0, 4], [False, 1], [0, 1, 2]):
            with self.subTest(replacement=replacement):
                graph = triangle()
                graph["edges"] = [replacement]
                with self.assertRaises(verify.VerificationError):
                    verify.verify_geometry(graph)
        graph = triangle()
        graph["edges"].append([1, 0])
        with self.assertRaisesRegex(verify.VerificationError, "duplicate"):
            verify.verify_geometry(graph)
        for bad in (1, 0.5, "0.5", "1/0", "nan"):
            graph = triangle()
            graph["vertices"][0]["a"][0] = bad
            with self.subTest(rational=bad), self.assertRaises(verify.VerificationError):
                verify.verify_geometry(graph)

    def test_unjustified_field_extension_rejected(self):
        # sqrt(3) cannot use the prime-5 ramification argument on which this
        # verifier's unique-pair representation is based.
        graph = triangle()
        graph["field"]["radicand"] = 3
        with self.assertRaisesRegex(verify.VerificationError, "radicands"):
            verify.verify_geometry(graph)


class ColoringTests(unittest.TestCase):
    def test_coloring(self):
        self.assertEqual(verify.verify_coloring(triangle(), {"k": 3, "colors": [0, 1, 2]})["k"], 3)

    def test_invalid_coloring(self):
        for coloring in (
            {"k": 2, "colors": [0, 1, 0]},
            {"k": 3, "colors": [0, 1]},
            {"k": 3, "colors": [0, 1, 3]},
            {"k": 3, "colors": [0, 1, True]},
        ):
            with self.subTest(coloring=coloring), self.assertRaises(verify.VerificationError):
                verify.verify_coloring(triangle(), coloring)

    def test_cnf_canonical_edge_order(self):
        graph = triangle()
        other = copy.deepcopy(graph)
        other["edges"] = [[2, 1], [2, 0], [1, 0]]
        self.assertEqual(verify.coloring_cnf(graph, 3), verify.coloring_cnf(other, 3))

    def test_cnf_equivalent_to_coloring_exhaustively(self):
        # Check all 2^(3*k) truth assignments, not only one-hot assignments.
        for k in (2, 3):
            lines = verify.coloring_cnf(triangle(), k).decode().splitlines()
            clauses = [list(map(int, line.split()))[:-1] for line in lines[1:]]
            satisfying = []
            for assignment in itertools.product((False, True), repeat=3*k):
                sat = all(any(assignment[abs(literal)-1] == (literal > 0) for literal in clause) for clause in clauses)
                if sat:
                    satisfying.append(assignment)
                    colors = []
                    for v in range(3):
                        selected = [c for c in range(k) if assignment[v*k+c]]
                        self.assertEqual(len(selected), 1)
                        colors.append(selected[0])
                    verify.verify_coloring(triangle(), {"k": k, "colors": colors})
            self.assertEqual(len(satisfying), 0 if k == 2 else 6)


class CertificateTests(unittest.TestCase):
    def test_non6_abstract_graph_with_bad_geometry_rejected(self):
        graph = triangle()
        graph["vertices"] = [point({0: i}) for i in range(7)]
        graph["edges"] = [list(edge) for edge in itertools.combinations(range(7), 2)]
        # K7 has no 6-coloring, but these collinear coordinates do not realize
        # its edges. The success gate must reject before reading a proof.
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_non6_certificate(graph, "none", "none", "none")

    def test_wrong_k_rejected_before_checker(self):
        with self.assertRaisesRegex(verify.VerificationError, "non-6"):
            verify.verify_non6_certificate(triangle(), "none", "none", "none", k=5)

    def test_wrong_cnf_rejected_before_checker(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            cnf, proof = root / "input.cnf", root / "proof.drat"
            cnf.write_bytes(verify.coloring_cnf(triangle(), 6) + b"0\n")
            proof.write_bytes(b"")
            with self.assertRaisesRegex(verify.VerificationError, "CNF bytes"):
                verify.verify_non6_certificate(triangle(), cnf, proof, sys.executable)

    def test_solver_status_is_not_proof(self):
        from subprocess import CompletedProcess
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            cnf, proof = root / "input.cnf", root / "proof.drat"
            cnf.write_bytes(verify.coloring_cnf(triangle(), 6))
            proof.write_bytes(b"0\n")
            for code, output in ((0, "s UNSATISFIABLE\n"), (1, "s VERIFIED\n"), (0, "c s VERIFIED\n")):
                with patch.object(verify.subprocess, "run", return_value=CompletedProcess([], code, output, "")):
                    with self.subTest(code=code, output=output), self.assertRaises(verify.VerificationError):
                        verify.verify_non6_certificate(triangle(), cnf, proof, sys.executable)


if __name__ == "__main__":
    unittest.main()
