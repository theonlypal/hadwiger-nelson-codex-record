"""Independent multiquadratic geometry and certificate-boundary tests."""

import itertools
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import verify_mq as verify


def point(x=None, y=None):
    result = {"x": ["0"]*16, "y": ["0"]*16}
    for axis, coefficients in (("x", x or {}), ("y", y or {})):
        for mask, coefficient in coefficients.items():
            result[axis][mask] = str(coefficient)
    return result


def triangle():
    return {
        "schema": "hn-multiquadratic-v1",
        "radicands": [2, 3, 5, 11],
        "vertices": [point(), point({0: 1}), point({0: "1/2"}, {2: "1/2"})],
        "edges": [[0, 1], [0, 2], [1, 2]],
        "metadata": {},
    }


class MultiquadraticGeometryTests(unittest.TestCase):
    def test_triangle(self):
        self.assertEqual(verify.verify_geometry(triangle()).edges, 3)

    def test_radical_unit_vectors(self):
        units = [
            point({0: "7/8"}, {6: "1/8"}),  # sqrt(15)
            point({0: "5/6"}, {8: "1/6"}),  # sqrt(11)
            point({1: "1/2"}, {1: "1/2"}),  # sqrt(2)
            point({1: "1/4", 3: "1/4"}, {3: "1/4", 1: "-1/4"}),
        ]
        for unit in units:
            graph = triangle()
            graph["vertices"] = [point(), unit]
            graph["edges"] = [[0, 1]]
            with self.subTest(unit=unit):
                self.assertEqual(verify.verify_geometry(graph).edges, 1)

    def test_independent_dense_norm_identity(self):
        # A translated 45-degree unit vector in a dense rational basis tests
        # subtraction, denominator normalization, and cancellation.
        graph = triangle()
        left_x = {mask: f"{mask+1}/30" for mask in range(16)}
        left_y = {mask: f"{mask-5}/42" for mask in range(16)}
        from fractions import Fraction
        right_x = dict(left_x); right_y = dict(left_y)
        right_x[1] = str(Fraction(left_x[1])+Fraction(1, 2))
        right_y[1] = str(Fraction(left_y[1])-Fraction(1, 2))
        graph["vertices"] = [point(left_x, left_y), point(right_x, right_y)]
        graph["edges"] = [[0, 1]]
        self.assertEqual(verify.verify_geometry(graph).edges, 1)

    def test_nonconstant_norm_terms_cannot_be_ignored(self):
        # x=(sqrt(2)+sqrt(3))/sqrt(5) has constant norm 1 but
        # squared norm 1 + 2*sqrt(6)/5, so it is not a unit vector.
        graph = triangle()
        graph["vertices"] = [point(), point({5: "1/5", 6: "1/5"})]
        graph["edges"] = [[0, 1]]
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_geometry(graph)

    def test_duplicate_points(self):
        graph = triangle()
        graph["vertices"].append(point({0: "2/2"}))
        with self.assertRaisesRegex(verify.VerificationError, "same exact point"):
            verify.verify_geometry(graph)

    def test_distinct_radicals(self):
        graph = triangle()
        graph["vertices"] = [point({1: 1}), point({2: 1}), point({3: 1})]
        graph["edges"] = []
        self.assertEqual(verify.verify_geometry(graph).vertices, 3)

    def test_wrong_radical_basis_rejected(self):
        for radicands in ([2, 3, 5, 15], [3, 2, 5, 11], [2, 3, 5, 11.0], [2, 3, 5]):
            graph = triangle()
            graph["radicands"] = radicands
            with self.subTest(radicands=radicands), self.assertRaises(verify.VerificationError):
                verify.verify_geometry(graph)

    def test_bad_rationals(self):
        for value in (1, 0.5, "0.5", "1/0", "nan"):
            graph = triangle()
            graph["vertices"][0]["x"][0] = value
            with self.subTest(value=value), self.assertRaises(verify.VerificationError):
                verify.verify_geometry(graph)

    def test_numerically_near_unit_rejected(self):
        graph = triangle()
        graph["vertices"][1]["x"][0] = "1000000000000000000000001/1000000000000000000000000"
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_geometry(graph)

    def test_cross_schema_geometry_rejected(self):
        graph = triangle()
        with self.assertRaises(verify.VerificationError):
            verify.common.verify_geometry(graph)
        graph["schema"] = "hn-exact-v1"
        with self.assertRaises(verify.VerificationError):
            verify.verify_geometry(graph)


class MultiquadraticCertificateTests(unittest.TestCase):
    def test_non6_k7_bad_geometry_rejected(self):
        graph = triangle()
        graph["vertices"] = [point({0: index}) for index in range(7)]
        graph["edges"] = [list(edge) for edge in itertools.combinations(range(7), 2)]
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_non6_certificate(graph, "none", "none", "none")

    def test_wrong_k_rejected(self):
        with self.assertRaisesRegex(verify.VerificationError, "non-6"):
            verify.verify_non6_certificate(triangle(), "none", "none", "none", k=5)

    def test_proper_coloring(self):
        self.assertEqual(verify.verify_coloring(triangle(), {"k": 3, "colors": [0, 1, 2]})["k"], 3)
        with self.assertRaisesRegex(verify.VerificationError, "monochromatic"):
            verify.verify_coloring(triangle(), {"k": 2, "colors": [0, 1, 0]})

    def test_canonical_cnf_matches_known_triangle(self):
        from tools import verify as common
        graph = triangle()
        abstract = dict(graph, schema="hn-exact-v1")
        self.assertEqual(verify.coloring_cnf(graph, 6), common.coloring_cnf(abstract, 6))

    def test_forged_cnf_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            cnf, proof = directory/"coloring.cnf", directory/"proof.drat"
            cnf.write_bytes(verify.coloring_cnf(triangle(), 6)+b"0\n")
            proof.write_bytes(b"0\n")
            with self.assertRaisesRegex(verify.VerificationError, "CNF bytes"):
                verify.verify_non6_certificate(triangle(), cnf, proof, sys.executable)

    def test_solver_status_not_accepted_as_proof(self):
        from subprocess import CompletedProcess
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            cnf, proof = directory/"coloring.cnf", directory/"proof.drat"
            cnf.write_bytes(verify.coloring_cnf(triangle(), 6)); proof.write_bytes(b"0\n")
            with patch.object(verify.common.subprocess, "run", return_value=CompletedProcess([], 0, "s UNSATISFIABLE\n", "")):
                with self.assertRaises(verify.VerificationError):
                    verify.verify_non6_certificate(triangle(), cnf, proof, sys.executable)


if __name__ == "__main__":
    unittest.main()
