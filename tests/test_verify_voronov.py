import copy
from fractions import Fraction
import itertools
from pathlib import Path
import tempfile
import unittest

from tools import verify_voronov as verify


def vertex(x=None, y=None):
    x = {} if x is None else x; y = {} if y is None else y
    return {"x": [str(x.get(i, 0)) for i in range(16)],
            "y": [str(y.get(i, 0)) for i in range(16)]}


def graph(vertices=None, edges=None):
    return {"schema": "hn-voronov-field-v1", "field": copy.deepcopy(verify.FIELD),
            "vertices": [vertex(), vertex({0: 1}), vertex({0: "1/2"}, {2: "1/2"})] if vertices is None else vertices,
            "edges": [[0, 1], [0, 2], [1, 2]] if edges is None else edges}


class VoronovGeometryTests(unittest.TestCase):
    def test_degree_proof_exact_signs(self):
        proof = verify.field_basis_certificate()
        self.assertEqual(proof["degree"], 16)
        self.assertEqual(proof["embedding_signs_h1_h7_product"],
                         [[1, 1, 1], [-1, 1, -1], [1, 1, 1], [1, -1, -1]])
        self.assertEqual(proof["h1_h7_product"], [-14, 48, 0, 0])

    def test_extension_product_basis(self):
        self.assertEqual(dict(verify.PRODUCTS[4][4]), {0: 5, 1: 3, 2: -1, 3: 3})
        self.assertEqual(dict(verify.PRODUCTS[8][8]), {0: 5, 1: 3, 2: 1, 3: -3})
        self.assertEqual(dict(verify.PRODUCTS[4][8]), {12: 1})
        self.assertEqual(dict(verify.PRODUCTS[12][12]), {0: -14, 1: 48})

    def test_exact_base_triangle(self):
        result = verify.verify_geometry(graph())
        self.assertEqual((result.vertices, result.edges, result.field_degree), (3, 3, 16))

    def test_units_using_each_extension(self):
        # Exact rational parametrization x=(1-h)/(1+h), y=2sqrt(h)/(1+h).
        # Coefficients were derived independently with symbolic rationalization.
        h1 = vertex({0: "-571/647", 1: "74/647", 2: "-214/1941", 3: "10/647"},
                    {4: "76/647", 5: "74/647", 6: "-214/1941", 7: "10/647"})
        h7 = vertex({0: "-571/647", 1: "74/647", 2: "214/1941", 3: "-10/647"},
                    {8: "76/647", 9: "74/647", 10: "214/1941", 11: "-10/647"})
        product = vertex({0: "-4413/4439", 1: "96/4439"}, {12: "26/4439", 13: "96/4439"})
        result = verify.verify_geometry(graph([vertex(), h1, h7, product], [[0, 1], [0, 2], [0, 3]]))
        self.assertEqual(result.edges, 3)

    def test_near_unit_edge_rejected_exactly(self):
        bad = graph(); bad["vertices"][1]["x"][0] = "1000000000000001/1000000000000000"
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_geometry(bad)

    def test_nonsquare_radical_is_not_falsely_unit(self):
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_geometry(graph([vertex(), vertex({4: 1})], [[0, 1]]))

    def test_equivalent_rational_duplicate_rejected(self):
        with self.assertRaisesRegex(verify.VerificationError, "same exact point"):
            verify.verify_geometry(graph([vertex({0: "1/2"}), vertex({0: "2/4"})], []))

    def test_bad_field_metadata_rejected(self):
        variants = []
        value = graph(); value["field"]["base_radicands"] = [3, 2]; variants.append(value)
        value = graph(); value["field"]["extension_squares"].reverse(); variants.append(value)
        value = graph(); value["field"]["extension_squares"][1][2] = True; variants.append(value)
        value = graph(); value["field"]["embedding"] = "other"; variants.append(value)
        for value in variants:
            with self.subTest(value=value["field"]), self.assertRaises(verify.VerificationError):
                verify.verify_geometry(value)

    def test_wrong_geometry_schema_rejected(self):
        value = graph(); value["schema"] = "hn-multiquadratic-v1"
        with self.assertRaisesRegex(verify.VerificationError, "Voronov geometry"):
            verify.verify_geometry(value)

    def test_malformed_coordinate_rejected(self):
        for coordinate in (["0"]*15, ["0"]*15+["1/0"], ["0"]*15+["sqrt(2)"], ["0"]*15+[1]):
            value = graph(); value["vertices"][0]["x"] = coordinate
            with self.subTest(coordinate=coordinate), self.assertRaises(verify.VerificationError):
                verify.verify_geometry(value)

    def test_bad_edges_rejected(self):
        for edges in ([[0, 0]], [[0, 3]], [[0, 1], [1, 0]]):
            with self.subTest(edges=edges), self.assertRaises(verify.VerificationError):
                verify.verify_geometry(graph(edges=edges))

    def test_proper_coloring_checked(self):
        self.assertEqual(verify.verify_coloring(graph(), {"k": 5, "colors": [0, 1, 2]}, 5)["edges"], 3)
        with self.assertRaises(verify.VerificationError):
            verify.verify_coloring(graph(), {"k": 5, "colors": [0, 1, 0]}, 5)

    def test_unrestricted_cnf_matches_geometry_agnostic_encoding(self):
        value = graph(); control = copy.deepcopy(value); control["schema"] = "hn-exact-v1"
        self.assertEqual(verify.coloring_cnf(value, 6), verify.common.coloring_cnf(control, 6))

    def test_public_gate_rejects_non6_count(self):
        with self.assertRaisesRegex(verify.VerificationError, "non-6"):
            verify.verify_non6_certificate(graph(), "missing", "missing", "missing", k=5)

    def test_k7_fake_geometry_cannot_reach_success_gate(self):
        value = graph([vertex({0: i}) for i in range(7)], [list(pair) for pair in itertools.combinations(range(7), 2)])
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_non6_certificate(value, "missing", "missing", "missing")

    def test_wrong_proof_cnf_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            cnf, proof = Path(directory)/"input.cnf", Path(directory)/"proof.drat"
            cnf.write_bytes(verify.coloring_cnf(graph(), 5)); proof.write_bytes(b"0\n")
            with self.assertRaises(verify.VerificationError):
                verify.verify_non6_certificate(graph(), cnf, proof, Path(__file__))


if __name__ == "__main__":
    unittest.main()
