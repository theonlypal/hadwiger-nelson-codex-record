import copy
import itertools
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools import verify_tower as verify
from tools import verify_mq


def point(parts=None):
    result = {"parts": [["0"]*12 for _ in range(4)]}
    for part, coefficients in (parts or {}).items():
        for power, value in coefficients.items():
            result["parts"][part][power] = str(value)
    return result


def triangle():
    return {"schema": "hn-cyclotomic-tower-v1", "field": {"cyclotomic_order": 42, "radicands": [5, 33]},
            "vertices": [point(), point({0: {0: 1}}), point({0: {7: 1}})],
            "edges": [[0, 1], [0, 2], [1, 2]], "metadata": {}}


class TowerGeometryTests(unittest.TestCase):
    def test_equilateral_triangle(self):
        self.assertEqual(verify.verify_geometry(triangle()).edges, 3)

    def test_two_roots_and_their_product(self):
        rotations = [
            point({0: {0: "7/8"}, 1: {0: "-1/8", 7: "1/4"}}),
            point({0: {0: "5/6"}, 2: {0: "-1/18", 7: "1/9"}}),
            point({0: {0: "35/48"}, 1: {0: "-5/48", 7: "5/24"},
                   2: {0: "-7/144", 7: "7/72"}, 3: {0: "-1/48"}}),
        ]
        for rotation in rotations:
            graph = triangle(); graph["vertices"] = [point(), rotation]; graph["edges"] = [[0, 1]]
            with self.subTest(rotation=rotation):
                self.assertEqual(verify.verify_geometry(graph).edges, 1)

    def test_bad_norm(self):
        graph = triangle(); graph["vertices"][1] = point({2: {0: 1}})
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_geometry(graph)

    def test_nonzero_cross_terms_rejected(self):
        # z=(1+sqrt5)/sqrt6 is unavailable directly; this intentionally
        # perturbs a genuine four-part unit vector in its product part.
        graph = triangle(); graph["vertices"][1] = point({0: {0: "35/48"}, 1: {0: "-5/48", 7: "5/24"},
                   2: {0: "-7/144", 7: "7/72"}, 3: {0: "1/48"}})
        with self.assertRaises(verify.VerificationError):
            verify.verify_geometry(graph)

    def test_duplicate_exact_points(self):
        graph = triangle(); graph["vertices"].append(point({0: {0: "2/2"}}))
        with self.assertRaisesRegex(verify.VerificationError, "same exact point"):
            verify.verify_geometry(graph)

    def test_parts_are_independent(self):
        graph = triangle(); graph["vertices"] = [point({i: {0: 1}}) for i in range(4)]; graph["edges"] = []
        self.assertEqual(verify.verify_geometry(graph).vertices, 4)

    def test_bad_part_shapes(self):
        for parts in ([["0"]*12]*3, [["0"]*11]*4, [0]*4):
            graph = triangle(); graph["vertices"][0]["parts"] = parts
            with self.subTest(parts=parts), self.assertRaises(verify.VerificationError):
                verify.verify_geometry(graph)

    def test_field_cannot_change(self):
        for field in ({"cyclotomic_order": 84, "radicands": [5, 33]},
                      {"cyclotomic_order": 42, "radicands": [5, 15]},
                      {"cyclotomic_order": 42, "radicands": [33, 5]}):
            graph = triangle(); graph["field"] = field
            with self.subTest(field=field), self.assertRaises(verify.VerificationError):
                verify.verify_geometry(graph)

    def test_existing_haugland_lift(self):
        path = Path(__file__).resolve().parents[1]/"artifacts/haugland_g3.json"
        if not path.exists():
            self.skipTest("Haugland control artifact not present")
        original = json.loads(path.read_text()); lifted = triangle()
        lifted["vertices"] = [{"parts": [p["a"], p["b"], ["0"]*12, ["0"]*12]} for p in original["vertices"]]
        lifted["edges"] = original["edges"]
        result = verify.verify_geometry(lifted)
        self.assertEqual((result.vertices, result.edges), (2131, 12530))

    def test_cross_schema_rejected(self):
        for other in (verify.common, verify_mq):
            with self.assertRaises(verify.VerificationError):
                other.verify_geometry(triangle())
        graph = triangle(); graph["schema"] = "hn-exact-v1"
        with self.assertRaises(verify.VerificationError):
            verify.verify_geometry(graph)


class TowerGateTests(unittest.TestCase):
    def test_wrong_k_rejected(self):
        with self.assertRaisesRegex(verify.VerificationError, "non-6"):
            verify.verify_non6_certificate(triangle(), "none", "none", "none", k=5)

    def test_invalid_k7_geometry_rejected(self):
        graph = triangle(); graph["vertices"] = [point({0: {0: i}}) for i in range(7)]
        graph["edges"] = [list(edge) for edge in itertools.combinations(range(7), 2)]
        with self.assertRaisesRegex(verify.VerificationError, "exact squared length"):
            verify.verify_non6_certificate(graph, "none", "none", "none")

    def test_coloring_and_cnf(self):
        graph = triangle()
        verify.verify_coloring(graph, {"k": 3, "colors": [0, 1, 2]}, 3)
        with self.assertRaisesRegex(verify.VerificationError, "monochromatic"):
            verify.verify_coloring(graph, {"k": 3, "colors": [0, 1, 0]}, 3)
        abstract = dict(graph, schema="hn-exact-v1")
        self.assertEqual(verify.coloring_cnf(graph, 6), verify.common.coloring_cnf(abstract, 6))


if __name__ == "__main__":
    unittest.main()
