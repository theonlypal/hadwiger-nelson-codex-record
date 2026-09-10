from pathlib import Path
import subprocess
import tempfile
import unittest

from tools.native_pair_diagnostic import ROOT, check_conditional_proof, conditional_cnf
from tools import verify_tower


DRIVER = ROOT/"sources/cadical/build/cadical-phase"
CHECKER = ROOT/"sources/drat-trim/drat-trim"


def path_graph():
    def vertex(n):
        parts = [["0"]*12 for _ in range(4)]; parts[0][0] = str(n)
        return {"parts": parts}
    return {"schema": "hn-cyclotomic-tower-v1",
            "field": {"cyclotomic_order": 42, "radicands": [5, 33]},
            "vertices": [vertex(0), vertex(1), vertex(2)], "edges": [[0, 1], [1, 2]]}


@unittest.skipUnless(DRIVER.exists() and CHECKER.exists(), "native checker tools must be built")
class ConditionalPairGateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name); self.graph = path_graph()
        verify_tower.verify_geometry(self.graph)
        verify_tower.verify_coloring(self.graph, {"k": 2, "colors": [0, 1, 0]}, 2)
        self.cnf, self.proof = self.root/"conditional.cnf", self.root/"proof.drat"
        _, expected, _ = conditional_cnf(self.graph, (0, 2), 2)
        self.cnf.write_bytes(expected)
        result = subprocess.run([str(DRIVER), "--unsat", "-w", str(self.root/"solution.out"),
                                 str(self.cnf), str(self.proof)], capture_output=True, timeout=10)
        self.assertEqual(result.returncode, 20)

    def test_actual_unit_path_conditional_proof(self):
        result, _, _ = check_conditional_proof(self.graph, (0, 2), self.cnf, self.proof, CHECKER, k=2)
        self.assertTrue(result["forced_equality_established"])
        self.assertFalse(result["goal_satisfied"])
        self.assertEqual(result["color_count"], 2)
        self.assertEqual(result["added_inequality_clauses"], [[-1, -5], [-2, -6]])

    def test_extra_constraint_rejected_before_checker(self):
        data = self.cnf.read_bytes().replace(b"p cnf 6 12\n", b"p cnf 6 13\n")+b"1 0\n"
        self.cnf.write_bytes(data)
        with self.assertRaisesRegex(ValueError, "exact independent regeneration"):
            check_conditional_proof(self.graph, (0, 2), self.cnf, self.proof, CHECKER, k=2)

    def test_wrong_pair_rejected_before_checker(self):
        with self.assertRaisesRegex(ValueError, "exact independent regeneration"):
            check_conditional_proof(self.graph, (0, 1), self.cnf, self.proof, CHECKER, k=2)

    def test_empty_proof_does_not_establish_pair(self):
        self.proof.write_bytes(b"")
        result, _, _ = check_conditional_proof(self.graph, (0, 2), self.cnf, self.proof, CHECKER, k=2)
        self.assertFalse(result["forced_equality_established"])
        self.assertFalse(result["goal_satisfied"])


if __name__ == "__main__":
    unittest.main()
