from pathlib import Path
import subprocess
import tempfile
import unittest

from hn.solve_external import DEFAULT_CHECKER, DEFAULT_SOLVER, parse_model, verify_full_cnf_model

DRIVER = DEFAULT_SOLVER.with_name('cadical-phase')


@unittest.skipUnless(DRIVER.exists(), 'Build the phase-hint API driver first')
class PhaseDriverTests(unittest.TestCase):
    def run_formula(self, folder, formula, phases):
        folder = Path(folder)
        cnf, phase = folder/'input.cnf', folder/'phases.txt'
        cnf.write_text(formula); phase.write_text(phases)
        proof, model = folder/'proof.drat', folder/'solution.out'
        result = subprocess.run([str(DRIVER), '--sat', '-w', str(model),
             str(cnf), str(proof), '--phase-file', str(phase)], capture_output=True, timeout=10)
        return result, cnf, proof, model

    def test_wrong_hint_does_not_become_a_constraint(self):
        with tempfile.TemporaryDirectory() as folder:
            result, cnf, proof, model = self.run_formula(folder, 'p cnf 1 1\n1 0\n', '-1 0\n')
            self.assertEqual(result.returncode, 10)
            assignment = parse_model(model, 1)
            self.assertTrue(assignment[1])
            self.assertEqual(verify_full_cnf_model(cnf, assignment), 1)

    def test_unsat_proof_checks_against_original_unmodified_cnf(self):
        with tempfile.TemporaryDirectory() as folder:
            result, cnf, proof, model = self.run_formula(folder, 'p cnf 1 2\n1 0\n-1 0\n', '1 0\n')
            self.assertEqual(result.returncode, 20)
            checked = subprocess.run([str(DEFAULT_CHECKER), str(cnf), str(proof)],
                                     capture_output=True, timeout=10)
            self.assertEqual(checked.returncode, 0)
            self.assertIn(b's VERIFIED', checked.stdout)

    def test_duplicate_or_trailing_hints_rejected(self):
        for phases in ('1 -1 0\n', '1 0 1\n', '1\n'):
            with self.subTest(phases=phases), tempfile.TemporaryDirectory() as folder:
                result, *_ = self.run_formula(folder, 'p cnf 1 1\n1 0\n', phases)
                self.assertNotEqual(result.returncode, 10)


if __name__ == '__main__':
    unittest.main()
