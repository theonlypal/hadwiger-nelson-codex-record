import json
from pathlib import Path
import tempfile
import unittest

from tools.verify_kempe_witnesses import replay
from tools.verify import _sha256


class CompactColoringWitnessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.graph = self.root/"graph.json"
        self.write(self.graph, {"schema": "hn-cyclotomic-tower-v1",
                   "vertices": [{}, {}, {}], "edges": [[0, 1], [1, 2]]})
        self.sha = _sha256(self.graph)
        self.model = self.root/"model-000"; self.model.mkdir()
        self.write(self.model/"baseline-coloring.json", {"k": 5,
                   "colors": [0, 1, 0], "graph_sha256": self.sha})
        self.swap = {"color_pair": [0, 2], "swapped_vertices": [0],
                     "affected_edge_incidences_checked": 1,
                     "separated_pair_indices": [0]}
        self.swaps = {"baseline": "baseline-coloring.json", "graph_sha256": self.sha,
                      "witnesses": [self.swap]}
        self.write(self.model/"swap-witnesses.json", self.swaps)
        self.write(self.root/"refutations.json", {"pairs": [[0, 2]],
                   "refutations_by_pair_index": {"0": {
                   "kind": "checked-compact-Kempe-coloring",
                   "witness": "model-000/swap-witnesses.json", "swap_index": 0}}})
        self.write(self.root/"unresolved-pairs.json", [])

    def write(self, path, value):
        path.write_text(json.dumps(value))

    def test_valid_compact_counterexample(self):
        result = replay(self.graph, self.root)
        self.assertEqual(result["pairs_disproved_forced_equal"], 1)
        self.assertFalse(result["goal_satisfied"])

    def test_monochromatic_changed_boundary_rejected(self):
        self.swap.update({"color_pair": [0, 1], "swapped_vertices": [0, 1],
                          "affected_edge_incidences_checked": 3})
        self.write(self.model/"swap-witnesses.json", self.swaps)
        with self.assertRaisesRegex(ValueError, "monochromatic"):
            replay(self.graph, self.root)

    def test_swap_that_does_not_separate_rejected(self):
        self.swap["swapped_vertices"] = [0, 2]
        self.swap["affected_edge_incidences_checked"] = 2
        self.write(self.model/"swap-witnesses.json", self.swaps)
        with self.assertRaisesRegex(ValueError, "does not separate"):
            replay(self.graph, self.root)

    def test_false_unresolved_partition_rejected(self):
        self.write(self.root/"unresolved-pairs.json", [[0, 2]])
        with self.assertRaisesRegex(ValueError, "exact complement"):
            replay(self.graph, self.root)


if __name__ == "__main__":
    unittest.main()
