import json
from pathlib import Path
import tempfile
import unittest

from hn.color_filter import check_functional
from tools.verify_tower import VerificationError


def triangle():
    points = [{'parts': [['0']*12 for _ in range(4)]} for _ in range(3)]
    points[1]['parts'][0][0] = '1'
    points[2]['parts'][0][7] = '1'
    return {'schema': 'hn-cyclotomic-tower-v1',
            'field': {'cyclotomic_order': 42, 'radicands': [5, 33]},
            'vertices': points, 'edges': [[0, 1], [0, 2], [1, 2]]}


class ColorFilterTests(unittest.TestCase):
    def check(self, graph, functional):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'graph.json'
            path.write_text(json.dumps(graph))
            return check_functional(path, functional, Path(folder)/'out')

    def test_full_exact_coloring_passes_but_never_target(self):
        functional = [1]+[0]*47; functional[7] = 2
        result = self.check(triangle(), functional)
        self.assertEqual(result['status'], 'FIVE_COLORING_VERIFIED')
        self.assertFalse(result['goal_satisfied'])
        self.assertEqual(result['geometry']['edges'], 3)

    def test_new_monochromatic_edge_rejects_rule(self):
        result = self.check(triangle(), [1]+[0]*47)
        self.assertEqual(result['status'], 'FUNCTIONAL_REJECTED')
        self.assertEqual(result['first_monochromatic_edges'], [[0, 2]])

    def test_nonunit_colored_graph_never_prunes(self):
        graph = triangle(); graph['vertices'][2]['parts'][0][7] = '2'
        functional = [1]+[0]*47; functional[7] = 2
        with self.assertRaises(VerificationError):
            self.check(graph, functional)

    def test_pole_is_inapplicable_not_obstruction(self):
        graph = triangle(); graph['vertices'][2]['parts'][0][7] = '1/5'
        result = self.check(graph, [1]+[0]*47)
        self.assertEqual(result['status'], 'NOT_APPLICABLE')
        self.assertFalse(result['goal_satisfied'])


if __name__ == '__main__':
    unittest.main()
