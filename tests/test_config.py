import unittest

from pbmc_project.config import KMEANS_RANGE, LOUVAIN_RESOLUTIONS, MARKER_SETS, format_resolution


class ConfigTest(unittest.TestCase):
    def test_parameter_ranges_match_project_plan(self):
        self.assertEqual(KMEANS_RANGE[0], 4)
        self.assertEqual(KMEANS_RANGE[-1], 12)
        self.assertEqual(LOUVAIN_RESOLUTIONS, (0.4, 0.8, 1.2))

    def test_marker_sets_have_supporting_genes(self):
        self.assertIn("T cells", MARKER_SETS)
        self.assertIn("CD3D", MARKER_SETS["T cells"])
        self.assertTrue(all(genes for genes in MARKER_SETS.values()))

    def test_louvain_resolution_key_format(self):
        self.assertEqual(format_resolution(0.8), "0_8")


if __name__ == "__main__":
    unittest.main()
