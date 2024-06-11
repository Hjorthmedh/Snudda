import unittest
import os

import numpy as np
from snudda.analyse.analyse_topology_activity import SnuddaAnalyseTopologyActivity


class TestTopologyAnalysis(unittest.TestCase):

    def setUp(self) -> None:
        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

    def test_basics(self):
        sa = SnuddaAnalyseTopologyActivity()

        with self.subTest("match_closest_unique"):

            spikes_a = np.array([1, 2, 3, 4, 5])
            spikes_b = np.array([1.001, 1.002, 1.995, 2.5, 4.020, 5])
            delta_t = 10e-3

            delta_a, delta_b = sa.match_closest_unique(spikes_a, spikes_b, delta_t=delta_t)

            idx_a = np.invert(np.isnan(delta_a))
            idx_b = np.invert(np.isnan(delta_b))

            ref_a = np.array([0.001, -0.005, np.nan, np.nan, 0])
            ref_b = np.array([-0.001, np.nan, 0.005, np.nan, np.nan, 0])

            self.assertTrue(np.allclose(delta_a[idx_a], ref_a[idx_a]))
            self.assertTrue((np.isnan(ref_a[np.isnan(delta_a)])).all())

            self.assertTrue(np.allclose(delta_b[idx_b], ref_b[idx_b]))
            self.assertTrue((np.isnan(ref_b[np.isnan(delta_b)])).all())

        with self.subTest("triggered"):

            spikes_a = np.array([1, 2, 3, 4, 5, 6])
            spikes_b = np.array([0, 0.1, 1.1, 1.2, 3.2, 4, 5.5, 10, 11])

            delta_t = sa.get_spike_triggered_deltas(spikes_a, spikes_b)

            self.assertTrue(np.isnan(delta_t[0:2]).all())
            self.assertTrue(np.allclose(delta_t[2:], np.array([0.1, 0.2, 0.2, 0, 0.5, 4, 5])))


if __name__ == '__main__':
    unittest.main()
