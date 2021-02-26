import os
import unittest
import numpy as np

from snudda.neuron_morphology import NeuronMorphology


class MyTestCase(unittest.TestCase):

    def setUp(self):

        #swc_file = os.path.join(os.path.dirname(__file__), "validation", "ballanddoublestick", "double.swc")

        swc_file = os.path.join(os.path.dirname(__file__), "validation", "striatum", "fs",
                                "str-fs-e161205_FS1-mMTC180800A-IDB-v20190312", "MTC180800A-IDB-cor-rep.swc")
        self.nm = NeuronMorphology(swc_filename=swc_file, load_morphology=True)

    def test_input_location(self, stage="dist_to_soma"):

        synapse_density = "(d > 100e-6)*1"
        rng = np.random.default_rng(123456)  #

        xyz, sec_id, sec_x, dist_to_soma = self.nm.dendrite_input_locations(synapse_density=synapse_density,
                                                                            rng=rng, num_locations=100)

        # Please note that num_locations currently does not guarantee 100 synpases when requesting it

        # 3e-6 due to compartment length sampled at 3 micrometers
        self.assertTrue((dist_to_soma > 100e-6 - 3e-6).all())

        self.assertEqual(xyz.shape[1], 3)
        self.assertEqual(len(dist_to_soma), xyz.shape[0])
        self.assertEqual(len(dist_to_soma), len(sec_id))
        self.assertEqual(len(sec_id), len(sec_x))

        # Repeat test but for smaller than 100
        synapse_density = "(d < 200e-6)*1"

        xyz, sec_id, sec_x, dist_to_soma = self.nm.dendrite_input_locations(synapse_density=synapse_density,
                                                                            rng=rng, num_locations=100)
        # 3e-6 due to compartment length sampled at 3 micrometers
        self.assertTrue((dist_to_soma < 200e-6 + 3e-6).all())

    def test_rand_rotation(self, stage="rand_rotation"):

        for idx in range(0, 100):
            rot_mat = self.nm.rand_rotation_matrix()
            self.assertAlmostEqual(np.linalg.det(rot_mat), 1, places=10)

    def test_clone(self, stage="clone"):
        new_nm = self.nm.clone()
        
        self.assertTrue((self.nm.dend == new_nm.dend).all())
        self.assertTrue((self.nm.axon == new_nm.axon).all())
        self.assertTrue((self.nm.soma == new_nm.soma).all())

        # Make sure that clone is a copy, and don't point back to same arrays
        ang = np.pi
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])

        new_nm.place(rotation=R_x, position=np.array([0, 0, 0]))

        self.assertTrue((np.abs(self.nm.dend[:, 0] - new_nm.dend[:, 0]) < 1e-6).all())
        self.assertTrue((np.abs(self.nm.dend[:, 1] + new_nm.dend[:, 1]) < 1e-6).all())
        self.assertTrue((np.abs(self.nm.dend[:, 2] + new_nm.dend[:, 2]) < 1e-6).all())

        self.assertTrue((np.abs(self.nm.axon[:, 0] - new_nm.axon[:, 0]) < 1e-6).all())
        self.assertTrue((np.abs(self.nm.axon[:, 1] + new_nm.axon[:, 1]) < 1e-6).all())
        self.assertTrue((np.abs(self.nm.axon[:, 2] + new_nm.axon[:, 2]) < 1e-6).all())

        new_nm.place(position=np.array([1, 2, 3]))
        self.assertTrue((np.abs(new_nm.soma[0,:3] - np.array([1, 2, 3])) < 1e-6).all())


if __name__ == '__main__':
    unittest.main()
