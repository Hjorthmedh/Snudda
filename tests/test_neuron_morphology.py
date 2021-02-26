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
        rng = np.random.default_rng()  #

        xyz, sec_id, sec_x, dist_to_soma = self.nm.dendrite_input_locations(synapse_density=synapse_density,
                                                                            rng=rng, num_locations=100)

        # !!! We did not always get 100 synapses!! check why

        # 3e-6 due to compartment length sampled at 3 micrometers
        self.assertTrue((dist_to_soma > 100e-6 - 3e-6).all())
        self.assertEqual(len(dist_to_soma), 100)
        self.assertEqual(xyz.shape[0], 100)
        self.assertEqual(xyz.shape[1], 3)
        self.assertEqual(len(sec_id), 100)
        self.assertEqual(len(sec_x), 100)

        # Repeat test but for smaller than 100
        synapse_density = "(d < 100e-6)*1"
        rng = np.random.default_rng()  #

        xyz, sec_id, sec_x, dist_to_soma = self.nm.dendrite_input_locations(synapse_density=synapse_density,
                                                                            rng=rng, num_locations=100)

        self.assertTrue((dist_to_soma < 100e-6 + 3e-6).all())



if __name__ == '__main__':
    unittest.main()
