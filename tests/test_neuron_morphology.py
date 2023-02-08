import os
import unittest
import numpy as np
import scipy

from snudda.neurons.neuron_morphology_extended import NeuronMorphologyExtended


class NeuronMorphologyExtendedTestCase(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        #swc_file = os.path.join(os.path.dirname(__file__), "validation", "ballanddoublestick", "double.swc")

        swc_file = os.path.join(os.path.dirname(__file__), "validation", "striatum", "fs",
                                "str-fs-e161205_FS1-mMTC180800A-IDB-v20190312", "MTC180800A-IDB-cor-rep.swc")
        self.nm = NeuronMorphologyExtended(swc_filename=swc_file, load_morphology=True)  #, use_cache=False

    def test_input_location(self, stage="dist_to_soma"):

        synapse_density = "(d > 100e-6)*1"
        rng = np.random.default_rng(123456)  #

        xyz, sec_id, sec_x, dist_to_soma = self.nm.dendrite_input_locations(synapse_density_str=synapse_density,
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

        xyz, sec_id, sec_x, dist_to_soma = self.nm.dendrite_input_locations(synapse_density_str=synapse_density,
                                                                            rng=rng, num_locations=100)
        # 3e-6 due to compartment length sampled at 3 micrometers
        self.assertTrue((dist_to_soma < 200e-6 + 3e-6).all())

#   -- rand_rotation is now moved to rotation.py
#
#     def test_rand_rotation(self, stage="rand_rotation"):
#
#         for idx in range(0, 100):
#             rot_mat = self.nm.rand_rotation_matrix()
#             self.assertAlmostEqual(np.linalg.det(rot_mat), 1, places=10)

    def test_clone(self, stage="clone"):
        new_nm = self.nm.clone()

        self.assertTrue((self.nm.morphology_data["neuron"].geometry == new_nm.morphology_data["neuron"].geometry).all())
        self.assertTrue((self.nm.morphology_data["neuron"].section_data == new_nm.morphology_data["neuron"].section_data).all())

        # Make sure that clone is a copy, and don't point back to same arrays
        ang = np.pi
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])

        new_nm.place(rotation=R_x, position=np.array([0, 0, 0]))

        self.assertTrue((np.abs(self.nm.morphology_data["neuron"].geometry[:, 0] - new_nm.morphology_data["neuron"].geometry[:, 0]) < 1e-6).all())
        self.assertTrue((np.abs(self.nm.morphology_data["neuron"].geometry[:, 1] + new_nm.morphology_data["neuron"].geometry[:, 1]) < 1e-6).all())
        self.assertTrue((np.abs(self.nm.morphology_data["neuron"].geometry[:, 2] + new_nm.morphology_data["neuron"].geometry[:, 2]) < 1e-6).all())

        new_nm2 = self.nm.clone()
        new_nm2.place(position=np.array([1, 2, 3]))
        self.assertTrue((np.abs(new_nm2.position - np.array([1, 2, 3])) < 1e-6).all())

    def test_cluster_synapses(self, stage="cluster_synapses"):

        cluster_spread = 30e-6
        n_synapses = 10

        cluster_sec_x, syn_coords, soma_dist = self.nm.cluster_synapses(sec_id=29, sec_x=0.5,
                                                                        count=n_synapses, distance=cluster_spread,
                                                                        rng=np.random.default_rng(20220125))

        self.assertEqual(len(cluster_sec_x), n_synapses)
        self.assertTrue(np.max(soma_dist) - np.min(soma_dist) <= cluster_spread)

        d = scipy.spatial.distance.pdist(syn_coords)
        self.assertTrue(np.all(d <= cluster_spread))
        self.assertTrue(np.any(d >= 10e-6))

        # TODO: Verify self.nm.sec_id_to_len[0] value ... ie segment length of soma?
        # xx = self.nm.compartment_length()
        # np.sum(self.nm.sec_id_to_len[1:])
        # ska vara samma som np.sum(xx)
        # och self.nm.sec_id_to_len[0] deals with the soma.

        if False:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax = self.nm.plot_neuron(plot_axon=False, plot_dendrite=True, axis=ax)
            ax.scatter(syn_coords[:, 0], syn_coords[:, 1], syn_coords[:, 2], s=5, c="red")
            fig.show()
            plt.pause(3)

            import pdb
            pdb.set_trace()


if __name__ == '__main__':
    unittest.main()
