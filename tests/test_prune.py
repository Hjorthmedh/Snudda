import os
import unittest

from snudda.place.create_cube_mesh import create_cube_mesh
from snudda.detect.detect import SnuddaDetect
from snudda.utils.load import SnuddaLoad
from snudda.neurons.neuron_morphology_extended import NeuronMorphologyExtended
from snudda.place.place import SnuddaPlace
import numpy as np

from snudda.detect.prune import SnuddaPrune


class TestPrune(unittest.TestCase):
    
    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        self.network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_prune3")

        create_cube_mesh(file_name=os.path.join(self.network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        config_file = os.path.join(self.network_path, "network-config.json")
        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        sp = SnuddaPlace(config_file=config_file, d_view=None, verbose=True)

        sp.parse_config()
        sp.write_data(position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        from snudda.utils.reposition_neurons import RepositionNeurons
        repos = RepositionNeurons(position_file=position_file)

        # Reposition the neurons so that we know how many synapses and where they will be located before pruning
        # OBS, these positions and rotations are not written to the HDF5 file, they are only in memory!
        neuron_positions = np.array([[0, 20, 0],  # Postsynaptiska
                                     [0, 40, 0],
                                     [0, 60, 0],
                                     [0, 80, 0],
                                     [0, 100, 0],
                                     [0, 120, 0],
                                     [0, 140, 0],
                                     [0, 160, 0],
                                     [0, 180, 0],
                                     [0, 200, 0],
                                     [20, 0, 0],  # Presynaptiska
                                     [40, 0, 0],
                                     [60, 0, 0],
                                     [80, 0, 0],
                                     [100, 0, 0],
                                     [120, 0, 0],
                                     [140, 0, 0],
                                     [160, 0, 0],
                                     [180, 0, 0],
                                     [200, 0, 0],
                                     [70, 0, 500],  # For gap junction check
                                     [110, 0, 500],
                                     [150, 0, 500],
                                     [190, 0, 500],
                                     [0, 70, 500],
                                     [0, 110, 500],
                                     [0, 150, 500],
                                     [0, 190, 500],
                                     ]) * 1e-6

        # TODO: Add potential for gap junctions also by having 5 + 5 neurons in other grid
        
        ang = -np.pi / 2
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])
        
        ang = np.pi / 2
        R_y = np.array([[np.cos(ang), 0, np.sin(ang)],
                        [0, 1, 0],
                        [-np.sin(ang), 0, np.cos(ang)]])

        for idx in range(0, 10):  # Post synaptic neurons
            repos.place(neuron_id=idx, position=neuron_positions[idx, :], rotation=R_x)

        for idx in range(10, 20):  # Presynaptic neurons
            repos.place(neuron_id=idx, position=neuron_positions[idx, :], rotation=R_y)

        for idx in range(24, 28):  # GJ neurons
            repos.place(neuron_id=idx, position=neuron_positions[idx, :], rotation=R_x)

        ang = np.pi / 2
        R_z = np.array([[np.cos(ang), -np.sin(ang), 0],
                        [np.sin(ang), np.cos(ang), 0],
                        [0, 0, 1]])

        for idx in range(20, 24):  # GJ neurons
            repos.place(neuron_id=idx, position=neuron_positions[idx, :], rotation=np.matmul(R_z, R_x))

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file,
                               save_file=save_file, rc=None,
                               hyper_voxel_size=120, verbose=True)

        self.sd.detect(restart_detection_flag=True)

        if False:
            self.sd.process_hyper_voxel(1)
            self.sd.plot_hyper_voxel(plot_neurons=True,
                                     fig_file_name=os.path.join(self.network_path,
                                                                "test_prune_figure.png"))

            import pdb
            pdb.set_trace()

    def test_prune(self):

        pruned_output = os.path.join(self.network_path, "network-synapses.hdf5")

        with self.subTest(stage="No-pruning"):

            sp = SnuddaPrune(network_path=self.network_path, config_file=None, verbose=True, keep_files=True)  # Use default config file
            sp.prune()
            sp = []

            # Load the pruned data and check it

            sl = SnuddaLoad(pruned_output)
            # TODO: Call a plot function to plot entire network with synapses and all

            self.assertEqual(sl.data["num_synapses"], (20*8 + 10*2)*2)  # Update, now AMPA+GABA, hence *2 at end

            # This checks that all synapses are in order
            # The synapse sort order is destID, sourceID, synapsetype (channel model id).

            syn = sl.data["synapses"][:sl.data["num_synapses"], :]
            syn_order = (syn[:, 1] * len(self.sd.neurons) + syn[:, 0]) * 12 + syn[:, 6]  # The 12 is maxChannelModelID
            self.assertTrue((np.diff(syn_order) >= 0).all())

            # Note that channel model id is dynamically allocated, starting from 10 (GJ have ID 3)
            # Check that correct number of each type
            self.assertEqual(np.sum(sl.data["synapses"][:, 6] == 10), 20*8 + 10*2)
            self.assertEqual(np.sum(sl.data["synapses"][:, 6] == 11), 20*8 + 10*2)

            self.assertEqual(sl.data["num_gap_junctions"], 4*4*4)
            gj = sl.data["gap_junctions"][:sl.data["num_gap_junctions"], :2]
            gj_order = gj[:, 1] * len(self.sd.neurons) + gj[:, 0]
            self.assertTrue((np.diff(gj_order) >= 0).all())

        with self.subTest(stage="load-testing"):
            sl = SnuddaLoad(pruned_output, verbose=True)

            # Try and load a neuron
            n = sl.load_neuron(neuron_id=0)
            self.assertTrue(type(n) == NeuronMorphologyExtended)

            syn_ctr = 0
            for s in sl.synapse_iterator(chunk_size=50):
                syn_ctr += s.shape[0]
            self.assertEqual(syn_ctr, sl.data["num_synapses"])
        
            gj_ctr = 0
            for gj in sl.gap_junction_iterator(chunk_size=50):
                gj_ctr += gj.shape[0]
            self.assertEqual(gj_ctr, sl.data["num_gap_junctions"])

            syn, syn_coords = sl.find_synapses(pre_id=14)
            self.assertTrue((syn[:, 0] == 14).all())
            self.assertEqual(syn.shape[0], 40)

            syn, syn_coords = sl.find_synapses(post_id=3)
            self.assertTrue((syn[:, 1] == 3).all())
            self.assertEqual(syn.shape[0], 36)

            cell_id_perm = sl.get_neuron_id_of_type("ballanddoublestick", random_permute=True, num_neurons=28)
            cell_id = sl.get_neuron_id_of_type("ballanddoublestick", random_permute=False)

            self.assertEqual(len(cell_id_perm), 28)
            self.assertEqual(len(cell_id), 28)
            
            for cid in cell_id_perm:
                self.assertTrue(cid in cell_id)

        # It is important merge file has synapses sorted with dest_id, source_id as sort order since during pruning
        # we assume this to be able to quickly find all synapses on post synaptic cell.
        # TODO: Also include the ChannelModelID in sorting check
        with self.subTest("Checking-merge-file-sorted"):

            for mf in ["temp/synapses-for-neurons-0-to-28-MERGE-ME.hdf5",
                       "temp/gap_junctions-for-neurons-0-to-28-MERGE-ME.hdf5",
                       "network-synapses.hdf5"]:

                merge_file = os.path.join(self.network_path, mf)

                sl = SnuddaLoad(merge_file, verbose=True)
                if "synapses" in sl.data:
                    syn = sl.data["synapses"][:sl.data["num_synapses"], :2]
                    syn_order = syn[:, 1] * len(self.sd.neurons) + syn[:, 0]
                    self.assertTrue((np.diff(syn_order) >= 0).all())

                if "gap_junctions" in sl.data:
                    gj = sl.data["gap_junctions"][:sl.data["num_gap_junctions"], :2]
                    gj_order = gj[:, 1] * len(self.sd.neurons) + gj[:, 0]
                    self.assertTrue((np.diff(gj_order) >= 0).all())

        with self.subTest("synapse-f1"):
            # Test of f1
            testing_config_file = os.path.join(self.network_path, "network-config-test-1.json")
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it

            sl = SnuddaLoad(pruned_output, verbose=True)
            # Setting f1=0.5 in config should remove 50% of GABA synapses, but does so randomly, for AMPA we used f1=0.9
            gaba_id = sl.data["connectivity_distributions"]["ballanddoublestick","ballanddoublestick"]["GABA"]["channel_model_id"]
            ampa_id = sl.data["connectivity_distributions"]["ballanddoublestick","ballanddoublestick"]["AMPA"]["channel_model_id"]

            n_gaba = np.sum(sl.data["synapses"][:, 6] == gaba_id)
            n_ampa = np.sum(sl.data["synapses"][:, 6] == ampa_id)

            # gaba: var(x) = n * p * (1-p), ie 180 * 0.5 * 0. 5 = 45 --> std(x) = 6.7, 2* std(x) = 13
            # ampa: 2*std(x) = 2*sqrt(180*0.9*0.1) = 8
            self.assertTrue((20*8 + 10*2)*0.5 - 13 < n_gaba < (20*8 + 10*2)*0.5 + 13)
            self.assertTrue((20*8 + 10*2)*0.9 - 8 < n_ampa < (20*8 + 10*2)*0.9 + 8)

        with self.subTest("synapse-softmax"):
            # Test of softmax
            testing_config_file = os.path.join(self.network_path, "network-config-test-2.json")  # Only GABA synapses in this config
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it
            sl = SnuddaLoad(pruned_output)
            # Softmax reduces number of synapses
            self.assertTrue(sl.data["num_synapses"] < 20*8 + 10*2)

        with self.subTest("synapse-mu2"):
            # Test of mu2
            testing_config_file = os.path.join(self.network_path, "network-config-test-3.json")
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it
            sl = SnuddaLoad(pruned_output)
            # With mu2 having 2 synapses means 50% chance to keep them, having 1 will be likely to have it removed
            self.assertTrue(20*8*0.5 - 13 < sl.data["num_synapses"] < 20*8*0.5 + 13)

        with self.subTest("synapse-a3"):
            # Test of a3
            testing_config_file = os.path.join(self.network_path, "network-config-test-4.json")
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it
            sl = SnuddaLoad(pruned_output)

            # a3=0.6 means 40% chance to remove all synapses between a pair
            self.assertTrue((20*8 + 10*2)*0.6 - 14 < sl.data["num_synapses"] < (20*8 + 10*2)*0.6 + 14)

        with self.subTest("synapse-distance-dependent-pruning"):
            # Testing distance dependent pruning
            testing_config_file = os.path.join(self.network_path, "network-config-test-5.json")
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it
            sl = SnuddaLoad(pruned_output)

            # "1*(d >= 100e-6)" means we remove all synapses closer than 100 micrometers
            print(f"num synapses : {sl.data['num_synapses']}")

            self.assertEqual(sl.data["num_synapses"], 20*5)
            self.assertTrue((sl.data["synapses"][:, 8] >= 100).all())  # Column 8 -- distance to soma in micrometers

        # TODO: Need to do same test for Gap Junctions also -- but should be same results, since same codebase
        with self.subTest("gap-junction-f1"):
            # Test of f1
            testing_config_file = os.path.join(self.network_path, "network-config-test-6.json")
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it

            sl = SnuddaLoad(pruned_output)
            # Setting f1=0.7 in config should remove 30% of gap junctions, but does so randomly
            self.assertTrue(64*0.7 - 10 < sl.data["num_gap_junctions"] < 64*0.7 + 10)

        with self.subTest("gap-junction-softmax"):
            # Test of softmax
            testing_config_file = os.path.join(self.network_path, "network-config-test-7.json")
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it
            sl = SnuddaLoad(pruned_output)
            # Softmax reduces number of synapses
            self.assertTrue(sl.data["num_gap_junctions"] < 16*2 + 10)

        with self.subTest("gap-junction-mu2"):
            # Test of mu2
            testing_config_file = os.path.join(self.network_path, "network-config-test-8.json")
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it
            sl = SnuddaLoad(pruned_output)
            # With mu2 having 4 synapses means 50% chance to keep them, having 1 will be likely to have it removed
            self.assertTrue(64*0.5 - 10 < sl.data["num_gap_junctions"] < 64*0.5 + 10)

        with self.subTest("gap-junction-a3"):
            # Test of a3
            testing_config_file = os.path.join(self.network_path, "network-config-test-9.json")
            sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
            sp.prune()

            # Load the pruned data and check it
            sl = SnuddaLoad(pruned_output, verbose=True)

            # a3=0.7 means 30% chance to remove all synapses between a pair
            self.assertTrue(64*0.7 - 10 < sl.data["num_gap_junctions"] < 64*0.7 + 10)

        if False:  # Distance dependent pruning currently not implemented for gap junctions
            with self.subTest("gap-junction-distance-dependent-pruning"):
                # Testing distance dependent pruning
                testing_config_file = os.path.join(self.network_path, "network-config-test-10.json")
                sp = SnuddaPrune(network_path=self.network_path, config_file=testing_config_file, verbose=True, keep_files=True)  # Use default config file
                sp.prune()

                # Load the pruned data and check it
                sl = SnuddaLoad(pruned_output, verbose=True)

                # "1*(d <= 120e-6)" means we remove all synapses further away than 100 micrometers
                self.assertEqual(sl.data["num_gap_junctions"], 2*4*4)
                self.assertTrue((sl.data["gap_junctions"][:, 8] <= 120).all())  # Column 8 -- distance to soma in micrometers


if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_prune
