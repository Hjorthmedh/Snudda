import os
import unittest

from snudda.create_cube_mesh import create_cube_mesh
from snudda.detect import SnuddaDetect
from snudda.load import SnuddaLoad
from snudda.place import SnuddaPlace
import numpy as np

from snudda.prune import SnuddaPrune


class TestPrune(unittest.TestCase):
    
    def setUp(self):

        os.chdir(os.path.dirname(__file__))
        self.network_path = os.path.join(os.path.dirname(__file__), "tests", "network_testing_prune2")

        create_cube_mesh(file_name=os.path.join(self.network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        config_file = os.path.join(self.network_path, "network-config.json")
        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        sp = SnuddaPlace(config_file=config_file, d_view=None)

        sp.read_config()
        sp.write_data(position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file,
                               save_file=save_file, rc=None,
                               hyper_voxel_size=120)

        # Reposition the neurons so we know how many synapses and where they will be located before pruning
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

        for idx, pos in enumerate(neuron_positions):
            self.sd.neurons[idx]["position"] = pos
        
        ang = -np.pi / 2
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])
        
        ang = np.pi / 2
        R_y = np.array([[np.cos(ang), 0, np.sin(ang)],
                        [0, 1, 0],
                        [-np.sin(ang), 0, np.cos(ang)]])

        for idx in range(0, 10):  # Post synaptic neurons
            self.sd.neurons[idx]["rotation"] = R_x
        
        for idx in range(10, 20):  # Presynaptic neurons
            self.sd.neurons[idx]["rotation"] = R_y

        for idx in range(24, 28):  # GJ neurons
            self.sd.neurons[idx]["rotation"] = R_x

        ang = np.pi / 2
        R_z = np.array([[np.cos(ang), -np.sin(ang), 0],
                        [np.sin(ang), np.cos(ang), 0],
                        [0, 0, 1]])

        for idx in range(20, 24):  # GJ neurons
            self.sd.neurons[idx]["rotation"] = np.matmul(R_z, R_x)

        self.sd.detect(restart_detection_flag=True)

        if False:
            self.sd.process_hyper_voxel(1)
            self.sd.plot_hyper_voxel(plot_neurons=True)
            # TODO: Check why soma is plotted in wrong place? Mistake with origo plotoffset?
            import pdb
            pdb.set_trace()

    def test_prune(self):

        work_log = os.path.join(self.network_path, "log", "network-detect-worklog.hdf5")
        pruned_output = os.path.join(self.network_path, "network-pruned-synapses.hdf5")

        with self.subTest(stage="No-pruning"):

            sp = SnuddaPrune(work_history_file=work_log, config_file=None)  # Use default config file
            sp.prune(pre_merge_only=False)
            sp = []

            # Load the pruned data and check it

            sl = SnuddaLoad(pruned_output)
            # TODO: Call a plot function to plot entire network with synapses and all

            self.assertEqual(sl.data["nSynapses"], 20*8 + 10*2)

            # This checks that all synapses are in order
            syn = sl.data["synapses"][:sl.data["nSynapses"], :2]
            syn_order = syn[:, 1] * len(self.sd.neurons) + syn[:, 0]
            self.assertTrue((np.diff(syn_order) >= 0).all())

            self.assertEqual(sl.data["nGapJunctions"], 4*4*4)
            gj = sl.data["gapJunctions"][:sl.data["nGapJunctions"], :2]
            gj_order = gj[:, 1] * len(self.sd.neurons) + gj[:, 0]
            self.assertTrue((np.diff(gj_order) >= 0).all())

        # It is important merge file has synapses sorted with dest_id, source_id as sort order since during pruning
        # we assume this to be able to quickly find all synapses on post synaptic cell.
        with self.subTest("Checking-merge-file-sorted"):
            merge_file = os.path.join(self.network_path, "network-putative-synapses-MERGED.hdf5")

            sl = SnuddaLoad(merge_file)
            syn = sl.data["synapses"][:sl.data["nSynapses"], :2]
            syn_order = syn[:, 1] * len(self.sd.neurons) + syn[:, 0]
            self.assertTrue((np.diff(syn_order) >= 0).all())

            gj = sl.data["gapJunctions"][:sl.data["nGapJunctions"], :2]
            gj_order = gj[:, 1] * len(self.sd.neurons) + gj[:, 0]
            self.assertTrue((np.diff(gj_order) >= 0).all())

        # Test of f1
        testing_config_file = os.path.join(self.network_path, "network-config-test-1.json")
        sp = SnuddaPrune(work_history_file=work_log, config_file=testing_config_file)  # Use default config file
        sp.prune(pre_merge_only=False)

        # Load the pruned data and check it

        sl = SnuddaLoad(pruned_output)
        # Setting f1=0.5 in config should remove 50% of synapses, but does so randomly
        self.assertTrue((20*8 + 10*2)*0.5 - 10 < sl.data["nSynapses"] < (20*8 + 10*2)*0.5 + 10)

        # Test of softmax
        testing_config_file = os.path.join(self.network_path, "network-config-test-2.json")
        sp = SnuddaPrune(work_history_file=work_log, config_file=testing_config_file)  # Use default config file
        sp.prune(pre_merge_only=False)

        # Load the pruned data and check it
        sl = SnuddaLoad(pruned_output)
        # Softmax reduces number of synapses
        self.assertTrue(sl.data["nSynapses"] < 20*8 + 10*2)

        # Test of mu2
        testing_config_file = os.path.join(self.network_path, "network-config-test-3.json")
        sp = SnuddaPrune(work_history_file=work_log, config_file=testing_config_file)  # Use default config file
        sp.prune(pre_merge_only=False)

        # Load the pruned data and check it
        sl = SnuddaLoad(pruned_output)
        # With mu2 having 2 synapses means 50% chance to keep them, having 1 will be likely to have it removed
        self.assertTrue(20*8*0.5 - 10 < sl.data["nSynapses"] < 20*8*0.5 + 10)

        # Test of a3
        testing_config_file = os.path.join(self.network_path, "network-config-test-4.json")
        sp = SnuddaPrune(work_history_file=work_log, config_file=testing_config_file)  # Use default config file
        sp.prune(pre_merge_only=False)

        # Load the pruned data and check it
        sl = SnuddaLoad(pruned_output)

        # a3=0.6 means 40% chance to remove all synapses between a pair
        self.assertTrue((20*8 + 10*2)*0.6 - 10 < sl.data["nSynapses"] < (20*8 + 10*2)*0.6 + 10)

        # Testing distance dependent pruning
        testing_config_file = os.path.join(self.network_path, "network-config-test-5.json")
        sp = SnuddaPrune(work_history_file=work_log, config_file=testing_config_file)  # Use default config file
        sp.prune(pre_merge_only=False)

        # Load the pruned data and check it
        sl = SnuddaLoad(pruned_output)

        # "1*(d >= 100e-6)" means we remove all synapses closer than 100 micrometers
        self.assertEqual(sl.data["nSynapses"], 20*6)
        self.assertTrue((sl.data["synapses"][:, 8] >= 100).all())  # Column 8 -- distance to soma in micrometers

        # TODO: Need to do same tests for Gap Junctions also -- but should be same results, since same codebase
        