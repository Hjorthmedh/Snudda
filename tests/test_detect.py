import os
import sys
import unittest
import numpy as np

from snudda.detect.detect import SnuddaDetect
from snudda.place.create_cube_mesh import create_cube_mesh
from snudda.place.place import SnuddaPlace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


class TestDetect(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_detect")

        create_cube_mesh(file_name=os.path.join(network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        config_file = os.path.join(network_path, "network-config.json")
        position_file = os.path.join(network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(network_path, "voxels", "network-putative-synapses.hdf5")

        #  TODO: If d_view is None code run sin serial, add test parallel
        sp = SnuddaPlace(config_file=config_file, d_view=None, verbose=True)
        sp.place()

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file,
                               save_file=save_file, rc=None,
                               hyper_voxel_size=130, verbose=True)

    def test_detect(self):

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
                                     [100, 100, -40],  # To get a gap junction
                                     ]) * 1e-6

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

        ang = -np.pi / 2
        R_gj = np.array([[np.cos(ang), 0, np.sin(ang)],
                         [0, 1, 0],
                         [-np.sin(ang), 0, np.cos(ang)]])

        for idx in range(0, 10):  # Post synaptic neurons
            self.sd.neurons[idx]["rotation"] = R_x

        for idx in range(10, 20):  # Presynaptic neurons
            self.sd.neurons[idx]["rotation"] = R_y

        self.sd.neurons[20]["rotation"] = R_gj

        self.sd.detect(restart_detection_flag=True)

        synapse_voxel_loc = self.sd.hyper_voxel_synapses[:self.sd.hyper_voxel_synapse_ctr, 2:5]
        synapse_coords = synapse_voxel_loc * self.sd.voxel_size + self.sd.hyper_voxel_origo

        fig_path = os.path.join("networks", "network_testing_detect", "figures")

        if not os.path.exists(fig_path):
            os.mkdir(fig_path)

        if False:  # Set to True to include plot
            self.sd.plot_hyper_voxel(plot_neurons=True, fig_file_name="touch-detection-validation")

        # TODO: Also add tests for soma-axon synapses

        # Check location and source, targets.
        with self.subTest(stage="gap_junction_check"):
            self.assertEqual(self.sd.hyper_voxel_gap_junction_ctr, 1)
            self.assertTrue((self.sd.hyper_voxel_gap_junctions[0, 0:2] == [4, 20]).all())
            self.assertTrue(self.compare_voxels_to_coordinates(self.sd.hyper_voxel_gap_junctions[0, 6:9],
                                                               np.array([100, 100, 0]) * 1e-6))
        with self.subTest(stage="synapses_check"):
            print(f"synapse ctr {self.sd.hyper_voxel_synapse_ctr}")

            self.assertEqual(self.sd.hyper_voxel_synapse_ctr, 101)

            for post_id in range(0, 10):
                for pre_id in range(10, 20):
                    self.assertEqual(self.check_neuron_pair_has_synapse(pre_id, post_id), 1)

            for post_id in range(0, 10):
                for pre_id in range(0, 10):
                    self.assertFalse(self.check_neuron_pair_has_synapse(pre_id, post_id))

            for post_id in range(10, 20):
                for pre_id in range(10, 20):
                    self.assertFalse(self.check_neuron_pair_has_synapse(pre_id, post_id))

        with self.subTest(stage="synapse_sorting_check"):
            syn = self.sd.hyper_voxel_synapses[:self.sd.hyper_voxel_synapse_ctr, :2]
            syn_order = syn[:, 1] * len(self.sd.neurons) + syn[:, 0]
            self.assertTrue((np.diff(syn_order) >= 0).all())

        with self.subTest(stage="gap_junction_sorting_check"):
            gj = self.sd.hyper_voxel_gap_junctions[:self.sd.hyper_voxel_gap_junction_ctr, :2]
            gj_order = gj[:, 1] * len(self.sd.neurons) + gj[:, 0]
            self.assertTrue((np.diff(gj_order) >= 0).all())

        with self.subTest(stage="synapse_conductance_check"):
            cond = self.sd.hyper_voxel_synapses[:self.sd.hyper_voxel_synapse_ctr, 11] * 1e-12
            self.assertTrue(0.8e-9 < np.mean(cond) < 1.4e-9)  # mean 1.1e-9
            self.assertTrue(0.5e-9 < np.std(cond) < 2.5e-9)  # std 1.5e-9

        with self.subTest(stage="gap_junction_conductance_check"):
            cond = self.sd.hyper_voxel_gap_junctions[:self.sd.hyper_voxel_gap_junction_ctr, 10] * 1e-12
            self.assertTrue(2e-10 < np.mean(cond) < 8e-10)  # mean 5e-10

            # Only one Gap junction, cant check std here
            # self.assertTrue(0.1e-10 < np.std(cond) < 2e-9)  # std 1e-10

        # We should probably store the matrix as unsigned.
        self.assertTrue((self.sd.hyper_voxel_synapses >= 0).all())
        self.assertTrue((self.sd.hyper_voxel_gap_junctions >= 0).all())

        with self.subTest(stage="resiz_matrix_check"):
            old = self.sd.hyper_voxel_synapses.copy()
            self.sd.resize_hyper_voxel_synapses_matrix()
            # Check all synapses are preserved when resizing
            self.assertTrue((self.sd.hyper_voxel_synapses[:old.shape[0], :] == old).all())
            # Check new rows are empty
            self.assertTrue((self.sd.hyper_voxel_synapses[old.shape[0]:, :] == 0).all())

        # These test drawing not essential to Snudda, quite slow.
        if False:
            with self.subTest(stage="export_voxel_vis"):
                self.sd.export_voxel_visualisation_csv(neuron_id=np.arange(0, 10))

            with self.subTest(stage="plot_hyper_voxel"):
                # Matplotlib is kind of slow
                self.sd.plot_neurons_in_hyper_voxel(neuron_id=np.arange(0, 10),
                                                    neuron_colour=np.zeros((10, 3)),
                                                    show_plot=False, dpi=90)

            with self.subTest(stage="example-draw"):
                # Just checking that the drawing works
                self.sd.test_voxel_draw()

        print("Checking detect done.")

    def check_neuron_pair_has_synapse(self, pre_neuron, post_neuron):

        connections = dict()

        for synapse_row in self.sd.hyper_voxel_synapses[0:self.sd.hyper_voxel_synapse_ctr, :]:

            loc = (synapse_row[0], synapse_row[1])
            if loc in connections:
                connections[loc] += 1
            else:
                connections[loc] = 1

        if (pre_neuron, post_neuron) in connections:
            # print(f"pre: {pre_neuron}, post: {post_neuron}, connections: {connections[(pre_neuron, post_neuron)]}")
            return connections[(pre_neuron, post_neuron)]
        else:
            # print(f"No connection between pre {pre_neuron} and post {post_neuron}")
            return False

    def convert_to_coordinates(self, voxel_xyz):
        return np.array(voxel_xyz) * self.sd.voxel_size + self.sd.hyper_voxel_origo

    def compare_voxels_to_coordinates(self, voxel_index, coordinates):
        return (np.abs(self.convert_to_coordinates(voxel_index) - np.array(coordinates)) < self.sd.voxel_size).all()

    def test_detect_lines(self):

        # Cases to test
        # id 0-9 connecting to id 10-19, connection angle is 0-45 degrees
        # id 20 is 4 micrometer and parallel to another dendrite, no intersection

        neuron_positions = np.array([[0, 0, 0],  # Postsynaptiska
                                     [10, 10, 15],
                                     [20, 20, 30],
                                     [30, 30, 45],
                                     [40, 40, 60],
                                     [50, 50, 75],
                                     [60, 60, 90],
                                     [70, 70, 105],
                                     [80, 80, 120],
                                     [90, 90, 135],
                                     [50, -100, 0],  # Presynaptiska
                                     [60, -90, 15],
                                     [70, -80, 30],
                                     [80, -70, 45],
                                     [90, -60, 60],
                                     [100, -50, 75],
                                     [110, -40, 90],
                                     [120, -30, 105],
                                     [130, -20, 120],
                                     [140, -10, 135],
                                     [230, 0, 4],  # 4 micrometers from first neuron
                                     ]) * 1e-6

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

        for idx, ang in zip(range(10, 20), np.linspace(0, -np.pi / 4, 10)):  # Presynaptic neurons

            R_z = np.array([[np.cos(ang), -np.sin(ang), 0],
                            [np.sin(ang), np.cos(ang), 0],
                            [0, 0, 1]])

            print(f"idx = {idx}, ang = {ang}")

            self.sd.neurons[idx]["rotation"] = np.matmul(R_z, R_y)

        ang = np.pi / 2
        R_z = np.array([[np.cos(ang), -np.sin(ang), 0],
                        [np.sin(ang), np.cos(ang), 0],
                        [0, 0, 1]])

        self.sd.neurons[20]["rotation"] = np.matmul(R_z, R_y)

        self.sd.detect(restart_detection_flag=True)

        synapse_voxel_loc = self.sd.hyper_voxel_synapses[:self.sd.hyper_voxel_synapse_ctr, 2:5]
        synapse_coords = synapse_voxel_loc * self.sd.voxel_size + self.sd.hyper_voxel_origo

        if False:  # Set to True to include plot
            self.sd.plot_hyper_voxel(plot_neurons=True, fig_file_name="axon_dend_intersection_angle_0_45")

        with self.subTest(stage="synapses_check"):
            print(f"synapse ctr {self.sd.hyper_voxel_synapse_ctr}")

            # Depending on angle and location of voxels,
            # we can occasionally get more than one synapse at an intersection
            self.assertTrue(self.sd.hyper_voxel_synapse_ctr >= 10)

            # Verify that all pairs are connected
            for post_id in range(0, 10):
                pre_id = post_id + 10
                self.assertTrue(1 <= self.check_neuron_pair_has_synapse(pre_id, post_id) <= 2)

    # TODO: Gör ett test som distriburerar alla sections till sina respektive hypervoxlar,
    #       sen tvinga en touch detection som ignorerar den infon och tar med allt
    #       och se om några section id finns med som inte finns med i listorna...


if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_detect
