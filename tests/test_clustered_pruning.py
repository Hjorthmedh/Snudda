import os
import sys
import unittest

import numpy as np
import scipy

from snudda.detect.detect import SnuddaDetect
from snudda.detect.prune import SnuddaPrune
from snudda.place.create_cube_mesh import create_cube_mesh
from snudda.place.place import SnuddaPlace
from snudda.utils.reposition_neurons import RepositionNeurons
from snudda.utils.load import SnuddaLoad


class TestClusteredPruning(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        self.network_path_unclustered = os.path.join(os.path.dirname(__file__),
                                                     "networks", "network_testing_unclustered_pruning")
        self.network_path_clustered = os.path.join(os.path.dirname(__file__),
                                                   "networks", "network_testing_clustered_pruning")

        self.create_network(self.network_path_clustered)
        self.create_network(self.network_path_unclustered)

    def create_network(self, network_path):

        create_cube_mesh(file_name=os.path.join(network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        config_file = os.path.join(network_path, "network-config.json")
        position_file = os.path.join(network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(network_path, "voxels", "network-putative-synapses.hdf5")

        #  TODO: If d_view is None code run sin serial, add test parallel
        sp = SnuddaPlace(config_file=config_file, d_view=None, verbose=True)
        sp.place()

        rpn = RepositionNeurons(position_file=position_file)
        # Post synaptic neurons
        rpn.place(0, position=[0, 0, 0], rotation=np.eye(3))
        rpn.place(1, position=[0, 50e-6, 0], rotation=np.eye(3))
        rpn.place(2, position=[0, 100e-6, 0], rotation=np.eye(3))
        rpn.place(3, position=[0, 150e-6, 0], rotation=np.eye(3))

        ang = np.pi / 2
        R_z = np.array([[np.cos(ang), -np.sin(ang), 0],
                        [np.sin(ang), np.cos(ang), 0],
                        [0, 0, 1]])

        # Pre synaptic neurons
        rpn.place(4, position=[200e-6, 0, 0], rotation=R_z)
        rpn.place(5, position=[200e-6, 50e-6, 0], rotation=R_z)
        rpn.place(6, position=[200e-6, 100e-6, 0], rotation=R_z)
        rpn.place(7, position=[200e-6, 150e-6, 0], rotation=R_z)

        rpn.close()

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        sd = SnuddaDetect(network_path=network_path, rc=None,
                          hyper_voxel_size=130, verbose=True)
        sd.detect()

        if False:
            sd.plot_hyper_voxel(plot_neurons=True, fig_file_name="touch-detection-clusters-validation")

        sp = SnuddaPrune(network_path=network_path, verbose=True)
        sp.prune()

    def cluster_compare(self):

        clustered_data_loader = SnuddaLoad(network_file=self.network_path_clustered, load_synapses=True)
        clustered_synapses = clustered_data_loader.data["synapses"]
        clustered_voxel_size = clustered_data_loader.data["voxel_size"]
        clustered_simulation_origo = clustered_data_loader.data["simulation_origo"]

        nonclustered_data_loader = SnuddaLoad(network_file=self.network_path_unclustered, load_synapses=True)
        nonclustered_synapses = nonclustered_data_loader.data["synapses"]
        nonclustered_voxel_size = nonclustered_data_loader.data["voxel_size"]
        nonclustered_simulation_origo = nonclustered_data_loader.data["simulation_origo"]

        for post_id, pre_id in zip([0, 1, 2, 3], [4, 5, 6, 7]):
            clustered_idx = np.where(np.logical_and(clustered_synapses[:, 0] == pre_id, clustered_synapses[:, 1] == post_id))[0]
            nonclustered_idx = np.where(np.logical_and(nonclustered_synapses[:, 0] == pre_id, nonclustered_synapses[:, 1] == post_id))[0]

            clustered_synapse_coords = clustered_synapses[clustered_idx, 2:5] * clustered_voxel_size + clustered_simulation_origo
            nonclustered_synapse_coords = nonclustered_synapses[nonclustered_idx, 2:5] * nonclustered_voxel_size + nonclustered_simulation_origo

            clustered_synapse_dist = scipy.spatial.distance.cdist(clustered_synapse_coords, clustered_synapse_coords)
            nonclustered_synapse_dist = scipy.spatial.distance.cdist(nonclustered_synapse_coords, nonclustered_synapse_coords)

            clustered_n_close = np.sum(clustered_synapse_dist < (10e-6 / clustered_voxel_size))
            nonclustered_n_close = np.sum(nonclustered_synapse_dist < (10e-6 / nonclustered_voxel_size))

            print(f"clustered = {clustered_n_close}, nonclustered: {nonclustered_n_close}")
            self.assertTrue(clustered_n_close > nonclustered_n_close * 1.5)

    def test_clustered(self):

        data_loader = SnuddaLoad(network_file=self.network_path_clustered, load_synapses=True)
        synapses = data_loader.data["synapses"]
        voxel_size = data_loader.data["voxel_size"]
        simulation_origo = data_loader.data["simulation_origo"]

        for post_id, pre_id in zip([0, 1, 2, 3], [4, 5, 6, 7]):
            idx = np.where(np.logical_and(synapses[:, 0] == pre_id, synapses[:, 1] == post_id))[0]
            coords = synapses[idx, 2:5]*voxel_size + simulation_origo

            max_dist = np.linalg.norm(np.max(coords, axis=0) - np.min(coords, axis=0))
            self.assertTrue(max_dist < 50e-6)

    def test_unclustered(self):

        data_loader = SnuddaLoad(network_file=self.network_path_unclustered, load_synapses=True)
        synapses = data_loader.data["synapses"]
        voxel_size = data_loader.data["voxel_size"]
        simulation_origo = data_loader.data["simulation_origo"]

        for post_id, pre_id in zip([0, 1, 2, 3], [4, 5, 6, 7]):
            idx = np.where(np.logical_and(synapses[:, 0] == pre_id, synapses[:, 1] == post_id))[0]
            coords = synapses[idx, 2:5]*voxel_size + simulation_origo

            max_dist = np.linalg.norm(np.max(coords, axis=0) - np.min(coords, axis=0))
            self.assertTrue(max_dist > 120e-6)


if __name__ == '__main__':
    unittest.main()
