import h5py
import numpy as np
import scipy.spatial

from snudda import SnuddaLoad
from collections import OrderedDict

# TODO: !!!  ORPHANED FILE, DELETE IN NEXT CLEANUP!!! !!!

class ReplaceMorphology:

    def __init__(self, original_network, new_network, new_snudda_data_path):

        """ This function takes a network-neuron-positions.hdf5 file and replaces the morphologies by switching
            SNUDDA_DATA path to new location, and remapping morphologies to point to new morphologies. """


class MergeSynapseFiles:

    def __init__(self, parent_network_A, parent_network_B, new_network, merge_type=None):

        """ Merging the synapses in parent_network_A and parent_network_B, writing the new synapses to new_network.

            There are different merge styles available:

                "addition" - The neurons in the two networks are part of different networks,
                             just add all neurons together. This will renumber all neuron_id:s.

                "structural_changes" - The neurons in the parent networks must be related in a 1-1 fashion.
                                       This can be used for example to degenerate an existing network by replacing
                                       morphologies with degenerated versions. Important that neuron location is
                                       identical in the two networks.

                                       Synapses on the old morphology are retained if the dendrites and axons remain.
                                       New synapses are added if the new morphology dendrites or axon did not exist
                                       in the old morphologies.

                                       !!! OBS, make sure coordinates are transformed for all synapses since simulation
                                       origo might be different





        """

        self.parent_network_file_A = parent_network_A
        self.parent_network_file_B = parent_network_B

        self.new_network_file = new_network
        self.merge_type = merge_type

        self.parent_network_A_loader = None
        self.parent_network_B_loader = None

        self.parent_network_A_data = None
        self.parent_network_B_data = None

        self.h5libver = "latest"
        self.h5driver = "sec2"

        # TODO: We need to start storing SNUDDA_DATA in the save file (meta),
        #  so the user does not need to specify it all the time

        # Workflow
        # 1. Generate parent_network_A as normal (init, place, detect, prune)
        # 2. Use network-neuron-position.hdf5 from parent_network_A and generate a network-neuron-position.hdf5
        #    for parent_network_B, then do (detect, prune) on that network.
        # 3. Merge synapse files together.
        #

        pass

    def load_parent_networks(self):

        self.parent_network_A_loader = SnuddaLoad(self.parent_network_file_A)
        self.parent_network_B_loader = SnuddaLoad(self.parent_network_file_B)

        self.parent_network_A_data = self.parent_network_A_loader.data
        self.parent_network_B_data = self.parent_network_B_loader.data

    def setup_new_network_file(self):

        """ Create initial hdf5 for new network. """

        self.new_network_file = h5py.File(self.new_network_file, "w", libver=self.h5libver, driver=self.h5driver)

        # TODO: Figure out what to do with the "config" data, should we merge it or just leave it blank?


        # Calculate new simulation origo

        voxel_size = self.parent_network_A_data["voxelSize"]
        assert voxel_size == self.parent_network_B_data["voxelSize"], f"Voxel size mismatch between the parent networks"

        sim_origo_A = self.parent_network_A_data["simulationOrig"]
        sim_origo_B = self.parent_network_B_data["simulationOrig"]

        sim_origo_new = np.min(sim_origo_A, sim_origo_B)

        coordinate_transform_A = sim_origo_new - sim_origo_A
        coordinate_transform_B = sim_origo_new - sim_origo_B

        voxel_transform_A = np.round(coordinate_transform_A / voxel_size)
        voxel_transform_B = np.round(coordinate_transform_B / voxel_size)

        meta_group = self.new_network_file.create_group("meta")
        meta_group.create_dataset("simulationOrigo", data=sim_origo_new)
        meta_group.create_dataset("voxelSize", data=voxel_size)

        self.new_network_file.copy("neurons", self.parent_network_file_B)

    def filter_synapses:

    def create_neuron_remapping(self):

        # 1. Identify neurons that are "the same neuron"
        # 2. Identify removed neurons
        # 3. Identify added neurons

        pass




    def remapping_neurons(self, correspondence_distance = 5e-6):

        parent_A_KDtree = scipy.spatial.cKDTree(self.parent_network_A_data["neuronPositions"])
        parent_B_KDtree = scipy.spatial.cKDTree(self.parent_network_B_data["neuronPositions"])

        mapping_new_to_A = OrderedDict()
        mapping_new_to_B = OrderedDict()
        mapping_A_to_new = OrderedDict()
        mapping_B_to_new = OrderedDict()

        neuron_counter = 0

        for neuron_B in self.parent_network_B_data["neurons"]:
            dist, closest_idx = parent_A_KDtree.query(neuron_B["position"])
            neuron_A = self.parent_network_A_data["neurons"][closest_idx]

            if self.check_matching_data(neuron_A=neuron_B, neuron_B=neuron_B):
                mapping_new_to_A[neuron_counter] = neuron_A["neuronID"]
                mapping_A_to_new[neuron_A["neuronID"]] = neuron_counter

            mapping_new_to_B[neuron_counter] = neuron_B["neuronID"]
            mapping_B_to_new[neuron_B["neuronID"]] = neuron_counter
            neuron_counter += 1

        return mapping_new_to_A, mapping_A_to_new, mapping_new_to_B, mapping_B_to_new


    def check_matching_data(self, neuron_A, neuron_B, max_diff = 5-6):

        if np.linalg.norm(neuron_A["position"] - neuron_B["position"]) < max_diff:
            return False

        if neuron_A["name"] != neuron_B["name"]:
            return False

        if np.linalg.norm(neuron_A["rotation"] - neuron_B["rotation"], ord=np.inf) > max_diff:
            return False

        return True

    def merge_neuron_lists(self, :):

        pass

    def merge_synapse_lists(self, ):

        pass