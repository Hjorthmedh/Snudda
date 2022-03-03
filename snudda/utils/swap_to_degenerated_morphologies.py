import json
import os

import h5py
import numpy as np
from scipy.spatial import KDTree

from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.load import SnuddaLoad
from snudda.utils.snudda_path import snudda_simplify_path, snudda_parse_path


class SwapToDegenerateMorphologies:

    def __init__(self, original_network_file, new_network_file,
                 original_snudda_data_dir, new_snudda_data_dir,
                 morphology_mapping_file):

        """ This code replaces the neuron morphologies in the original network with user provided degenerated copies
            of the neurons. The synapses that are on removed dendritic or axonal branches will also be removed.
            The section ID and section X is also updated to match the new degenerated morphologies.

            Args:
                original_network_file (str) : Path to input network-synapses.hdf5
                new_network_file (str) : Path to output network-synapses.hdf5
                original_snudda_data_dir (str) : Path to SNUDDA_DATA for original network
                new_snudda_data_dir (str) : Path to SNUDDA_DATA for new network
                morphology_mapping_file (str) : JSON file that maps original morphology_key to new morphology_key
            """

        self.original_network_file = original_network_file
        self.new_network_file = new_network_file
        self.new_hdf5 = None

        self.original_snudda_data_dir = original_snudda_data_dir
        self.new_snudda_data_dir = new_snudda_data_dir

        self.original_network_loader = SnuddaLoad(self.original_network_file)
        self.original_hdf5 = self.original_network_loader.hdf5_file
        self.original_data = self.original_network_loader.data

        with open(morphology_mapping_file, "r") as f:
            self.morphology_mapping = json.load(f)

        self.morphology_map = dict()

    def write_new_file(self):

        self.new_hdf5 = h5py.File(self.new_network_file)
        h5py.Group.copy(source=self.original_hdf5["meta"], dest=self.new_hdf5["meta"])
        network_group = self.new_hdf5.create_group("network")
        h5py.Group.copy(source=self.original_hdf5["network/neurons"], dest=self.new_hdf5["network/neurons"])

        # Update parameter keys and morphology keys
        for idx, neuron_id in enumerate(self.new_hdf5["networks/neurons/neuronID"]):
            assert idx == neuron_id, "There should be no gaps in numbering."
            param_key, morph_key, neuron_path = self.find_morpology(neuron_id)
            self.new_hdf5["networks/neurons/neuronID/parameterKey"][idx] = param_key
            self.new_hdf5["networks/neurons/neuronID/morphologyKey"][idx] = morph_key
            self.new_hdf5["networks/neurons/neuronID/neuronPath"][idx] = neuron_path

    def synapse_iterator(self, data_type=None):

        """ Each iteration will return the synapses between one pair of neurons. """

        data_loc = "network/synapses"
        if data_type is not None and data_type == "gapJunctions"
            data_loc = "network/gapJunctions"

        start_idx = 0
        next_idx = 0
        num_synapses = self.original_hdf5[data_loc][()]

        assert num_synapses == self.original_hdf5[data_loc].shape[0]

        last_pre = self.original_hdf5[data_loc][0, 0]
        last_post = self.original_hdf5[data_loc][0, 1]

        while next_idx < num_synapses:
            while next_idx < num_synapses \
                    and self.original_hdf5[data_loc][next_idx, 0] == last_pre \
                    and self.original_hdf5[data_loc][next_idx, 1] == last_post:
                next_idx += 1

            if next_idx < num_synapses:
                last_pre = self.original_hdf5[data_loc][next_idx, 0]
                last_post = self.original_hdf5[data_loc][next_idx, 1]

            yield self.original_hdf5[data_loc][start_idx:next_idx, :]

    def gap_junction_iterator(self):

        """ Each iteration will return the gap junctions between one pair of neurons. """

        return self.synapse_iterator(data_type="gapJunctions")

    def filter_synapses(self):
        # First version, will keep all synapses in memory to write a more efficient file
        new_synapses = np.zeros(self.original_hdf5["network/synapses"].shape,
                                dtype=self.original_hdf5["network/synapses"].dtype)
        syn_ctr = 0

        for synapses in self.synapse_iterator():
            new_syn = self.filter_synapses_helper(synapses)
            new_synapses[syn_ctr:syn_ctr+new_syn.shape[0]] = new_syn
            syn_ctr += new_syn.shape[0]

        self.new_hdf5["network"].create_dataset("synapses", data=new_synapses, compression="lzf")

        num_synapses = np.zeros((1,), dtype=np.uint64)
        self.new_hdf5["network"].create_dataset("nSynapses", data=new_synapses.shape[0], dtype=np.uint64)

    def filter_gap_junctions(self):
        # First version, will keep all synapses in memory to write a more efficient file
        new_gap_junctions = np.zeros(self.original_hdf5["network/gapJunctions"].shape,
                                     dtype=self.original_hdf5["network/gapJunctions"].dtype)
        gj_ctr = 0

        for gj in self.gap_junction_iterator():
            new_gj = self.filter_gap_junctions_helper(gj)
            new_gap_junctions[gj_ctr:gj_ctr+new_gap_junctions.shape[0]] = new_gj
            gj_ctr += new_gj.shape[0]

        self.new_hdf5["network"].create_dataset("gapJunctions", data=new_gap_junctions, compression="lzf")

        num_gap_junctions = np.zeros((1,), dtype=np.uint64)
        self.new_hdf5["network"].create_dataset("nGapJunctions", data=new_gap_junctions.shape[0], dtype=np.uint64)

    def get_morphology(self, neuron_id):

        param_key = self.new_hdf5["networks/neurons/neuronID/parameterKey"][neuron_id]
        morph_key = self.new_hdf5["networks/neurons/neuronID/morphologyKey"][neuron_id]
        neuron_path = self.new_hdf5["networks/neurons/neuronID/neuronPath"][neuron_id]

        neuron_prototype = NeuronPrototype(neuron_path=neuron_path)
        neuron = neuron_prototype.get_morphology(parameter_key=param_key, morph_key=morph_key)
        neuron.place(self.original_data["neuronPositions"][neuron_id, :])

        return neuron


    def filter_synapses_helper(self, synapses, max_dist=5e-6):
        pre_id = synapses[0, 0]
        post_id = synapses[0, 1]

        keep_synapses = np.ones((synapses.shape[0],), dtype=bool)

        # Sanity check that we only have synapses between one pair of neurons
        assert (synapses[:, 0] == pre_id).all()
        assert (synapses[:, 1] == post_id).all()

        pre_neuron = self.get_morphology(pre_id)
        post_neuron = self.get_morphology(post_id)

        if pre_neuron.axon.size > 0:
            # If there is an axon, use it to filter
            pre_axon = KDTree(pre_neuron.axon)
        else:
            pre_axon = None

        post_dend = KDTree(post_neuron.dend)

        for idx, syn in enumerate(synapses):
            syn_coord = syn[2:5] * self.original_data["voxelSize"] + self.original_data["simulationOrigo"]

            closest_post_point = post_dend.querty(syn_coord)
            if np.linalg.norm(syn_coord - closest_post_point) > max_dist:
                keep_synapses[idx] = False

            if pre_axon and keep_synapses[idx]:
                closest_pre_point = pre_axon.querty(syn_coord)
                if np.linalg.norm(syn_coord - closest_pre_point) > max_dist:
                    keep_synapses[idx] = False

        return synapses[keep_synapses, :]

    def filter_gap_junctions_helper(self, gap_junctions, max_dist=5e-6):

        # Gap junctions are bidirectional, but call them pre and post out of habit
        pre_id = gap_junctions[0, 0]
        post_id = gap_junctions[0, 1]

        keep_gap_junctions = np.ones((gap_junctions.shape[0],), dtype=bool)

        assert (gap_junctions[:, 0] == pre_id).all()
        assert (gap_junctions[:, 1] == post_id).all()

        pre_neuron = self.get_morphology(pre_id)
        post_neuron = self.get_morphology(post_id)

        pre_dend = KDTree(pre_neuron.dend)
        post_dend = KDTree(post_neuron.dend)

        for idx, gj in enumerate(gap_junctions):
            gj_coord = gj[6:9] * self.original_data["voxelSize"] + self.original_data["simulationOrigo"]

            closest_pre_point = pre_dend.querty(gj_coord)
            closest_post_point = post_dend.querty(gj_coord)

            if np.linalg.norm(gj_coord - closest_pre_point) > max_dist:
                keep_gap_junctions[idx] = False

            if np.linalg.norm(gj_coord - closest_post_point) > max_dist:
                keep_gap_junctions[idx] = False

        return gap_junctions[keep_gap_junctions, :]


    # TODO!! We also need to filter external input

    def find_morpology(self, neuron_id):

        """ Given neuron_id it looks up what was old morphology, then looks for a new morphology with the same name
            and finds the corresponding morphology_key. It also tries to pick a parameter_key """

        orig_morph_key = self.original_data["neurons"][neuron_id]["morphologyKey"]
        orig_param_key = self.original_data["neurons"][neuron_id]["parameterKey"]

        # We assume the new morphology is in the same relative path, but using a different SNUDDA_DATA
        orig_neuron_path = self.original_data["neurons"][neuron_id]["neuronPath"]
        orig_simple_path = snudda_simplify_path(orig_neuron_path, self.original_snudda_data_dir)
        new_neuron_path = snudda_parse_path(orig_simple_path, self.new_snudda_data_dir)

        # Here we assume there is a meta.json file
        with open(os.path.join(orig_neuron_path, "meta.json"), "r") as f:
            orig_meta_info = json.load(f)

        with open(os.path.join(new_neuron_path, "meta.json"), "r") as f:
            new_meta_info = json.load(f)

        # Find the parameter_key and morphology_key corresponding to the morphology used. Note that there might be
        # multiple parameter_key:s valid for a morphology, so then we randomly pick one of those

        orig_morph_name = orig_meta_info[orig_param_key][orig_morph_key]["morphology"]

        possible_keys = []

        for param_key, param_data in new_meta_info.items():
            for morph_key, morph_data in param_data():
                morph_name = morph_data["morphology"]
                if orig_morph_name == morph_name:
                    possible_keys.append((param_key, morph_key))

        assert len(possible_keys) > 0, \
            f"No morphology matching for {orig_morph_name}, unable to pick parameter_key and morphology_key"

        # Pick one of the key pairs
        idx = np.random.randint(low=0, high=len(possible_keys))
        new_param_key, new_morph_key = possible_keys[idx]

        return new_param_key, new_morph_key, new_neuron_path

    def create_morphology_mapping(self,
                                  original_neuron_path,
                                  original_parameter_key,
                                  original_morphology_key,
                                  new_neuron_path,
                                  new_parameter_key,
                                  new_morphology_key):

        # For each morphology pair, create a lookup for section ID
        #

        # TODO: Fortsätt här. Skapa self.morphology_map
        #       [org_morph_key][new_morph_key][section_id][section_x]
        #       ---> new_section_id , new_section_x    eller None om ej kvar

    # Idea:
    # 1. Use a KD tree for each neurons dendrites and axons (where available)
    # 2. For each synapse, check if it is still close to an existing point on the presynaptic axon and
    #    postsynaptic dendrite.
    # 3. Remap the section ID and section X for the remaining
    # 4. Update the config file
