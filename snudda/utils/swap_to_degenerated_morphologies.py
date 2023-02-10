import json
import os

import h5py
import numpy as np
from scipy.spatial import cKDTree

from argparse import ArgumentParser, RawTextHelpFormatter

from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.load import SnuddaLoad
from snudda.utils.snudda_path import snudda_simplify_path, snudda_parse_path


class SwapToDegeneratedMorphologies:

    def __init__(self, original_network_file, new_network_file,
                 original_snudda_data_dir, new_snudda_data_dir,
                 original_input_file=None, new_input_file=None,
                 filter_axon=False):

        """ This code replaces the neuron morphologies in the original network with user provided degenerated copies
            of the neurons. The synapses that are on removed dendritic will also be removed.
            The section ID and section X is also updated to match the new degenerated morphologies.

            Args:
                original_network_file (str) : Path to input network-synapses.hdf5
                new_network_file (str) : Path to output network-synapses.hdf5
                original_snudda_data_dir (str) : Path to SNUDDA_DATA for original network
                new_snudda_data_dir (str) : Path to SNUDDA_DATA for new network
            """

        self.original_network_file = original_network_file
        self.new_network_file = new_network_file
        self.new_hdf5 = None

        self.original_snudda_data_dir = original_snudda_data_dir
        self.new_snudda_data_dir = new_snudda_data_dir

        assert self.original_snudda_data_dir != self.new_snudda_data_dir, \
            f"SNUDDA_DATA_DIR should be different for WT and degenerated"

        self.original_input_file = original_input_file
        self.new_input_file = new_input_file

        self.original_network_loader = SnuddaLoad(self.original_network_file, load_synapses=False)
        self.old_hdf5 = self.original_network_loader.hdf5_file
        self.old_data = self.original_network_loader.data

        self.morphology_map = dict()
        self.section_lookup = dict()
        self.neuron_cache_key = dict()
        self.neuron_cache_id = dict()
        self.key_lookup = dict()

        self.create_section_lookup()
        self.kd_tree_cache = dict()

        self.filter_axon = filter_axon
        self.old_simulation_origo = self.old_hdf5["meta/simulationOrigo"][()]
        self.old_voxel_size = self.old_hdf5["meta/voxelSize"][()]

        self.has_axon_density = np.zeros((self.original_network_loader.data["nNeurons"],), dtype=bool)
        for neuron_id, _ in enumerate(self.original_network_loader.data["neurons"]):
            if "axonDensity" in self.original_network_loader.data["neurons"][neuron_id] \
                    and self.original_network_loader.data["neurons"][neuron_id]["axonDensity"] is not None:
                self.has_axon_density[neuron_id] = True

    def close(self):

        if self.new_hdf5:
            self.new_hdf5.close()

        if self.old_hdf5:
            self.old_hdf5.close()

    def write_new_network_file(self):

        if not os.path.isdir(os.path.dirname(self.new_network_file)):
            print(f"Creating directory {os.path.dirname(self.new_network_file)}")
            os.mkdir(os.path.dirname(self.new_network_file))

        print(f"Writing new network to {self.new_network_file}")
        self.new_hdf5 = h5py.File(self.new_network_file, "w")
        self.old_hdf5.copy(source=self.old_hdf5["meta"], dest=self.new_hdf5)

        if len(self.new_snudda_data_dir) > len(self.original_snudda_data_dir):
            del self.new_hdf5["meta/snuddaData"]
            self.new_hdf5["meta"].create_dataset("snuddaData", data=self.new_snudda_data_dir)
        else:
            self.new_hdf5["meta/snuddaData"][()] = self.new_snudda_data_dir
        network_group = self.new_hdf5.create_group("network")
        self.old_hdf5.copy(source=self.old_hdf5["network/neurons"], dest=self.new_hdf5["network"])

        # Update parameter keys and morphology keys
        for idx, neuron_id in enumerate(self.new_hdf5["network/neurons/neuronID"]):
            assert idx == neuron_id, "There should be no gaps in numbering."
            param_key, morph_key, neuron_path, param_id, morph_id = self.find_morpology(neuron_id)
            self.new_hdf5[f"network/neurons/parameterKey"][idx] = param_key
            self.new_hdf5[f"network/neurons/morphologyKey"][idx] = morph_key
            # self.new_hdf5[f"network/neurons/parameterID"][idx] = param_id
            # self.new_hdf5[f"network/neurons/morphologyID"][idx] = morph_id

        self.filter_synapses(filter_axon=self.filter_axon)
        self.filter_gap_junctions()

    def synapse_iterator(self, data_type=None, load_synapses=True, synapses=None):

        """ Each iteration will return the synapses between one pair of neurons.

            Normally loads the data from the hdf5 file, but if synapses is given,
            then it will use the provided matrix instead. The synapses argument is used by
            post_degeneration_pruning in swap_to_degenerated_morphology.

        """

        start_idx = 0
        next_idx = 0

        if synapses is None:
            data_loc = "network/synapses"
            if data_type is not None and data_type == "gapJunctions":
                data_loc = "network/gapJunctions"
                num_synapses = self.old_hdf5["network/nGapJunctions"][()]
            else:
                num_synapses = self.old_hdf5["network/nSynapses"][()]

            # Load synapses into memory
            if load_synapses:
                print("Loading synapses into memory.")
                synapses = self.old_hdf5[data_loc][()].copy()
            else:
                print("Reading synapses from disk (slower)")
                synapses = self.old_hdf5[data_loc]

            assert num_synapses == synapses.shape[0]

        else:
            assert data_type is None, (f"If synapses is given, the data will not be loaded from the hdf5 file, "
                                       f"leave data_type as None")
            num_synapses = synapses.shape[0]

        if num_synapses > 0:
            last_pre = synapses[0, 0]
            last_post = synapses[0, 1]

        while next_idx < num_synapses:
            while next_idx < num_synapses \
                    and synapses[next_idx, 0] == last_pre \
                    and synapses[next_idx, 1] == last_post:
                next_idx += 1
                if next_idx % 1000000 == 0:
                    print(f"{next_idx} / {num_synapses}")

            if next_idx < num_synapses:
                last_pre = synapses[next_idx, 0]
                last_post = synapses[next_idx, 1]

            yield synapses[start_idx:next_idx, :]
            start_idx = next_idx

        print(f"{next_idx} / {num_synapses}")

    def gap_junction_iterator(self, load_gap_junctions=True):

        """ Each iteration will return the gap junctions between one pair of neurons. """

        return self.synapse_iterator(data_type="gapJunctions", load_synapses=load_gap_junctions)

    def filter_synapses(self, filter_axon=False):
        # First version, will keep all synapses in memory to write a more efficient file
        new_synapses = np.zeros(self.old_hdf5["network/synapses"].shape,
                                dtype=self.old_hdf5["network/synapses"].dtype)
        syn_ctr = 0

        for synapses in self.synapse_iterator():

            new_syn = self.filter_synapses_helper(synapses, filter_axon=filter_axon)
            new_synapses[syn_ctr:syn_ctr + new_syn.shape[0]] = new_syn
            syn_ctr += new_syn.shape[0]

        self.new_hdf5["network"].create_dataset("synapses", data=new_synapses[:syn_ctr, :], compression="lzf")

        num_synapses = np.zeros((1,), dtype=np.uint64)
        self.new_hdf5["network"].create_dataset("nSynapses", data=syn_ctr, dtype=np.uint64)

        print(f"Keeping {self.new_hdf5['network/nSynapses'][()]} "
              f"out of {self.old_hdf5['network/nSynapses'][()]} synapses "
              f"({self.new_hdf5['network/nSynapses'][()] / self.old_hdf5['network/nSynapses'][()]*100:.3f} %)")

    def filter_gap_junctions(self):
        # First version, will keep all synapses in memory to write a more efficient file
        new_gap_junctions = np.zeros(self.old_hdf5["network/gapJunctions"].shape,
                                     dtype=self.old_hdf5["network/gapJunctions"].dtype)
        gj_ctr = 0

        for gj in self.gap_junction_iterator():
            new_gj = self.filter_gap_junctions_helper(gj)
            new_gap_junctions[gj_ctr:gj_ctr + new_gj.shape[0]] = new_gj
            gj_ctr += new_gj.shape[0]

        self.new_hdf5["network"].create_dataset("gapJunctions", data=new_gap_junctions[:gj_ctr, :], compression="lzf")

        num_gap_junctions = np.zeros((1,), dtype=np.uint64)
        self.new_hdf5["network"].create_dataset("nGapJunctions", data=gj_ctr, dtype=np.uint64)

        try:
            print(f"Keeping {self.new_hdf5['network/nGapJunctions'][()]} "
                  f"out of {self.old_hdf5['network/nGapJunctions'][()]} gap junctions "
                  f"({self.new_hdf5['network/nGapJunctions'][()] / max(1, self.old_hdf5['network/nGapJunctions'][()])*100:.3f} %)")
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

    def get_morphology(self, neuron_id=None, hdf5=None, neuron_path=None, snudda_data=None,
                       parameter_key=None, morphology_key=None,
                       neuron_cache_id=None, neuron_cache_key=None):

        """ If neuron_id is given, that neuron will be loaded."""

        if neuron_cache_id is None:
            neuron_cache_id = self.neuron_cache_id

        if neuron_cache_key is None:
            neuron_cache_key = self.neuron_cache_key

        if neuron_id is not None and neuron_id in neuron_cache_id:
            return neuron_cache_id[neuron_id]

        if neuron_path is not None:
            neuron_path = snudda_parse_path(neuron_path, snudda_data)

        if (neuron_path is not None and parameter_key is not None and morphology_key is not None
                and (neuron_path, parameter_key, morphology_key) in self.neuron_cache_key):
            return self.neuron_cache_key[(neuron_path, parameter_key, morphology_key)]

        if neuron_id is not None and hdf5 is not None:
            neuron_path = snudda_parse_path(SnuddaLoad.to_str(hdf5["network/neurons/neuronPath"][neuron_id]),
                                            snudda_data=snudda_data)
            parameter_key = SnuddaLoad.to_str(hdf5["network/neurons/parameterKey"][neuron_id])
            morphology_key = SnuddaLoad.to_str(hdf5["network/neurons/morphologyKey"][neuron_id])
            neuron_name = SnuddaLoad.to_str(hdf5["network/neurons/name"][neuron_id])
        else:
            neuron_name = None

        # assert neuron_path is not None and parameter_key is not None and morphology_key is not None, \
        #     "Either provide neuron_id, hdf5 or the three neuron_path, parameter_key and morphology_key"

        if neuron_id is not None and hdf5:
            pos = hdf5["network/neurons/position"][neuron_id, :]
            rot = hdf5["network/neurons/rotation"][neuron_id, :].reshape(3, 3)
        else:
            pos = None
            rot = None

        # TODO: Check if we need to pass snudda_data here!!!
        neuron_prototype = NeuronPrototype(neuron_path=neuron_path, neuron_name=neuron_name, snudda_data=snudda_data)
        neuron = neuron_prototype.clone(parameter_key=parameter_key, morphology_key=morphology_key,
                                        position=pos, rotation=rot)

        if neuron_id is not None:
            neuron_cache_id[neuron_id] = neuron
        else:
            neuron_cache_key[(neuron_path, parameter_key, morphology_key)] = neuron

        return neuron

    def filter_synapses_helper(self, synapses, filter_axon=False):

        """ This currently only takes DENDRITE degeneration into account.
            AXON degeneration is not taken into account for synapse removal.
        """

        pre_id = synapses[0, 0]
        post_id = synapses[0, 1]

        assert (synapses[:, 0] == pre_id).all()
        assert (synapses[:, 1] == post_id).all()

        old_sec_id = synapses[:, 9]
        old_sec_x = synapses[:, 10]

        keep_idx, new_sec_id, new_sec_x \
            = self.remap_sections_helper(neuron_id=post_id, old_sec_id=old_sec_id, old_sec_x=old_sec_x/1000.0)

        edited_synapses = synapses.copy()

        edited_synapses[:, 9] = new_sec_id
        edited_synapses[:, 10] = new_sec_x * 1000

        if filter_axon:
            filtered_synapses = self.filter_axonal_synapses_helper(edited_synapses[keep_idx, :])
        else:
            filtered_synapses = edited_synapses[keep_idx, :]

        return filtered_synapses

    def filter_axonal_synapses_helper(self, synapses, max_dist=5.41e-6):

        """ Filter the synapses that have the axon degeneration, presynaptic neurons without axons are ignored. """

        pre_id = synapses[:, 0]
        synapse_voxels = synapses[:, 2:5]

        synapse_coordinates = self.old_simulation_origo + self.old_voxel_size*synapse_voxels

        keep_idx = np.ones((synapses.shape[0],), dtype=bool)

        loaded_pid = None
        axon_tree = None

        # Loop through all the synapses, if they have an axon check that there is an axonal point close to synapse
        for idx, (pid, coord) in enumerate(zip(pre_id, synapse_coordinates)):

            if pid != loaded_pid:
                morph = self.get_morphology(neuron_id=pid, hdf5=self.new_hdf5, snudda_data=self.new_snudda_data_dir)
                loaded_pid = pid

                # Has axon?
                if 2 in morph.morphology_data["neuron"].sections and len(morph.morphology_data["neuron"].sections[2]) > 0:
                    axon_tree = self.get_kd_tree(morph, "axon")

                    if self.has_axon_density[pid]:
                        axon_tree = None
                        # print(f"Warning: Axon and axonal density specified for neuron {pid}")
                        # print(f"Ignoring morphology --- for now, will change behaviour in the future.")
                else:
                    axon_tree = None

            if axon_tree is None:
                # No pre-synaptic axon exists for this neuron (no info, so keep synapse)
                continue

            closest_dist, closest_point = axon_tree.query(coord)

            # Here we assume one query per point, so closest_dist should be a scalar
            if closest_dist > max_dist:
                keep_idx[idx] = 0

        return synapses[keep_idx, :].copy()

    def filter_gap_junctions_helper(self, gap_junctions):

        # Gap junctions are bidirectional, but call them pre and post out of habit
        pre_id = gap_junctions[0, 0]
        post_id = gap_junctions[0, 1]

        assert (gap_junctions[:, 0] == pre_id).all()
        assert (gap_junctions[:, 1] == post_id).all()

        old_pre_sec_id = gap_junctions[:, 2]
        old_pre_sec_x = gap_junctions[:, 4]
        old_post_sec_id = gap_junctions[:, 3]
        old_post_sec_x = gap_junctions[:, 5]

        keep_idx_pre, new_pre_sec_id, new_pre_sec_x \
            = self.remap_sections_helper(neuron_id=pre_id, old_sec_id=old_pre_sec_id, old_sec_x=old_pre_sec_x/1000)

        keep_idx_post, new_post_sec_id, new_post_sec_x \
            = self.remap_sections_helper(neuron_id=post_id, old_sec_id=old_post_sec_id, old_sec_x=old_post_sec_x/1000)

        keep_idx = np.intersect1d(keep_idx_pre, keep_idx_post)

        edited_gap_junctions = gap_junctions.copy()
        edited_gap_junctions[:, 2] = new_pre_sec_id
        edited_gap_junctions[:, 3] = new_post_sec_id
        edited_gap_junctions[:, 4] = new_pre_sec_x * 1000
        edited_gap_junctions[:, 5] = new_post_sec_x * 1000

        return edited_gap_junctions[keep_idx, :]

    def find_old_morphology(self, neuron_id):

        if neuron_id not in self.key_lookup:
            orig_morph_key = SnuddaLoad.to_str(self.old_data["neurons"][neuron_id]["morphologyKey"])
            orig_param_key = SnuddaLoad.to_str(self.old_data["neurons"][neuron_id]["parameterKey"])
            orig_neuron_path = SnuddaLoad.to_str(self.old_data["neurons"][neuron_id]["neuronPath"])

            self.key_lookup[neuron_id] = orig_param_key, orig_morph_key, orig_neuron_path
        else:
            orig_param_key, orig_morph_key, orig_neuron_path = self.key_lookup[neuron_id]

        return orig_param_key, orig_morph_key, orig_neuron_path

    def find_morpology(self, neuron_id):

        """ Given neuron_id it looks up what was old morphology, then looks for a new morphology with the same name
            and finds the corresponding morphology_key. It also tries to pick a parameter_key """

        orig_morph_key = SnuddaLoad.to_str(self.old_data["neurons"][neuron_id]["morphologyKey"])
        orig_param_key = SnuddaLoad.to_str(self.old_data["neurons"][neuron_id]["parameterKey"])

        # We assume the new morphology is in the same relative path, but using a different SNUDDA_DATA
        orig_simple_path = SnuddaLoad.to_str(self.old_data["neurons"][neuron_id]["neuronPath"])
        orig_neuron_path = snudda_parse_path(orig_simple_path, os.path.realpath(self.original_snudda_data_dir))
        # new_neuron_path = snudda_parse_path(orig_simple_path, os.path.realpath(self.new_snudda_data_dir))
        new_neuron_path = orig_neuron_path.replace(os.path.realpath(self.original_snudda_data_dir),
                                                   os.path.realpath(self.new_snudda_data_dir))

        if orig_morph_key is None or orig_morph_key == '':
            # Only a single morphology

            original_morphology_id = 0  # self.old_data["neurons"][neuron_id]["morphologyID"]
            original_parameter_id = 0  # self.old_data["neurons"][neuron_id]["parameterID"]

            return '', '', new_neuron_path, original_parameter_id, original_morphology_id

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
            for morph_key, morph_data in param_data.items():
                morph_name = morph_data["morphology"]
                if orig_morph_name == morph_name:
                    possible_keys.append((param_key, morph_key))

        assert len(possible_keys) > 0, \
            f"No morphology matching for {orig_morph_name}, unable to pick parameter_key and morphology_key"

        # Pick one of the key pairs
        idx = np.random.randint(low=0, high=len(possible_keys))
        new_param_key, new_morph_key = possible_keys[idx]

        # We also need parameter_id and morphology_id
        parameter_id = np.where([x == new_param_key for x in new_meta_info.keys()])[0][0]
        morphology_id = np.where([x == new_morph_key for x in new_meta_info[new_param_key].keys()])[0][0]

        return new_param_key, new_morph_key, new_neuron_path, parameter_id, morphology_id

    def get_sec_location(self, coords, neuron_path, snudda_data,
                         parameter_key, morphology_key, max_dist=5.41e-6):

        raise DeprecationWarning("This function is no longer used. It is based on old NeuronMorphology -- REMOVE?")
        assert False, "Do not run this!"

        morph = self.get_morphology(neuron_path=neuron_path,
                                    parameter_key=parameter_key,
                                    morphology_key=morphology_key,
                                    snudda_data=snudda_data)

        dend = self.get_kd_tree(morph, "dend")
        sec_id = np.zeros((dend.shape[0],), dtype=int)
        sec_x = np.zeros((dend.shape[0],))

        coord_to_sec_id_x = dict()
        for link, sec_id, sec_x in zip(morph.dend_links, morph.dend_sec_id, morph.dend_sec_x):
            coord = morph.dend[link[1], :]

            coord_to_sec_id_x[coord] = (sec_id, sec_x[1])

        for idx, coord in enumerate(coords):
            closest_dist, closest_point = dend.query(coord)

            if closest_dist <= max_dist:
                syn_sec_id, syn_sec_x = coord_to_sec_id_x[morph.dend[closest_point, :3]]
                sec_id[idx] = syn_sec_id
                sec_x[idx] = syn_sec_x
            else:
                sec_id[idx] = np.nan
                sec_x[idx] = np.nan

        return sec_id, sec_x

    def write_new_input_file(self, remap_removed_input=False, remapped_fraction=1.0):

        if self.original_input_file is None:
            print("No input file supplied, no new input file created")
            return

        if not os.path.isdir(os.path.dirname(self.new_input_file)):
            os.mkdir(os.path.dirname(self.new_input_file))

        print(f"Writing new input data to {self.new_input_file}")

        old_input = h5py.File(self.original_input_file, "r")
        new_input = h5py.File(self.new_input_file, "w")
        old_input.copy(source=old_input["config"], dest=new_input)

        new_input.create_group("input")

        for neuron in old_input["input"].keys():
            neuron_group = new_input["input"].create_group(neuron)

            old_n = 0
            new_n = 0
            remap_n = 0

            for input_type in old_input["input"][neuron].keys():
                # Note: This code assumes these are real neuron and not just virtual neurons.
                #       it does not handle the virtual neuron case here.
                old_input_data = old_input["input"][neuron][input_type]

                keep_idx, new_sec_id, new_sec_x \
                    = self.remap_sections_helper(neuron_id=int(neuron),
                                                 old_sec_id=old_input_data.attrs["sectionID"],
                                                 old_sec_x=old_input_data.attrs["sectionX"])

                if len(keep_idx) == 0 and not (remap_removed_input and remapped_fraction > 0):
                    continue

                input_group = neuron_group.create_group(input_type)

                if remap_removed_input:

                    # We need to find new positions for input marked as removed
                    morph = self.get_morphology(neuron_id=int(neuron), hdf5=self.new_hdf5,
                                                snudda_data=self.new_snudda_data_dir)

                    n_remap = len(old_input_data.attrs["sectionID"]) - len(keep_idx)
                    idx_remap = sorted(list(set(np.arange(0, len(old_input_data.attrs["sectionID"]))) - set(keep_idx)))

                    if remapped_fraction < 1.0:
                        n_remap = int(np.round(len(idx_remap) * remapped_fraction))
                        idx_remap = sorted(list(np.random.permutation(idx_remap)[:n_remap]))

                    try:
                        synapse_density = SnuddaLoad.to_str(old_input_data.attrs["synapseDensity"])
                    except:
                        import traceback
                        print(traceback.format_exc())
                        import pdb
                        pdb.set_trace()

                    xyz, sec_id, sec_x, dist_to_soma = morph.dendrite_input_locations(synapse_density_str=synapse_density,
                                                                                      num_locations=n_remap,
                                                                                      rng=self.rng,
                                                                                      cluster_size=1,
                                                                                      cluster_spread=None)

                    old_n += old_input_data['spikes'].shape[0]
                    new_n += len(keep_idx)
                    remap_n += len(idx_remap)

                    keep_idx2 = sorted(list(set(keep_idx).union(set(idx_remap))))

                    # Same spikes as before
                    spike_set = input_group.create_dataset("spikes", data=old_input_data["spikes"][keep_idx2, :],
                                                           compression="gzip", dtype=np.float32)
                    spike_set.attrs["nSpikes"] = old_input_data["spikes"].attrs["nSpikes"][keep_idx2].astype(np.int32)

                    # New locations for the remapped synapses
                    new_sec_id[idx_remap] = sec_id
                    input_group.attrs["sectionID"] = new_sec_id[keep_idx2].astype(np.int16)

                    new_sec_x[idx_remap] = sec_x
                    input_group.attrs["sectionX"] = new_sec_x[keep_idx2].astype(np.float16)

                    input_group.attrs["parameterID"] = old_input_data.attrs["parameterID"][keep_idx2].astype(np.int)

                    updated_dist = old_input_data.attrs["distanceToSoma"].copy()
                    updated_dist[idx_remap] = dist_to_soma
                    input_group.attrs["distanceToSoma"] = updated_dist[keep_idx2].astype(np.float16)

                    for data_name in ["freq", "correlation", "jitter", "synapseDensity", "start", "end", "conductance",
                                      "populationUnitID", "populationUnitSpikes", "generator",
                                      "modFile", "parameterFile", "parameterList"]:

                        if data_name in old_input_data.attrs:
                            input_group.attrs[data_name] = old_input_data.attrs[data_name]
                        elif data_name in old_input_data["spikes"].attrs:
                            input_group["spikes"].attrs[data_name] = old_input_data["spikes"].attrs[data_name]
                        elif data_name in old_input_data:
                            old_input_data.copy(source=old_input_data[data_name], dest=input_group)

                else:
                    input_group.create_dataset("spikes", data=old_input_data["spikes"][keep_idx, :],
                                               compression="gzip", dtype=np.float32)
                    input_group["spikes"].attrs["nSpikes"] = old_input_data["spikes"].attrs["nSpikes"][keep_idx]
                    input_group.attrs["sectionID"] = new_sec_id[keep_idx]
                    input_group.attrs["sectionX"] = new_sec_x[keep_idx]
                    input_group.attrs["parameterID"] = old_input_data.attrs["parameterID"][keep_idx]
                    input_group.attrs["distanceToSoma"] = old_input_data.attrs["distanceToSoma"][keep_idx]

                    for data_name in ["freq", "correlation", "jitter", "synapseDensity", "start", "end", "conductance",
                                      "populationUnitID", "populationUnitSpikes", "generator",
                                      "modFile", "parameterFile", "parameterList"]:

                        if data_name in old_input_data.attrs:
                            input_group.attrs[data_name] = old_input_data.attrs[data_name]
                        elif data_name in old_input_data["spikes"].attrs:
                            input_group["spikes"].attrs[data_name] = old_input_data["spikes"].attrs[data_name]
                        elif data_name in old_input_data:
                            old_input_data.copy(source=old_input_data[data_name], dest=input_group)

                    old_n += old_input_data['spikes'].shape[0]
                    new_n += len(keep_idx)

            print(f"Processed input to {self.old_data['neurons'][int(neuron)]['name']} ({neuron}), "
                  f"keeping {new_n} out of {old_n} inputs (plus remapping {remap_n} inputs)"
                  f"({new_n / max(old_n,1) * 100 :.2f} %) ({remap_n / max(old_n,1)*100 :.2f} % remapped)")

        self.close()

    def get_kd_tree(self, neuron, tree_type, kd_tree_cache=None):

        if kd_tree_cache is None:
            kd_tree_cache = self.kd_tree_cache

        morph_type_lookup = {"axon": 2, "dend": 3}

        if (neuron, tree_type) not in kd_tree_cache:
            morph_type = morph_type_lookup[tree_type]

            idx = np.where(neuron.morphology_data["neuron"].section_data[:, 2] == morph_type)[0]
            coords = neuron.morphology_data["neuron"].geometry[idx, :3]

            if coords.size > 0:
                kd_tree_cache[(neuron, tree_type)] = cKDTree(coords)
            else:
                kd_tree_cache[(neuron, tree_type)] = None

        return kd_tree_cache[(neuron, tree_type)]

    def create_section_lookup(self):

        for neuron_id in self.old_hdf5["network/neurons/neuronID"]:
            self.create_section_lookup_helper(neuron_id)

    # Note this does not take AXON degeneration into account, need to do that separately
    def create_section_lookup_helper(self, neuron_id):

        old_param_key, old_morph_key, old_path = self.find_old_morphology(neuron_id=neuron_id)

        new_param_key, new_morph_key, new_path, _, _ = self.find_morpology(neuron_id=neuron_id)

        if (old_param_key, old_morph_key, old_path) in self.section_lookup:
            # We already have the morphology in the section_lookup
            return

        old_morph = self.get_morphology(parameter_key=old_param_key, morphology_key=old_morph_key,
                                        neuron_path=old_path, snudda_data=self.original_snudda_data_dir)
        new_morph = self.get_morphology(parameter_key=new_param_key, morphology_key=new_morph_key,
                                        neuron_path=new_path, snudda_data=self.new_snudda_data_dir)
        coord_to_sec_id_x = dict()
        old_dend_idx = np.where(old_morph.morphology_data["neuron"].section_data[:, 2] == 3)[0]
        for idx in old_dend_idx:
            coord = (old_morph.morphology_data["neuron"].geometry[idx, :3] * 1e9).astype(int)
            assert (coord[0], coord[1], coord[2]) not in coord_to_sec_id_x, \
                f"Coordinates {(coord[0], coord[1], coord[2])} already exists in {coord_to_sec_id_x}"

            old_sec_id = old_morph.morphology_data["neuron"].section_data[idx, 0]
            old_sec_x = old_morph.morphology_data["neuron"].section_data[idx, 1] / 1e3
            coord_to_sec_id_x[coord[0], coord[1], coord[2]] = (old_sec_id, old_sec_x)

        old_to_new_sec_id = dict()
        old_sec_x_list = dict()
        new_sec_x_list = dict()
        # We just need to find the maximal old sec_x still present,
        # that value will map to sec_x 1.0 in new (stored as int sec_x*1000)

        new_dend_idx = np.where(new_morph.morphology_data["neuron"].section_data[:, 2] == 3)[0]
        for idx in new_dend_idx:
            coord = (new_morph.morphology_data["neuron"].geometry[idx, :3] * 1e9).astype(int)
            new_sec_id = new_morph.morphology_data["neuron"].section_data[idx, 0]
            new_sec_x = new_morph.morphology_data["neuron"].section_data[idx, 1] / 1e3

            try:
                old_sec_id, old_sec_x = coord_to_sec_id_x[coord[0], coord[1], coord[2]]
            except:
                import traceback
                print(traceback.format_exc())
                print(f"Coordinate point in new file missing in old file.\n"
                      f"Old morphology: {old_morph.swc_filename}\n"
                      f"New morphology: {new_morph.swc_filename}")
                import pdb
                pdb.set_trace()

            # TODO: What happens if a branch looses one of its side branches, will then those two sections
            #       be merged into one? This code will then complain.

            if old_sec_id in old_to_new_sec_id:
                assert old_to_new_sec_id[old_sec_id] == new_sec_id, \
                    (f"Old sec_id {old_sec_id} maps to multiple " 
                     f"new sec_id {new_sec_id} and {old_to_new_sec_id[old_sec_id]}")
            else:
               old_to_new_sec_id[old_sec_id] = new_sec_id

            if old_sec_id not in old_sec_x_list:
                old_sec_x_list[old_sec_id] = [old_sec_x]
                new_sec_x_list[old_sec_id] = [new_sec_x]
            else:
                old_sec_x_list[old_sec_id].append(old_sec_x)
                new_sec_x_list[old_sec_id].append(new_sec_x)

        new_to_old_sec_id = dict()
        for old_sec_id, new_sec_id in old_to_new_sec_id.items():
            if new_sec_id not in new_to_old_sec_id:
                new_to_old_sec_id[new_sec_id] = [old_sec_id]
            else:
                new_to_old_sec_id[new_sec_id].append(old_sec_id)

        for old_sec_id, new_sec_id in old_to_new_sec_id.items():
            # Is smallest section_x not zero?
            if old_sec_x_list[old_sec_id][0] > 0:
                # Are there multiple sections mapping to this new section_id?
                if len(new_to_old_sec_id[new_sec_id]) == 1 or old_sec_id == min(new_to_old_sec_id[new_sec_id]):
                    # Either just one section mapping to the new section_id
                    # or this is the smallest old section_id mapping to the new section_id
                    old_sec_x_list[old_sec_id].insert(0, 0)
                    new_sec_x_list[old_sec_id].insert(0, 0)

        # Soma ID now -1, updated mapping
        neuron_section_lookup = {-1: (-1, np.array([0, 1]), np.array([0, 1]))}  # Add SOMA mapping. SecX 0-1 --> ID 0-1

        for old_sec_id in old_to_new_sec_id.keys():

            assert (np.diff(old_sec_x_list[old_sec_id]) > 0).all()
            assert (np.diff(new_sec_x_list[old_sec_id]) > 0).all()

            neuron_section_lookup[old_sec_id] = (old_to_new_sec_id[old_sec_id],
                                                 np.array(old_sec_x_list[old_sec_id]),
                                                 np.array(new_sec_x_list[old_sec_id]))

        if True:
            assert len(neuron_section_lookup) > 3, (f"Section lookup has few elements. Does morphologies match?"
                                                     f"\nOld = {old_path, old_param_key, old_morph_key}"
                                                     f"\nNew = {new_path, new_param_key, new_morph_key}")

            # Check that all new_coords exist in the old_coords list.

        # import pdb
        # pdb.set_trace()

        self.section_lookup[old_param_key, old_morph_key, old_path] = neuron_section_lookup

    def remap_sections_helper(self, neuron_id, old_sec_id, old_sec_x):

        """
            Remaps sections for neuron_id

            This does currently NOT take AXON DEGENERATION into account

            Args:
                neuron_id (int) : Neuron ID
                old_sec_id (np.array) : Section ID of old section
                old_sec_x (np.array) : Section X of old section (float 0-1)

            Returns:
                keep_idx (np.array) : Indexes to keep
                new_sec_id (np.array) : Section ID of new section, positions are np.nan if not mapped
                new_sec_x (np.array) : Section X of new section, positions are np.nan if not mapped

        """

        keep_idx = np.zeros(old_sec_id.shape, dtype=bool)
        new_sec_id = np.full(old_sec_id.shape, -9999999, dtype=int)
        new_sec_x = np.full(old_sec_id.shape, np.nan)

        old_param_key, old_morph_key, old_path = self.find_old_morphology(neuron_id=neuron_id)

        lookup = self.section_lookup[old_param_key, old_morph_key, old_path]

        for idx, (old_id, old_x) in enumerate(zip(old_sec_id, old_sec_x)):

            if old_id in lookup:
                new_id = lookup[old_id][0]

                # We know range of old x that are still valid, and new x they map to. If outside, returns nan
                try:
                    new_x = np.interp(old_x, lookup[old_id][1], lookup[old_id][2], right=np.nan, left=np.nan)
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

                if not np.isnan(new_x):
                    assert 0 <= new_x <= 1, f"Out of range new_x={new_x} (required 0-1)"

                    keep_idx[idx] = True
                    new_sec_id[idx] = new_id
                    new_sec_x[idx] = new_x

        return np.where(keep_idx)[0], new_sec_id, new_sec_x


def cli():
    parser = ArgumentParser(
        description="Replace WT morphologies with PD morphologies, remove orphaned synapses and input",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument("original_network_path", type=str)
    parser.add_argument("new_network_path", type=str)
    parser.add_argument("original_snudda_path", type=str)
    parser.add_argument("new_snudda_path", type=str)
    parser.add_argument("--remap", action="store_true")

    args = parser.parse_args()

    if not os.path.isdir(args.new_network_path):
        os.mkdir(args.new_network_path)

    original_network_file = os.path.join(args.original_network_path, "network-synapses.hdf5")
    new_network_file = os.path.join(args.new_network_path, "network-synapses.hdf5")
    original_input_file = os.path.join(args.original_network_path, "input-spikes.hdf5")
    new_input_file = os.path.join(args.new_network_path, "input-spikes.hdf5")

    swap = SwapToDegeneratedMorphologies(original_network_file=original_network_file,
                                         new_network_file=new_network_file,
                                         original_snudda_data_dir=args.original_snudda_path,
                                         new_snudda_data_dir=args.new_snudda_path,
                                         original_input_file=original_input_file,
                                         new_input_file=new_input_file)
    swap.write_new_network_file()
    swap.write_new_input_file(remap_removed_input=args.remap)

    swap.close()


if __name__ == "__main__":
    cli()
