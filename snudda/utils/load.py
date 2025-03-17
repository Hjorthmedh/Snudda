#!/usr/bin/env python3
import os
import sys
import json
import timeit

import numpy as np

from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.numpy_encoder import NumpyEncoder
import scipy.sparse as sparse
from scipy.spatial import distance_matrix


class SnuddaLoad(object):
    """
    Load data from network-neuron-positions.hdf5 or network-neuron-synapses.hdf5 into python dictionary
    """

    ############################################################################

    def __init__(self, network_file, snudda_data=None, load_synapses=True, verbose=False):

        """
        Constructor

        Args:
            network_file (str) : Data file to load
            snudda_data (str, optional) : Snudda Data path, if you want to override the one specified in the hdf5 file
            load_synapses (bool, optional) : Whether to read synapses into memory, or keep them on disk (this keeps file open)
            verbose (bool, optional) : Print more info during execution

        """

        # This variable will only be set if the synapses are not kept in
        # memory so that we can access them later, otherwise the hdf5 file is
        # automatically closed
        self.hdf5_file = None
        self.verbose = verbose
        self.config = None

        self.network_file = None
        self.snudda_data = snudda_data

        if network_file:
            alt_file = os.path.join(network_file, "network-synapses.hdf5")

            if os.path.isdir(network_file):
                if not os.path.isfile(alt_file):
                    raise ValueError(f"Network path {network_file} specified, but no file {alt_file}")
                network_file = alt_file

            self.data = self.load_hdf5(network_file, load_synapses)
        else:
            self.data = None

    ############################################################################

    def close(self):

        """ Close hdf5 data file. """

        if self.hdf5_file:
            self.hdf5_file.close()
            self.hdf5_file = None

    ############################################################################

    def __del__(self):

        if self.hdf5_file is not None:
            try:
                self.hdf5_file.close()
            except:
                print("Unable to close HDF5, already closed?")

    ############################################################################

    @staticmethod
    def to_str(data_str):

        """ Helper function to convert data to string. """

        if data_str is None:
            return None

        if type(data_str) in [bytes, np.bytes_]:
            return data_str.decode()

        # Warn the user if they accidentally call to_str on an int or something else
        assert type(data_str) == str, f"to_str is used on strings or bytes (that are converted to str), " \
                                      f"received {type(data_str)} -- {data_str}"

        return data_str

    ############################################################################

    def load_hdf5(self, network_file, load_synapses=True, load_morph=False):

        """
        Load data from hdf5 file.

        Args:
            network_file (str) : Network file to load data from
            load_synapses (bool) : Load synapses into memory, or read on demand from file (keeps file open)
            load_morph

        Returns:
            data (dictionary) : Dictionary with data.

    Data format. The top level of the data hierarchy are "meta", "morphologies", "network".

        "meta" hierarchy:
            "slurm_id" (int) : Run ID of Slurm job
            "axon_stump_id_flag" (bool) : Should axon be replaced by a axon stump
            "config" : Config data
            "config_file" : Config file
            "connectivity_distributions" : Information about connections allowed, including pruning information
            "hyper_voxel_ids" : List of all hyper voxels IDs
            "hyper_voxel_size" (int, int, int): Number of hyper voxels along each dimension x, y, z
            "hyper_voxel_width" (float, float, float) : Size of hypervoxel in meters (x,y,z)
            "num_hyper_voxels" (int) : Number of hyper voxels
            "position_file" : Neuron position file
            "simulation_origo" (float, float, float) : Simulation origo (x, y, z) in SI-units (m)
            "voxel_size" (float) : Voxel size in meters

        "morphologies" heirarchy:
            List of neuron names, contains the location of the morphologies
            (either SWC file or directory with multiple SWC files)

        "network" hierarchy:
            "gap_junctions" : Gap junction data matrix (see below for format)
            "num_gap_junctions" (int) : Number of gap junctions
            "num_synapses" (int) : Number of synapses
            "neurons" : Neuron data structure (see below for format)
            "synapses" : Synapse data matrix (see below for format)

    Neuron data format:


    Synapse data format (column within parenthesis):

        0: source_cell_id, 1: dest_cell_id, 2: voxel_x, 3: voxel_y, 4: voxel_z,
        5: hyper_voxel_id, 6: channel_model_id,
        7: source_axon_soma_dist (not SI scaled 1e6, micrometers),
        8: dest_dend_soma_dist (not SI scalled 1e6, micrometers)
        9: dest_seg_id, 10: dest_seg_x (int 0 - 1000, SONATA wants float 0.0-1.0)
        11: conductance (int, not SI scaled 1e12, in pS)
        12: parameterID

    Note on parameterID:

        If there are n parameter sets for the particular synapse type, then
        the ID to use is parameterID % n, this way we can reuse connectivity
        if we add more synapse parameter sets later.

    Gap junction format (column in parentheses):

        0: source_cell_id, 1: dest_cell_id, 2: source_seg_id, 3: dest_seg_id,
        4: source_seg_x, 5: dest_seg_x, 6: voxel_x, 7: voxel_y, 8: voxel_z,
        9: hyper_voxel_id, 10: conductance (integer, in pS)

        """

        assert not load_morph, "load_hdf5: load_morph=True currently disabled, does not handle morph variations"

        self.network_file = network_file

        if self.verbose:
            print(f"Loading {network_file}")

        start_time = timeit.default_timer()
        data = dict([])

        # Save a reference to the name of the loaded network file
        data["network_file"] = self.network_file

        # Blender notebook has hdf5 library/header file mismatch, so importing this only where needed
        # This allows us to use fake_load.py in snudda.utils

        import h5py
        f = h5py.File(network_file, 'r')

        # We need to keep f open if load_synapses = False, using "with" would close file

        if "config" in f["meta"]:
            if self.verbose:
                print("Loading config data from HDF5")
            # data["config"] = SnuddaLoad.to_str(f["meta/config"][()])
            self.config = json.loads(f["meta/config"][()])
            data["config"] = self.config

        # Added so this code can also load the position file, which
        # does not have the network group yet
        if "network/synapses" in f:
            try:
                data["num_neurons"] = f["network/neurons/neuron_id"].shape[0]
                data["neuron_id"] = f["network/neurons/neuron_id"][()]
            except:
                if "neuronID" in f["network/neurons"]:
                    print(f"\nERROR: Old version of the network detected!\n"
                          f"Please regenerate your Snudda network: {self.network_file}\n\n")

                import traceback
                print(traceback.format_exc())

                print("Type \"exit()\" to exit.")

                import pdb
                pdb.set_trace()

            if "num_synapses" in f["network"]:
                if f["network/num_synapses"].shape == ():
                    data["num_synapses"] = f["network/num_synapses"][()]
                else:
                    data["num_synapses"] = f["network/num_synapses"][0]

                if data["num_synapses"] != f["network/synapses"].shape[0]:
                    print(f"Expected {data['num_synapses']} synapses, found {f['network/synapses'].shape[0]} synapse rows")
                    data["num_synapses"] = f["network/synapses"].shape[0]
            else:
                data["num_synapses"] = f["network/synapses"].shape[0]

            if "num_gap_junctions" in f["network"]:
                if f["network/num_gap_junctions"].shape == ():
                    data["num_gap_junctions"] = f["network/num_gap_junctions"][()]
                else:
                    data["num_gap_junctions"] = f["network/num_gap_junctions"][0]

                if data["num_gap_junctions"] != f["network/gap_junctions"].shape[0]:
                    print(f"Expected {data['num_gap_junctions']} gap junctions, "
                          f"found {f['network/gap_junctions'].shape[0]} gap junction rows")
                    data["num_gap_junctions"] = f["network/gap_junctions"].shape[0]
            else:
                data["num_gap_junctions"] = f["network/gap_junctions"].shape[0]

            if data["num_synapses"] > 100e6:
                print(f"Found {data['num_synapses']} synapses (too many!), not loading them into memory!")
                load_synapses = False

            if "network/hyper_voxel_ids" in f:
                data["hyper_voxel_ids"] = f["network/hyper_voxel_ids"][()]

            if "meta/hyper_voxel_size" in f:
                data["hyper_voxel_size"] = f["meta/hyper_voxel_size"][()]

            if "meta/hyper_voxel_width" in f:
                data["hyper_voxel_width"] = f["meta/hyper_voxel_width"][()]

            if "meta/num_hyper_voxels" in f:
                data["num_hyper_voxels"] = f["meta/num_hyper_voxels"][()]

            if load_synapses:
                # 0: source_cell_id, 1: dest_cell_id, 2: voxel_x, 3: voxel_y, 4: voxel_z,
                # 5: hyper_voxel_id, 6: channel_model_id,
                # 7: source_axon_soma_dist (not SI scaled 1e6, micrometers),
                # 8: dest_dend_soma_dist (not SI scalled 1e6, micrometers)
                # 9: dest_seg_id, 10: dest_seg_x (int 0 - 1000, SONATA wants float 0.0-1.0)
                # 11: conductance (int, not SI scaled 1e12, in pS)
                # 12: parameterID
                #
                # Note on parameterID:
                # If there are n parameter sets for the particular synapse type, then
                # the ID to use is parameterID % n, this way we can reuse connectivity
                # if we add more synapse parameter sets later.

                data["synapses"] = f["network/synapses"][:]
                data["gap_junctions"] = f["network/gap_junctions"][:]

                # !!! Convert from voxel idx to coordinates
                if f["network/synapses"].shape[0] > 0:
                    data["synapse_coords"] = f["network/synapses"][:, 2:5] * f["meta/voxel_size"][()] \
                                            + f["meta/simulation_origo"][()]
                else:
                    data["synapse_coords"] = np.zeros((3, 0))
            else:
                # Point the data structure to the synapses and gap junctions on file
                # This will be slower, and only work while the file is open
                data["synapses"] = f["network/synapses"]
                data["gap_junctions"] = f["network/gap_junctions"]

                # We need to keep f alive, since we did not load synapses into
                # the memory
                self.hdf5_file = f

        else:
            data["num_neurons"] = f["network/neurons/neuron_id"].shape[0]
            assert data["num_neurons"] == f["network/neurons/neuron_id"][-1] + 1, \
                "Internal error, something fishy with number of neurons found"

        data["config_file"] = SnuddaLoad.to_str(f["meta/config_file"][()])

        if "meta/position_file" in f:
            data["position_file"] = SnuddaLoad.to_str(f["meta/position_file"][()])

        if "meta/slurm_id" in f:
            if type(f["meta/slurm_id"][()]) in [bytes, np.bytes_]:
                data["slurm_id"] = int(f["meta/slurm_id"][()].decode())
            else:
                data["slurm_id"] = int(f["meta/slurm_id"][()])

        else:
            if self.verbose:
                print("No slurm_id set, using -1")
            data["slurm_id"] = -1

        if "meta/simulation_origo" in f:
            data["simulation_origo"] = f["meta/simulation_origo"][()]

        if "meta/voxel_size" in f:
            data["voxel_size"] = f["meta/voxel_size"][()]

        if "meta/axon_stump_id_flag" in f:
            data["axon_stump_id_flag"] = f["meta/axon_stump_id_flag"][()]

        if "meta/snudda_data" in f:
            data["snudda_data"] = SnuddaLoad.to_str(f["meta/snudda_data"][()])

            if self.snudda_data is None:
                self.snudda_data = data["snudda_data"]

        data["neurons"] = self.extract_neurons(f)

        # This is for old format, update for new format
        if "parameters" in f:
            # print("Parameters found, loading")
            data["synapse_range"] = f["parameters/synapse_range"][()]
            data["gap_junction_range"] = f["parameters/gap_junction_range"][()]
            data["min_synapse_spacing"] = f["parameters/min_synapse_spacing"][()]

        data["neuron_positions"] = f["network/neurons/position"][()]
        data["name"] = [SnuddaLoad.to_str(x) for x in f["network/neurons/name"][()]]

        if "population_unit_id" in f["network/neurons"]:
            data["population_unit"] = f["network/neurons/population_unit_id"][()]
        else:
            if self.verbose:
                print("No Population Units detected.")
            data["population_unit"] = np.zeros(data["num_neurons"], dtype=int)

        # TODO: Remove this, or make it able to handle multiple morphologies for each neuron_name,
        #  ie when morphologies is given as a dir
        if load_morph and "morphologies" in f:
            data["morph"] = dict([])

            for morph_name in f["morphologies"].keys():
                data["morph"][morph_name] = {"swc": f["morphologies"][morph_name]["swc"][()],
                                             "location": f["morphologies"][morph_name]["location"][()]}

        data["connectivity_distributions"] = dict([])

        if "connectivity_distributions" in f["meta"]:
            orig_connectivity_distributions = \
                json.loads(SnuddaLoad.to_str(f["meta/connectivity_distributions"][()]))

            for keys in orig_connectivity_distributions:
                (pre_type, post_type) = keys.split("$$")
                data["connectivity_distributions"][pre_type, post_type] \
                    = orig_connectivity_distributions[keys]

        if "synapses" in data and self.verbose:
            if "gap_junctions" in data:
                print(f"Loading {len(data['neurons'])} neurons with {data['num_synapses']} synapses"
                      f" and {data['num_gap_junctions']} gap junctions")
            else:
                print(f"Loading {len(data['neurons'])} neurons with {data['synapses'].shape[0]} synapses")

        if self.verbose:
            print(f"Load done. {timeit.default_timer() - start_time:.1f}")

        if load_synapses:
            f.close()
        else:
            self.hdf5_file = f

        return data

    ############################################################################

    @staticmethod
    def gather_extra_axons(hdf5_file):

        extra_axons = dict()

        if "extra_axons" in hdf5_file["network/neurons"]:
            for neuron_id, axon_name, position, rotation, swc_file \
                    in zip(hdf5_file["network/neurons/extra_axons/parent_neuron"],
                           hdf5_file["network/neurons/extra_axons/name"],
                           hdf5_file["network/neurons/extra_axons/position"],
                           hdf5_file["network/neurons/extra_axons/rotation"],
                           hdf5_file["network/neurons/extra_axons/morphology"]):

                if neuron_id not in extra_axons:
                    extra_axons[neuron_id] = dict()

                extra_axons[neuron_id][axon_name] = dict()
                extra_axons[neuron_id][axon_name]["position"] = position.copy()
                extra_axons[neuron_id][axon_name]["rotation"] = rotation.copy().reshape(3, 3)
                extra_axons[neuron_id][axon_name]["morphology"] = SnuddaLoad.to_str(swc_file)

        return extra_axons

    ############################################################################

    @staticmethod
    def extract_neurons(hdf5_file):

        """
        Helper function to extract neuron data from hdf5 file and put it in a dictionary.

        Args:
            hdf5_file : hdf5 file object

        Returns:
            List containing neurons as dictionary elements.

        """

        neurons = []
        extra_axons = SnuddaLoad.gather_extra_axons(hdf5_file=hdf5_file)

        for name, neuron_id, hoc, pos, rot, virtual, vID, \
            axon_density_type, axon_density, axon_density_radius, \
            axon_density_bounds_xyz, \
            morph, neuron_path, \
            parameter_key, morphology_key, modulation_key, population_unit_id, \
            reaction_diffusion_file \
                in zip(hdf5_file["network/neurons/name"][:],
                       hdf5_file["network/neurons/neuron_id"][:],
                       hdf5_file["network/neurons/hoc"][:],
                       hdf5_file["network/neurons/position"][()],
                       hdf5_file["network/neurons/rotation"][()],
                       hdf5_file["network/neurons/virtual_neuron"][:],
                       hdf5_file["network/neurons/volume_id"][:],
                       hdf5_file["network/neurons/axon_density_type"][:],
                       hdf5_file["network/neurons/axon_density"][:],
                       hdf5_file["network/neurons/axon_density_radius"][:],
                       hdf5_file["network/neurons/axon_density_bounds_xyz"][:],
                       hdf5_file["network/neurons/morphology"][:],
                       hdf5_file["network/neurons/neuron_path"][:],
                       hdf5_file["network/neurons/parameter_key"][:],
                       hdf5_file["network/neurons/morphology_key"][:],
                       hdf5_file["network/neurons/modulation_key"][:],
                       hdf5_file["network/neurons/population_unit_id"][:],
                       hdf5_file["network/neurons/reaction_diffusion_file"][:]
                       ):

            n = dict([])

            n["name"] = SnuddaLoad.to_str(name)

            if morph is not None:
                n["morphology"] = SnuddaLoad.to_str(morph)

            # Naming convention is TYPE_X, where XX is a number starting from 0
            n["type"] = n["name"].split("_")[0]

            n["neuron_id"] = neuron_id
            n["volume_id"] = SnuddaLoad.to_str(vID)
            n["hoc"] = SnuddaLoad.to_str(hoc)
            n["neuron_path"] = SnuddaLoad.to_str(neuron_path)

            n["position"] = pos.copy()
            n["rotation"] = rot.copy().reshape(3, 3)
            n["virtual_neuron"] = virtual

            if len(axon_density_type) > 0:
                n["axon_density_type"] = SnuddaLoad.to_str(axon_density_type)
            else:
                n["axon_density_type"] = None

            if len(axon_density) > 0:
                n["axon_density"] = SnuddaLoad.to_str(axon_density)
            else:
                n["axon_density"] = None

            if n["axon_density_type"] == "xyz":
                n["axon_density_bounds_xyz"] = axon_density_bounds_xyz
            else:
                n["axon_density_bounds_xyz"] = None

            n["axon_density_radius"] = axon_density_radius

            # If the code fails here, use snudda/utils/upgrade_old_network_file.py to upgrade your old data files
            par_key = SnuddaLoad.to_str(parameter_key)
            morph_key = SnuddaLoad.to_str(morphology_key)
            mod_key = SnuddaLoad.to_str(modulation_key)
            n["parameter_key"] = par_key if len(par_key) > 0 else None
            n["morphology_key"] = morph_key if len(morph_key) > 0 else None
            n["modulation_key"] = mod_key if len(mod_key) > 0 else None
            n["reaction_diffusion_file"] = SnuddaLoad.to_str(reaction_diffusion_file) if len(reaction_diffusion_file) > 0 else None

            n["population_unit"] = population_unit_id

            if neuron_id in extra_axons:
                n["extra_axons"] = extra_axons[neuron_id].copy()

            neurons.append(n)

        return neurons

    ############################################################################

    def load_config_file(self):

        """ Load config data from JSON file. """

        if self.config is None:
            config_file = self.data["config_file"]
            self.config = json.load(open(config_file, 'r'))

    ############################################################################

    def synapse_iterator(self, chunk_size=1000000, data_type=None):

        """
        Iterates through all synapses in chunks (default 1e6 synapses).

        Args:
            chunk_size (int) : Number of synapses per chunk
            data_type (string) : "synapses" (default) or "gap_junctions"

        Returns:
            Iterator over the synapses
        """

        if data_type is None:
            data_type = "synapses"

        # data_type is "synapses" or "gap_junctions"
        assert data_type in ["synapses", "gap_junctions"]

        num_rows = self.data[data_type].shape[0]
        if num_rows == 0:
            # No synapses
            return

        chunk_size = min(num_rows, chunk_size)
        num_steps = int(np.ceil(num_rows / chunk_size))

        row_start = 0

        for row_end in np.linspace(chunk_size, num_rows, num_steps, dtype=int):
            synapses = self.data[data_type][row_start:row_end, :]
            row_start = row_end

            yield synapses

    ############################################################################

    def gap_junction_iterator(self, chunk_size=1000000):

        """
        Iterates through all gap junctions in chunks (default 1e6 gap junctions).

        Args:
            chunk_size (int) : Number of gap junctions per chunk

        Returns:
            Iterator over the gap junctions
        """

        return self.synapse_iterator(chunk_size=chunk_size, data_type="gap_junctions")

    ############################################################################

    # Helper methods for sorting

    @staticmethod
    def _row_eval_post_pre(row, num_neurons):
        return row[1] * num_neurons + row[0]

    @staticmethod
    def _row_eval_post(row, num_neurons):
        return row[1]

    ############################################################################

    def find_synapses_slow(self, pre_id, n_max=1000000):

        """
        Returns subset of synapses.

        Args:
            pre_id (int) : Pre-synaptic neuron ID (can also be a list)
            n_max (int) : Maximum number of synapses to return

        Returns:
            Subset of synapse matrix, synapse coordinates

        """

        if self.verbose:
            print(f"Finding synapses originating from {pre_id}, this is slow")

        synapses = np.zeros((n_max, 13), dtype=np.int32)
        syn_ctr = 0

        if np.issubdtype(type(pre_id), np.integer):
            for syn_list in self.synapse_iterator():
                for syn in syn_list:
                    if syn[0] == pre_id:
                        synapses[syn_ctr, :] = syn
                        syn_ctr += 1

        else:
            for syn_list in self.synapse_iterator():
                for syn in syn_list:
                    if syn[0] in pre_id:
                        synapses[syn_ctr, :] = syn
                        syn_ctr += 1

        synapse_coords = synapses[:, 2:5][:syn_ctr, :] * self.data["voxel_size"] + self.data["simulation_origo"]

        return synapses[:syn_ctr, :], synapse_coords

    ############################################################################

    # Either give preID and postID, or just postID

    def find_synapses(self, pre_id=None, post_id=None, silent=True, return_index=False):

        """
        Returns subset of synapses.

        Args:
            pre_id (int) : Pre-synaptic neuron ID
            post_id (int) : Post-synaptic neuron ID
            silent (bool) : Work quietly or verbosely

        Returns:
            Subset of synapse matrix.

        """

        if self.data["synapses"].shape[0] == 0:
            if not silent:
                print("No synapses in network")
            if return_index:
                return None, None, None
            else:
                return None, None

        if post_id is None:
            assert return_index is False, "You must specify pre_id and post_id if return_index is True"
            return self.find_synapses_slow(pre_id=pre_id)

        assert post_id is not None, "Must specify at least postID"

        num_rows = self.data["synapses"].shape[0]
        num_neurons = len(self.data["neuron_id"])

        if pre_id is None:
            row_eval = self._row_eval_post
            val_target = post_id
        else:
            row_eval = self._row_eval_post_pre
            val_target = post_id * num_neurons + pre_id

        idx_a1 = 0
        idx_a2 = num_rows - 1

        idx_found = None

        # We use idxA1 and idxA2 as upper and lower range within which we
        # hope to find one of the synapses. Once we found a synapse row
        # we go up and down in matrix to find the range of the synapses
        # matching the requested condition. This works because the matrix is
        # sorted on postID, and then preID if postID matches

        if row_eval(self.data["synapses"][idx_a1, :], num_neurons) == val_target:
            idx_found = idx_a1

        if row_eval(self.data["synapses"][idx_a2, :], num_neurons) == val_target:
            idx_found = idx_a2

        # -1 since if idx_a1 and idx_a2 are one apart, we have checked all values
        while idx_a1 < idx_a2 - 1 and idx_found is None:

            idx_next = int(np.round((idx_a1 + idx_a2) / 2))
            val_next = row_eval(self.data["synapses"][idx_next, :], num_neurons)

            if val_next < val_target:
                idx_a1 = idx_next
            elif val_next > val_target:
                idx_a2 = idx_next
            else:
                # We found a hit
                idx_found = idx_next
                break

        if idx_found is None:
            # No synapses found
            if self.verbose:
                print("No synapses found")

            if return_index:
                return None, None, None
            else:
                return None, None

        # Find start of synapse range
        idx_b1 = idx_found
        val_b1 = row_eval(self.data["synapses"][idx_b1 - 1, :], num_neurons)

        while val_b1 == val_target and idx_b1 > 0:
            idx_b1 -= 1
            val_b1 = row_eval(self.data["synapses"][idx_b1 - 1, :], num_neurons)

        # Find end of synapse range
        idx_b2 = idx_found

        if idx_b2 + 1 < self.data["synapses"].shape[0]:
            val_b2 = row_eval(self.data["synapses"][idx_b2 + 1, :], num_neurons)

            while val_b2 == val_target and idx_b2 + 1 < self.data["synapses"].shape[0]:
                idx_b2 += 1
                val_b2 = row_eval(self.data["synapses"][idx_b2 + 1, :], num_neurons)

        synapses = self.data["synapses"][idx_b1:idx_b2 + 1, :].copy()

        if not silent and self.verbose:
            print(f"Synapse range, first {idx_b1}, last {idx_b2}")
            print(f"{synapses}")

        # Calculate coordinates
        synapse_coords = synapses[:, 2:5] * self.data["voxel_size"] + self.data["simulation_origo"]

        if return_index:
            return synapses, synapse_coords, np.arange(idx_b1, idx_b2+1)
        else:
            return synapses, synapse_coords

    ############################################################################

    def get_neuron_population_units(self, neuron_id=None, return_set=False):

        if neuron_id is not None:
            neuron_population_units = self.data["population_unit"][neuron_id].flatten().copy()
        else:
            neuron_population_units = self.data["population_unit"].flatten().copy()

        if return_set:
            return set(neuron_population_units)
        else:
            return neuron_population_units

    def get_neuron_types(self, neuron_id=None, return_set=False):

        if neuron_id is not None:
            neuron_types = [self.data["neurons"][x]["type"] for x in neuron_id]
        else:
            neuron_types = [x["type"] for x in self.data["neurons"]]

        if return_set:
            return set(neuron_types)
        else:
            return neuron_types

    ############################################################################

    # Returns neuron_id of all neurons of neuron_type
    # OBS, random_permute is not using a controled rng, so not affected by random seed set

    def get_neuron_id_of_type(self, neuron_type, num_neurons=None, random_permute=False, volume=None,
                              include_virtual=True):

        """
        Find all neuron ID of a specific neuron type.

        Args:
            neuron_type (string) : Neuron type (e.g. "FS")
            num_neurons (int) : Maximum number of neurons to return
            random_permute (bool) : Shuffle the resulting neuron IDs?
            volume (string) : volume_id containing neurons (default None -- all neurons of type)

        Returns:
            List of neuron ID of specified neuron type

        """

        neuron_id = np.array([x["neuron_id"] for x in self.data["neurons"]
                              if (neuron_type is None or x["type"] == neuron_type)
                              and (volume is None or x["volume_id"] == volume)
                              and (include_virtual or not x["virtual_neuron"])])

        assert not random_permute or num_neurons is not None, "random_permute is only valid when num_neurons is given"

        if num_neurons is not None:
            if random_permute:
                # Do not use this if you have a simulation with multiple
                # workers... they might randomize differently, and you might
                # fewer neurons in total than you wanted
                keep_idx = np.random.permutation(len(neuron_id))

                if len(keep_idx) > num_neurons:
                    keep_idx = keep_idx[:num_neurons]

                neuron_id = neuron_id[keep_idx]
            else:
                neuron_id = neuron_id[:num_neurons]

            if len(neuron_id) < num_neurons:
                if self.verbose:
                    print(f"get_neuron_id_of_type: wanted {num_neurons} only got {len(neuron_id)} "
                          f"neurons of type {neuron_type}")

        # Double check that all of the same type (or neuron_type is None)
        assert neuron_type is None or np.array([self.data["neurons"][x]["type"] == neuron_type for x in neuron_id]).all()
        assert volume is None or np.array([self.data["neurons"][x]["volume_id"] == volume for x in neuron_id]).all()

        return neuron_id

    def get_neuron_id(self, include_virtual=True):

        """ Returns all neuron_id, if include_virtual is set (default) virtual neurons are also included."""

        neuron_id = np.array([x["neuron_id"] for x in self.data["neurons"] if include_virtual or not x["virtual_neuron"]])

        return neuron_id

    def get_neuron_id_with_name(self, neuron_name, include_virtual=True):

        """
        Find neuron ID of neurons with a given name.

        Args:
            neuron_name (str): Name of neurons (e.g. "dSPN_0")
            include_virtual (bool): Should virtual neurons also be included?

        Returns:
            List of neuron ID
        """

        neuron_id = np.array([x["neuron_id"] for x in self.data["neurons"]
                              if x["name"] == neuron_name and (include_virtual or not x["virtual_neuron"])])

        return neuron_id

    def get_population_unit_members(self, population_unit, num_neurons=None, random_permute=False):

        """
        Returns neuron ID of neurons belonging to a specific population unit.

        Args:
            population_unit (int) : Population unit ID
            num_neurons (int) : Number of neurons to return (None = all hits)
            random_permute (bool) : Randomly shuffle neuron IDs?

        Returns:
            List of neuron ID belonging to population unit.

        """

        neuron_id = np.where(self.data["population_unit"] == population_unit)[0]

        if num_neurons:
            if random_permute:
                neuron_id = np.random.permutation(neuron_id)

            if len(neuron_id) > num_neurons:
                neuron_id = neuron_id[:num_neurons]

        # Just double check
        assert (self.data["population_unit"][neuron_id] == population_unit).all()

        return neuron_id

    def load_neuron(self, neuron_id):

        """
        Loads a specific neuron. Returns a NeuronMorphology object.

        Args:
            neuron_id (int): Neuron ID

        Returns:
            NeuronMorphology object.

        """

        neuron_prototype = NeuronPrototype(neuron_path=self.data["neurons"][neuron_id]["neuron_path"],
                                           snudda_data=self.data["snudda_data"],
                                           neuron_name=self.data["neurons"][neuron_id]["name"])
        neuron_object = neuron_prototype.clone(parameter_key=self.data["neurons"][neuron_id]["parameter_key"],
                                               morphology_key=self.data["neurons"][neuron_id]["morphology_key"],
                                               modulation_key=self.data["neurons"][neuron_id]["modulation_key"],
                                               position=self.data["neurons"][neuron_id]["position"],
                                               rotation=self.data["neurons"][neuron_id]["rotation"])
        return neuron_object

    def iter_neuron_id(self):

        for x, v in enumerate(self.data["neurons"]):
            assert x == v["neuron_id"], \
                f"Neuron at position {x} has neuron_id {v['neuron_id']} (should be same)"
            yield x

    def get_neuron_keys(self, neuron_id):
        n = self.data["neurons"][neuron_id]
        return n["parameter_key"], n["morphology_key"], n["modulation_key"]

    def get_neuron_params(self, neuron_id):

        neuron_path = self.data["neurons"][neuron_id]["neuron_path"]
        parameter_key = self.data["neurons"][neuron_id]["parameter_key"]
        parameter_file = os.path.join(neuron_path, "parameters.json")
        with open(parameter_file, "r") as f:
            parameter_data = json.load(f)

        param_data = parameter_data[parameter_key]

        if "modulation_key" in self.data["neurons"][neuron_id]:
            modulation_key = self.data["neurons"][neuron_id]["modulation_key"]
            modulation_file = os.path.join(neuron_path, "modulation.json")
            with open(modulation_file, "r") as f:
                modulation_data = json.load(f)
            mod_data = modulation_data[modulation_key]
        else:
            mod_data = None

        return param_data, mod_data

    def find_gap_junctions(self, neuron_id, n_max=1000000, return_index=False):

        """ Find gap junctions associated with neuron_id

        Args:
            neuron_id (int) : Neuron with gap junction (can also be a list)
            n_max (int) : Maximum number of gap junctions to return
            return_index (bool): Should third return value index be present

        Returns:
            Subset of gap junction matrix, gap junction coordinates

        """

        if self.verbose:
            print(f"Finding gap junctions connecting neuron {neuron_id}")

        gap_junctions = np.zeros((n_max, 11), dtype=np.int32)
        gj_ctr = 0
        gj_index = 0
        gj_index_list = np.zeros((n_max,), dtype=int)

        if np.issubdtype(type(neuron_id), np.integer):
            for gj_list in self.gap_junction_iterator():
                for gj in gj_list:
                    if gj[0] == neuron_id or gj[1] == neuron_id:
                        gap_junctions[gj_ctr, :] = gj
                        gj_ctr += 1
                        gj_index_list[gj_ctr] = gj_index
                    gj_index += 1
        else:
            for gj_list in self.gap_junction_iterator():
                for gj in gj_list:
                    if gj[0] in neuron_id or gj[1] in neuron_id:
                        gap_junctions[gj_ctr, :] = gj
                        gj_ctr += 1
                        gj_index_list[gj_ctr] = gj_index
                    gj_index += 1

        gj_coords = gap_junctions[:, 6:9][:gj_ctr, :] * self.data["voxel_size"] + self.data["simulation_origo"]

        if return_index:
            return gap_junctions[:gj_ctr, :], gj_coords, gj_index_list[:gj_ctr]
        else:
            return gap_junctions[:gj_ctr, :], gj_coords

    def get_centre_neurons_iterator(self, n_neurons=None, neuron_type=None, centre_point=None, max_distance=None,
                                    return_distance=True, include_virtual=False):

        """ Return neuron id:s, starting from the centre most and moving outwards

        Args:
            n_neurons (int) : Number of neurons to return, None = all available
            neuron_type (str) : Type of neurons to return, None = all available
            centre_point (np.array) : x,y,z of centre position, None = auto detect centre
        """

        if centre_point is None:
            centre_point = np.mean(self.data["neuron_positions"], axis=0)

        dist_to_centre = np.linalg.norm(self.data["neuron_positions"] - centre_point, axis=-1)
        idx = np.argsort(dist_to_centre)

        neuron_ctr = 0

        for neuron_id in idx:
            if not include_virtual and self.data["neurons"][neuron_id]["virtual_neuron"]:
                # Ignore virtual neurons
                continue

            if neuron_type is not None and self.data["neurons"][neuron_id]["type"] != neuron_type:
                # Wrong neuron type
                continue

            if max_distance is not None and dist_to_centre[neuron_id] > max_distance:
                # Stop iterator if max distance is reached
                return

            if return_distance:
                yield neuron_id, dist_to_centre[neuron_id]
            else:
                yield neuron_id

            neuron_ctr += 1

            if n_neurons is not None and neuron_ctr >= n_neurons:
                # Stop iterator if n_neurons are delivered
                return

    ############################################################################

    def create_connection_matrix(self, sparse_matrix=True, connection_type="synapses"):

        if sparse_matrix:
            connection_matrix = sparse.lil_matrix((self.data["num_neurons"], self.data["num_neurons"]), dtype=np.int16)
        else:
            connection_matrix = np.zeros((self.data["num_neurons"], self.data["num_neurons"]), dtype=np.int16)

        for syn_row in self.data[connection_type]:
            connection_matrix[syn_row[0], syn_row[1]] += 1

        return connection_matrix

    def find_neighbours(self, neuron_id, connection_matrix=None, exclude_parent=True):

        if connection_matrix is None:
            connection_matrix = self.create_connection_matrix(sparse_matrix=False, connection_type="synapses")

        pre_neighbours = set(np.where(np.sum(connection_matrix[:, list(neuron_id)], axis=1))[0])
        post_neighbours = set(np.where(np.sum(connection_matrix[list(neuron_id), :], axis=0))[0])

        if exclude_parent:
            parent_id = set(neuron_id)
            pre_neighbours -= parent_id
            post_neighbours -= parent_id

        return pre_neighbours, post_neighbours

    def find_neighbours_gap_junctions(self, neuron_id, connection_matrix=None, exclude_parent=True):

        if connection_matrix is None:
            connection_matrix = self.create_connection_matrix(sparse_matrix=False, connection_type="gap_junctions")

        # This matrix should be symmetric!
        pre_neighbours = set(np.where(np.sum(connection_matrix[:, list(neuron_id)], axis=1))[0])
        post_neighbours = set(np.where(np.sum(connection_matrix[list(neuron_id), :], axis=0))[0])

        neighbours = pre_neighbours | post_neighbours

        if exclude_parent:
            parent_id = set(neuron_id)
            neighbours -= parent_id

        return neighbours

    def create_distance_matrix(self, neuron_id=None, pre_id=None, post_id=None):

        if neuron_id is not None and pre_id is not None and post_id is not None:
            raise ValueError("Specify either neuron_id or the two parameters pre_id and post_id.")

        if (pre_id is None) ^ (post_id is None):
            raise ValueError("pre_id and post_id must both either be specified, or neither")

        pos = self.data["neuron_positions"]

        if neuron_id is not None:
            pos = pos[neuron_id, :]
            dist_matrix = distance_matrix(pos, pos)

        elif pre_id is not None and post_id is not None:
            dist_matrix = distance_matrix(pos[pre_id, :], pos[post_id, :])

        else:
            dist_matrix = distance_matrix(pos, pos)

        return dist_matrix

    def print_all_synapse_counts_per_type(self):

        synapse_types = sorted(list(self.get_neuron_types(return_set=True)))
        connection_matrix = self.create_connection_matrix()

        for pre_type in synapse_types:
            for post_type in synapse_types:
                count = self.count_synapses_per_type(pre_neuron_type=pre_type,
                                                     post_neuron_type=post_type,
                                                     connection_matrix=connection_matrix)
                if count > 0:
                    print(f"{pre_type} -> {post_type}: {count} synapses")

    def count_synapses_per_type(self, pre_neuron_type=None, post_neuron_type=None, connection_matrix=None):

        """
            Args:
                pre_neuron_type (list): List of neuron types
                post_neuron_type (list): List of neuron types, e.g. ["dSPN", "iSPN"]
                connection_matrix

        """

        if connection_matrix is None:
            connection_matrix = self.create_connection_matrix()

        pre_mask = np.zeros((connection_matrix.shape[0]), dtype=bool)
        post_mask = np.zeros((connection_matrix.shape[1]), dtype=bool)

        if pre_neuron_type is None:
            pre_neuron_type = self.get_neuron_types(return_set=True)
        elif type(pre_neuron_type) != list:
            pre_neuron_type = [pre_neuron_type]

        if post_neuron_type is None:
            post_neuron_type = self.get_neuron_types(return_set=True)
        elif type(post_neuron_type) != list:
            post_neuron_type = [post_neuron_type]

        for neuron_type in pre_neuron_type:
            pre_id = self.get_neuron_id_of_type(neuron_type)
            pre_mask[pre_id] = True

        for neuron_type in post_neuron_type:
            post_id = self.get_neuron_id_of_type(neuron_type)
            post_mask[post_id] = True

        synapse_count = np.sum(connection_matrix[pre_id, :][:, post_id])

        return synapse_count

    def count_incoming_connections(self, neuron_type):

        neuron_id = self.get_neuron_id_of_type(neuron_type)
        neuron_id_mask = np.zeros((self.data["num_neurons"],), dtype=bool)
        neuron_id_mask[neuron_id] = True

        synapse_count = 0
        gap_junction_count = 0

        for synapses in self.synapse_iterator():
            synapse_count += np.sum(neuron_id_mask[synapses[:, 1]])

        for gap_junctions in self.gap_junction_iterator():
            gap_junction_count += np.sum(neuron_id_mask[gap_junctions[:, 0]])
            gap_junction_count += np.sum(neuron_id_mask[gap_junctions[:, 1]])

        gap_junction_count /= 2

        return synapse_count, gap_junction_count


def snudda_load_cli():
    """ Command line parser for SnuddaLoad script """

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Load snudda network file (hdf5)")
    parser.add_argument("network_file", help="Network file (hdf5)", type=str)
    parser.add_argument("--listN", help="Lists neurons in network", action="store_true")
    parser.add_argument("--listT", type=str, help="List neurons of type, --listT ? list the types.", default=None)
    parser.add_argument("--listPre", help="List pre synaptic neurons", type=int)
    parser.add_argument("--listPost", help="List post synaptic neurons (slow)", type=int)
    parser.add_argument("--listGJ", help="List gap junctions (slow)", type=int)
    parser.add_argument("--listTotalIncoming", help="List number of total incoming connections to neuron type",
                        type=str, default=None)
    parser.add_argument("--keepOpen", help="This prevents loading of synapses to memory, and keeps HDF5 file open",
                        action="store_true")
    parser.add_argument("--detailed", help="More information", action="store_true")
    parser.add_argument("--voxels", help="Voxel information", action="store_true")
    parser.add_argument("--centre", help="List n neurons in centre (-1 = all)", type=int)
    parser.add_argument("--listParam", help="List parameters for neuron_id", type=int)
    parser.add_argument("--countSyn", help="Count synapses per type", action="store_true")

    args = parser.parse_args()

    if args.keepOpen:
        load_synapses = False
    else:
        load_synapses = True

    nl = SnuddaLoad(args.network_file, load_synapses=load_synapses, verbose=True)

    if args.listN:
        print("Neurons in network: ")

        if args.detailed:
            for nid, name, pos, par_key, morph_key, mod_key, neuron_path, pop_id, virt_flag \
                    in [(x["neuron_id"], x["name"], x["position"],
                         x["parameter_key"], x["morphology_key"], x["modulation_key"],
                         x["neuron_path"], x["population_unit"], x["virtual_neuron"])
                        for x in nl.data["neurons"]]:
                print(f"{nid} : {name}{' [virtual]' if virt_flag else ''}, ({pos[0]:.6f}, {pos[1]:.6f}, {pos[2]:.6f}) "
                      f"pop_id {pop_id}, par_key {par_key}, morph_key {morph_key}, neuron_path: {neuron_path}")
        else:
            for nid, name, pos, pid, virt_flag \
                in [(x["neuron_id"], x["name"], x["position"], x["population_unit"], x["virtual_neuron"])
                                        for x in nl.data["neurons"]]:
                print(f"{nid} {'V' if virt_flag else ':'} {name} [{pid}], ({pos[0]:.6f}, {pos[1]:.6f}, {pos[2]:.6f})")

    if args.listT is not None:
        if args.listT == "?":
            print("List neuron types in network:")

            n_types = np.unique([x["type"] for x in nl.data["neurons"]])
            for nt in n_types:
                num = len([x["type"] for x in nl.data["neurons"] if x["type"] == nt])
                print(f"{nt} ({num} total)")

        else:
            print(f"Neurons of type {args.listT}:")
            n_of_type = [(x["neuron_id"], x["name"]) for x in nl.data["neurons"]
                         if x["type"] == args.listT]
            for nid, name in n_of_type:
                print("%d : %s" % (nid, name))

    if args.listPre is not None:
        print(f"List neurons pre-synaptic to neuron_id = {args.listPre} "
              f"({nl.data['neurons'][args.listPre]['name']})")
        synapses, synapse_coords = nl.find_synapses(post_id=args.listPre)

        if synapses is None:
            print("No pre synaptic neurons were found.")
        else:
            print(f"The neuron receives {synapses.shape[0]} synapses")
            pre_id = np.unique(synapses[:, 0])

            for nid, name in [(x["neuron_id"], x["name"]) for x in nl.data["neurons"] if x["neuron_id"] in pre_id]:
                n_syn = np.sum(synapses[:, 0] == nid)
                print(f"{nid} : {name} ({n_syn} synapses)")

                if args.detailed:
                    idx = np.where(synapses[:, 0] == nid)[0]
                    for i in idx:
                        print(f" -- SecID {synapses[i, 9]}, SecX {synapses[i, 10] * 1e-3:.2f}, "
                              f"Soma dist: {synapses[i,8]:.1f} μm, "
                              f"Coord: ({synapse_coords[i, 0]*1e6:.1f}, "
                              f"{synapse_coords[i, 1]*1e6:.1f}, "
                              f"{synapse_coords[i, 2]*1e6:.1f}) μm, "
                              f"Cond: {synapses[i, 11] * 1e-3:.2f} nS")

                    print("")

                if args.voxels:
                    idx = np.where(synapses[:, 0] == nid)[0]
                    for i in idx:
                        print(f" -- voxels: {synapses[i, 2:5]}, hyper id: {synapses[i, 5]}, "
                              f" sec id {synapses[i, 9]}, sec x {synapses[i, 10] * 1e-3:.2f}, "
                              f"Soma dist: {synapses[i,8]:.1f} μm, ")
                    print("")

    if args.listPost is not None:
        print(f"List neurons post-synaptic to neuron_id = {args.listPost}"
              f" ({nl.data['neurons'][args.listPost]['name']}):")
        synapses, synapse_coords = nl.find_synapses(pre_id=args.listPost)
        print(f"The neuron makes {synapses.shape[0]} synapses on other neurons")

        if synapses is None:
            print("No post synaptic targets found.")
        else:
            post_id = np.unique(synapses[:, 1])

            for nid, name in [(x["neuron_id"], x["name"]) for x in nl.data["neurons"] if x["neuron_id"] in post_id]:
                n_syn = np.sum(synapses[:, 1] == nid)
                print(f"{nid} : {name} ({n_syn} synapses)")

                if args.detailed:
                    idx = np.where(synapses[:, 1] == nid)[0]
                    for i in idx:
                        print(f" -- SecID {synapses[i, 9]}, SecX {synapses[i, 10] * 1e-3:.2f}, "
                              f"Soma dist: {synapses[i,8]:.1f} μm, "
                              f"Coord: ({synapse_coords[i, 0]*1e6:.1f}, "
                              f"{synapse_coords[i, 1]*1e6:.1f}, "
                              f"{synapse_coords[i, 2]*1e6:.1f}) μm, "
                              f"Cond: {synapses[i, 11] * 1e-3:.2f} nS")

                    print("")

    if args.listGJ is not None:
        print(f"List gap junctions of neuron_id = {args.listGJ}"
              f" ({nl.data['neurons'][args.listGJ]['name']})")
        gap_junctions, gap_junction_coords = nl.find_gap_junctions(neuron_id=args.listGJ)

        if gap_junctions.shape[0] == 0:
            print("No gap junctions on neuron.")
        else:
            connected_id = set(gap_junctions[:, 0]).union(gap_junctions[:, 1])
            connected_id.remove(args.listGJ)

            for nid, name in [(x["neuron_id"], x["name"]) for x in nl.data["neurons"] if x["neuron_id"] in connected_id]:
                n_gj = np.sum(gap_junctions[:, 0] == nid) + np.sum(gap_junctions[:, 1] == nid)
                print(f"{nid} : {name} ({n_gj} gap junctions)")

                if args.detailed:
                    idx1 = np.where(gap_junctions[:, 0] == nid)[0]
                    for i in idx1:
                        print(f" -- SecID {gap_junctions[i, 2]}, SecX {gap_junctions[i, 4] * 1e-3:.3f}, "
                              f"Coord: ({gap_junction_coords[i, 0]*1e6:.1f}, "
                              f"{gap_junction_coords[i, 1]*1e6:.1f}, "
                              f"{gap_junction_coords[i, 2]*1e6:.1f}) μm, "
                              f"Cond: {gap_junctions[i, 10] * 1e-3:.3f} nS")

                    idx2 = np.where(gap_junctions[:, 1] == nid)[0]
                    for i in idx2:
                        print(f" -- SecID {gap_junctions[i, 3]}, SecX {gap_junctions[i, 5] * 1e-3:.3f}, "
                              f"Coord: ({gap_junction_coords[i, 0]*1e6:.1f}, "
                              f"{gap_junction_coords[i, 1]*1e6:.1f}, "
                              f"{gap_junction_coords[i, 2]*1e6:.1f}) μm, "
                              f"Cond: {gap_junctions[i, 10] * 1e-3:.3f} nS")

    if args.listParam is not None:
        param_data, mod_data = nl.get_neuron_params(neuron_id=args.listParam)
        print(f"Neuron ID {args.listParam}\n"
              f"parameters: {json.dumps(param_data, indent=4, cls=NumpyEncoder)}\n\n"
              f"modulation: {json.dumps(mod_data, indent=4, cls=NumpyEncoder)}\n")

    if args.centre:
        if args.centre < 0:
            n_neurons = None
        else:
            n_neurons = args.centre

        for neuron_id, centre_dist in nl.get_centre_neurons_iterator(n_neurons):
            pos = nl.data["neurons"][neuron_id]["position"] * 1e6
            name = nl.data["neurons"][neuron_id]["name"]
            print(f"{neuron_id} {name}: ({pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f}) μm, distance to centre {centre_dist*1e6:.1f} μm")

    if args.countSyn:
        nl.print_all_synapse_counts_per_type()

    if args.listTotalIncoming:
        incoming_to_type = args.listTotalIncoming

        synapse_count, gap_junction_count = nl.count_incoming_connections(neuron_type=incoming_to_type)
        print(f"All neurons of type {incoming_to_type} receive in total {synapse_count:.0f} synapses, "
              f"and have {gap_junction_count:.0f} gap junctions in total.")


if __name__ == "__main__":
    snudda_load_cli()
