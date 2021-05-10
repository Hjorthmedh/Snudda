import numpy as np
import timeit
import json
import os


class SnuddaLoad(object):

    ############################################################################

    def __init__(self, network_file, load_synapses=True, verbose=False):

        # This variable will only be set if the synapses are not kept in
        # memory so we can access them later, otherwise the hdf5 file is
        # automatically closed
        self.hdf5_file = None
        self.verbose = verbose
        self.config = None

        self.network_file = None

        if network_file:
            self.data = self.load_hdf5(network_file, load_synapses)
        else:
            self.data = None

    ############################################################################

    def close(self):
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
        if type(data_str) in [bytes, np.bytes_]:
            return data_str.decode()

        # Warn the user if they accidentally call to_str on an int or something else
        assert type(data_str) == str, "to_str is used on strings or bytes (that are converted to str)"

        return data_str

    ############################################################################

    def load_hdf5(self, network_file, load_synapses=True, load_morph=True):

        self.network_file = network_file

        if self.verbose:
            print(f"Loading {network_file}")

        start_time = timeit.default_timer()
        data = dict([])

        # Blender notebook has hdf5 library/header file mismatch, so importing this only where needed
        # This allows us to use fake_load.py in snudda.utils

        import h5py
        f = h5py.File(network_file, 'r')

        # We need to keep f open if load_synapses = False, using "with" would close file

        if "config" in f:
            if self.verbose:
                print("Loading config data from HDF5")
            data["config"] = SnuddaLoad.to_str(f["config"][()])
            self.config = json.loads(f["config"][()])

        # Added so this code can also load the position file, which
        # does not have the network group yet
        if "network/synapses" in f:
            data["nNeurons"] = f["network/neurons/neuronID"].shape[0]
            data["neuronID"] = f["network/neurons/neuronID"][()]

            try:
                data["nSynapses"] = f["network/nSynapses"][0]
            except:
                data["nSynapses"] = f["network/synapses"].shape[0]

            try:
                data["nGapJunctions"] = f["network/nGapJunctions"][0]
            except:
                data["nGapJunctions"] = f["network/gapJunctions"].shape[0]

            if data["nSynapses"] > 100e6:
                print(str(data["nSynapses"]) + " synapses, which is a lot, not loading them into memory!")
                load_synapses = False

            if "network/hyperVoxelIDs" in f:
                data["hyperVoxelIDs"] = f["network/hyperVoxelIDs"][()]

            if load_synapses:
                # 0: sourceCellID, 1: destCellID, 2: voxelX, 3: voxelY, 4: voxelZ,
                # 5: hyperVoxelID, 6: channelModelID,
                # 7: sourceAxonSomaDist (not SI scaled 1e6, micrometers),
                # 8: destDendSomaDist (not SI scalled 1e6, micrometers)
                # 9: destSegID, 10: destSegX (int 0 - 1000, SONATA wants float 0.0-1.0)
                # 11: conductance (int, not SI scaled 1e12, in pS)
                # 12: parameterID
                #
                # Note on parameterID:
                # If there are n parameter sets for the particular synapse type, then
                # the ID to use is parameterID % n, this way we can reuse connectivity
                # if we add more synapse parameter sets later.

                data["synapses"] = f["network/synapses"][:]
                data["gapJunctions"] = f["network/gapJunctions"][:]

                # !!! Convert from voxel idx to coordinates
                data["synapseCoords"] = f["network/synapses"][:, 2:5] \
                                        * f["meta/voxelSize"][()] \
                                        + f["meta/simulationOrigo"][()]
            else:
                # Point the data structure to the synapses and gap junctions on file
                # This will be slower, and only work while the file is open
                data["synapses"] = f["network/synapses"]
                data["gapJunctions"] = f["network/gapJunctions"]

                # We need to keep f alive, since we did not load synapses into
                # the memory
                self.hdf5_file = f

        else:
            data["nNeurons"] = f["network/neurons/neuronID"].shape[0]
            assert data["nNeurons"] == f["network/neurons/neuronID"][-1] + 1, \
                "Internal error, something fishy with number of neurons found"

        data["configFile"] = SnuddaLoad.to_str(f["meta/configFile"][()])

        if "meta/positionFile" in f:
            data["positionFile"] = SnuddaLoad.to_str(f["meta/positionFile"][()])

        if "meta/SlurmID" in f:
            if type(f["meta/SlurmID"][()]) in [bytes, np.bytes_]:
                data["SlurmID"] = int(f["meta/SlurmID"][()].decode())
            else:
                data["SlurmID"] = int(f["meta/SlurmID"][()])

        else:
            if self.verbose:
                print("No SlurmID set, using -1")
            data["SlurmID"] = -1

        if "meta/simulationOrigo" in f:
            data["simulationOrigo"] = f["meta/simulationOrigo"][()]

        if "meta/voxelSize" in f:
            data["voxelSize"] = f["meta/voxelSize"][()]

        if "meta/axonStumpIDFlag" in f:
            data["axonStumpIDFlag"] = f["meta/axonStumpIDFlag"][()]

        data["neurons"] = self.extract_neurons(f)

        # This is for old format, update for new format
        if "parameters" in f:
            # print("Parameters found, loading")
            data["synapseRange"] = f["parameters/synapseRange"][()]
            data["gapJunctionRange"] = f["parameters/gapJunctionRange"][()]
            data["minSynapseSpacing"] = f["parameters/minSynapseSpacing"][()]

        data["neuronPositions"] = f["network/neurons/position"][()]
        data["name"] = [SnuddaLoad.to_str(x) for x in f["network/neurons/name"][()]]

        if "populationUnitID" in f["network/neurons"]:
            data["populationUnit"] = f["network/neurons/populationUnitID"][()]
        else:
            if self.verbose:
                print("No Population Units detected.")
            data["populationUnit"] = np.zeros(data["nNeurons"], dtype=int)

        if load_morph and "morphologies" in f:
            data["morph"] = dict([])

            for morph_name in f["morphologies"].keys():
                data["morph"][morph_name] = {"swc": f["morphologies"][morph_name]["swc"][()],
                                             "location": f["morphologies"][morph_name]["location"][()]}

        data["connectivityDistributions"] = dict([])

        if "connectivityDistributions" in f["meta"]:
            orig_connectivity_distributions = \
                json.loads(SnuddaLoad.to_str(f["meta/connectivityDistributions"][()]))

            for keys in orig_connectivity_distributions:
                (pre_type, post_type) = keys.split("$$")
                data["connectivityDistributions"][pre_type, post_type] \
                    = orig_connectivity_distributions[keys]

        if "synapses" in data:
            if "gapJunctions" in data:
                print(f"Loading {len(data['neurons'])} neurons with {data['nSynapses']} synapses"
                      f" and {data['nGapJunctions']} gap junctions")
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
    def extract_synapse_coords(gathered_synapses):

        syn_coords = dict([])

        for row in gathered_synapses:
            syn_coords[row[8]] = row[0:8]

        return syn_coords

    ############################################################################

    def extract_neurons(self, hdf5_file):

        neurons = []

        for name, neuron_id, hoc, pos, rot, dend_radius, axon_radius, virtual, vID, \
            axon_density_type, axon_density, axon_density_radius, \
            axon_density_bounds_xyz, \
            morph, parameter_id, modulation_id \
                in zip(hdf5_file["network/neurons/name"][:],
                       hdf5_file["network/neurons/neuronID"][:],
                       hdf5_file["network/neurons/hoc"][:],
                       hdf5_file["network/neurons/position"][()],
                       hdf5_file["network/neurons/rotation"][()],
                       hdf5_file["network/neurons/maxDendRadius"][:],
                       hdf5_file["network/neurons/maxAxonRadius"][:],
                       hdf5_file["network/neurons/virtualNeuron"][:],
                       hdf5_file["network/neurons/volumeID"][:],
                       hdf5_file["network/neurons/axonDensityType"][:],
                       hdf5_file["network/neurons/axonDensity"][:],
                       hdf5_file["network/neurons/axonDensityRadius"][:],
                       hdf5_file["network/neurons/axonDensityBoundsXYZ"][:],
                       hdf5_file["network/neurons/morphology"][:],
                       hdf5_file["network/neurons/parameterID"][:],
                       hdf5_file["network/neurons/modulationID"][:]):

            n = dict([])

            n["name"] = SnuddaLoad.to_str(name)

            if morph is not None:
                n["morphology"] = SnuddaLoad.to_str(morph)

            # Naming convention is TYPE_X, where XX is a number starting from 0
            n["type"] = n["name"].split("_")[0]

            n["neuronID"] = neuron_id
            n["volumeID"] = SnuddaLoad.to_str(vID)
            n["hoc"] = SnuddaLoad.to_str(hoc)

            n["position"] = pos
            n["rotation"] = rot.reshape(3, 3)
            n["maxDendRadius"] = dend_radius
            n["maxAxonRadius"] = axon_radius
            n["virtualNeuron"] = virtual

            if len(axon_density_type) > 0:
                n["axonDensityType"] = SnuddaLoad.to_str(axon_density_type)
            else:
                n["axonDensityType"] = None

            if len(axon_density) > 0:
                n["axonDensity"] = SnuddaLoad.to_str(axon_density)
            else:
                n["axonDensity"] = None

            if n["axonDensityType"] == "xyz":
                n["axonDensityBoundsXYZ"] = axon_density_bounds_xyz
            else:
                n["axonDensityBoundsXYZ"] = None

            n["axonDensityRadius"] = axon_density_radius

            n["parameterID"] = parameter_id
            n["modulationID"] = modulation_id

            neurons.append(n)

        return neurons

    ############################################################################

    def load_config_file(self):

        if self.config is None:
            config_file = self.data["configFile"]
            self.config = json.load(open(config_file, 'r'))

    ############################################################################

    def load_neuron(self, neuron_id):

        neuron_info = self.data["neurons"][neuron_id]
        self.load_config_file()

        prototype_info = self.config["Neurons"][neuron_info["name"]]

        from snudda.neurons.neuron_morphology import NeuronMorphology
        neuron = NeuronMorphology(name=neuron_info["name"],
                                  position=neuron_info["position"],
                                  rotation=neuron_info["rotation"],
                                  swc_filename=prototype_info["morphology"],
                                  mech_filename=prototype_info["mechanisms"],
                                  load_morphology=True)

        return neuron

    ############################################################################

    def synapse_iterator(self, chunk_size=1000000, data_type="synapses"):

        # data_type is "synapses" or "gapJunctions"
        assert data_type in ["synapses", "gapJunctions"]

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

        return self.synapse_iterator(chunk_size=chunk_size, data_type="gapJunctions")

    ############################################################################

    @staticmethod
    def _row_eval_post_pre(row, num_neurons):
        return row[1] * num_neurons + row[0]

    @staticmethod
    def _row_eval_post(row, num_neurons):
        return row[1]

    ############################################################################

    def find_synapses_slow(self, pre_id, n_max=1000000):

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

        synapse_coords = synapses[:, 2:5][:syn_ctr, :] * self.data["voxelSize"] + self.data["simulationOrigo"]

        return synapses[:syn_ctr, :], synapse_coords

    ############################################################################

    # Either give preID and postID, or just postID

    def find_synapses(self, pre_id=None, post_id=None, silent=True):

        if post_id is None:
            return self.find_synapses_slow(pre_id=pre_id)

        assert post_id is not None, "Must specify at least postID"

        num_rows = self.data["synapses"].shape[0]
        num_neurons = len(self.data["neuronID"])

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
        synapse_coords = synapses[:, 2:5] * self.data["voxelSize"] + self.data["simulationOrigo"]

        return synapses, synapse_coords

    ############################################################################

    # Returns cellID of all neurons of neuronType
    # OBS, random_permute is not using a controled rng, so not affected by random seed set

    def get_cell_id_of_type(self, neuron_type, num_neurons=None, random_permute=False):

        cell_id = [x["neuronID"] for x in self.data["neurons"] if x["type"] == neuron_type]

        assert not random_permute or num_neurons is not None, "random_permute is only valid when num_neurons is given"

        if num_neurons is not None:
            if random_permute:
                # Do not use this if you have a simulation with multiple
                # workers... they might randomize differently, and you might
                # fewer neurons in total than you wanted
                keep_idx = np.random.permutation(len(cell_id))

                if len(keep_idx) > num_neurons:
                    keep_idx = keep_idx[:num_neurons]

                cell_id = np.array([cell_id[x] for x in keep_idx])
            else:
                cell_id = np.array([cell_id[x] for x in range(num_neurons)])

            if len(cell_id) < num_neurons:
                if self.verbose:
                    print(f"get_cell_id_of_type: wanted {num_neurons} only got {len(cell_id)} "
                          f"neurons of type {neuron_type}")

        # Double check that all of the same type
        assert np.array([self.data["neurons"][x]["type"] == neuron_type for x in cell_id]).all()

        return cell_id

    def get_population_unit_members(self, population_unit, num_neurons=None, random_permute=False):

        cell_id = np.where(self.data["populationUnit"] == population_unit)[0]

        if num_neurons:
            if random_permute:
                cell_id = np.random.permutation(cell_id)

            if len(cell_id) > num_neurons:
                cell_id = cell_id[:num_neurons]

        # Just double check
        assert (self.data["populationUnit"][cell_id] == population_unit).all()

        return cell_id

    ############################################################################


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Load snudda network file (hdf5)")
    parser.add_argument("networkFile", help="Network file (hdf5)", type=str)
    parser.add_argument("--listN", help="Lists neurons in network",
                        action="store_true")
    parser.add_argument("--listT", type=str,
                        help="List neurons of type, --listT ? list the types.",
                        default=None)
    parser.add_argument("--listPre", help="List pre synaptic neurons",
                        type=int)
    parser.add_argument("--listPost", help="List post synaptic neurons (slow)",
                        type=int)
    parser.add_argument("--keepOpen", help="This prevents loading of synapses to memory, and keeps HDF5 file open",
                        action="store_true")

    args = parser.parse_args()

    if args.keepOpen:
        loadSynapses = False
    else:
        loadSynapses = True

    nl = SnuddaLoad(args.networkFile, load_synapses=loadSynapses)

    if args.listN:
        print("Neurons in network: ")

        for nid, name, pos in [(x["neuronID"], x["name"], x["position"]) for x in nl.data["neurons"]]:
            print("%d : %s  (x: %f, y: %f, z: %f)" % (nid, name, pos[0], pos[1], pos[2]))

    if args.listT is not None:
        if args.listT == "?":
            print("List neuron types in network:")

            nTypes = np.unique([x["type"] for x in nl.data["neurons"]])
            for nt in nTypes:
                num = len([x["type"] for x in nl.data["neurons"] if x["type"] == nt])
                print(f"{nt} ({num} total)")

        else:
            print(f"Neurons of type {args.listT}:")
            nOfType = [(x["neuronID"], x["name"]) for x in nl.data["neurons"]
                       if x["type"] == args.listT]
            for nid, name in nOfType:
                print("%d : %s" % (nid, name))

    if args.listPre is not None:
        print(f"List neurons pre-synaptic to neuronID = {args.listPre} "
              f"({nl.data['neurons'][args.listPre]['name']})")
        synapses = nl.find_synapses(post_id=args.listPre)

        if synapses[0] is None:
            # Nothing to display
            os.sys.exit(0)

        preID = np.unique(synapses[0][:, 0])

        for nid, name in [(x["neuronID"], x["name"]) for x in nl.data["neurons"]
                          if x["neuronID"] in preID]:
            nSyn = np.sum(synapses[0][:, 0] == nid)
            print("%d : %s (%d synapses)" % (nid, name, nSyn))

    if args.listPost is not None:
        print("List neurons post-synaptic to neuronID = " + str(args.listPost)
              + " (" + str(nl.data["neurons"][args.listPost]["name"]) + ")")
        synapses = nl.find_synapses(pre_id=args.listPost)
        postID = np.unique(synapses[0][:, 1])

        for nid, name in [(x["neuronID"], x["name"]) for x in nl.data["neurons"]
                          if x["neuronID"] in postID]:
            nSyn = np.sum(synapses[0][:, 1] == nid)
            print("%d : %s (%d synapses)" % (nid, name, nSyn))

        # List neurons of network

    # syn = nl.findSynapses(22,5)

    # syn2 = nl.findSynapses(postID=5)

    # cellID = nl.getCellIDofType(neuronType="FSN")
