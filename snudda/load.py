import numpy as np
import h5py
import timeit
import json
import os
from glob import glob


class Snuddaload(object):

    ############################################################################

    def __init__(self, network_file, loadSynapses=True):

        if network_file == "last":
            network_file = self.find_latest_file()

        # This variable will only be set if the synapses are not kept in
        # memory so we can access them later, otherwise the hdf5 file is
        # automatically closed
        self.hdf5File = None

        self.config = None
        self.data = self.load_HDF5(network_file, loadSynapses)
        self.network_file = network_file

    ############################################################################

    def __del__(self):

        if self.hdf5File is not None:
            try:
                self.hdf5File.close()
            except:
                print("Unable to close HDF5, alread closed?")

    ############################################################################

    def load_HDF5(self, network_file, load_synapses=True, load_morph=True):
        print("Loading " + network_file)

        start_time = timeit.default_timer()
        data = dict([])

        f = h5py.File(network_file, 'r')

        # with h5py.File(network_file,'r') as f:
        if True:  # Need f open when loadSynapses = False, "with" doesnt work then

            if "config" in f:
                print("Loading config data from HDF5")
                data["config"] = f["config"][()]
                self.config = json.loads(f["config"][()])

            # Added so this code can also load the position file, which
            # does not have the network group yet
            if "network/synapses" in f:
                data["nNeurons"] = f["network/neurons/neuronID"].shape[0]
                data["nSynapses"] = f["network/synapses"].shape[0]
                data["nGapJunctions"] = f["network/gapJunctions"].shape[0]

                if data["nSynapses"] > 100e6:
                    print(str(data["nSynapses"]) + \
                          " synapses, which is a lot, not loading them into memory!")
                    load_synapses = False

                # Deprecated ??
                if "network/GJIDoffset" in f:
                    data["GJIDoffset"] = f["network/GJIDoffset"][()]

                if "network/hyperVoxelIDs" in f:
                    data["hyperVoxelIDs"] = f["network/hyperVoxelIDs"][()]

                if load_synapses:
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
                    self.hdf5File = f

                    # data["origSynapseCoords"] = f["network/origSynapseCoords"][:]
                    # gatheredSynapses = f["network/origGJCoords"][:]
                    # data["origGJCoords"] = self.extractSynapseCoords(gatheredSynapses)
            else:
                data["nNeurons"] = f["network/neurons/neuronID"].shape[0]
                assert data["nNeurons"] == f["network/neurons/neuronID"][-1] + 1, \
                    "Internal error, something fishy with number of neurons found"

            config_file = f["meta/configFile"][()]
            if type(config_file) == bytes:
                config_file = config_file.decode()
            data["configFile"] = config_file

            if "meta/positionFile" in f:
                position_file = f["meta/positionFile"][()]

                if type(position_file) == bytes:
                    position_file = position_file.decode()

                data["positionFile"] = position_file

            if "meta/SlurmID" in f:
                if type(f["meta/SlurmID"][()]) == bytes:
                    data["SlurmID"] = int(f["meta/SlurmID"][()].decode())
                else:
                    data["SlurmID"] = int(f["meta/SlurmID"][()])

            else:
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

            if "nPopulationUnits" in f["network/neurons"]:
                data["nPopulationUnits"] = f["network/neurons/nPopulationUnits"][()]
                data["populationUnit"] = f["network/neurons/populationUnitID"][()]
                data["populationUnitPlacementMethod"] = f["network/neurons/populationUnitPlacementMethod"][()]
            else:
                print("No Population Units detected.")
                data["nPopulationUnits"] = 0
                data["populationUnit"] = np.zeros(data["nNeurons"])
                data["populationUnitPlacementMethod"] = "none"

            if load_morph and "morphologies" in f:
                data["morph"] = dict([])

                for name in f["morphologies"].keys():
                    data["morph"][name] = {"swc":
                                               f["morphologies"][name]["swc"][()],
                                           "location":
                                               f["morphologies"][name]["location"][()]}

            data["connectivityDistributions"] = dict([])
            # data["connectivityDistributionsGJ"] = dict([])

            if "connectivityDistributions" in f["meta"]:
                orig_connectivity_distributions = \
                    json.loads(f["meta/connectivityDistributions"][()])

                for keys in orig_connectivity_distributions:
                    (preType, postType) = keys.split("$$")
                    data["connectivityDistributions"][preType, postType] \
                        = orig_connectivity_distributions[keys]

            #      if("connectivityDistributionsGJ" in f["meta"]):
            #        origConnectivityDistributionsGJ = \
            #          json.loads(f["meta/connectivityDistributionsGJ"][()])
            #
            #        for keys in origConnectivityDistributionsGJ:
            #          (preType,postType) = keys.split("$$")
            #          data["connectivityDistributionsGJ"][preType,postType] \
            #            = origConnectivityDistributionsGJ[keys]

            if "synapses" in data:
                if "gapJunctions" in data:
                    print(str(len(data["neurons"])) + " neurons with " \
                          + str(data["synapses"].shape[0]) + " synapses" \
                          + " and " + str(data["gapJunctions"].shape[0]) \
                          + " gap junctions")
                else:
                    print(str(len(data["neurons"])) + " neurons with " \
                          + str(data["synapses"].shape[0]) + " synapses")

            print("Load done. " + str(timeit.default_timer() - start_time))

        if load_synapses:
            f.close()
        else:
            self.hdf5File = f

        return data

    ############################################################################

    def extract_synapse_coords(self, gatheredSynapses):

        syn_coords = dict([])

        for row in gatheredSynapses:
            syn_coords[row[8]] = row[0:8]

        return syn_coords

    ############################################################################

    def extract_neurons(self, HDF5_file):

        if "parameterID" not in HDF5_file["network/neurons"]:
            return self.extract_neurons_OLD(HDF5_file)

        neurons = []

        for name, neuron_id, hoc, pos, rot, dend_radius, axon_radius, virtual, vID, \
            axon_density_type, axon_density, axon_density_radius, \
            axon_density_bounds_xyz, \
            morph, parameter_id, modulation_id \
                in zip(HDF5_file["network/neurons/name"][:],
                       HDF5_file["network/neurons/neuronID"][:],
                       HDF5_file["network/neurons/hoc"][:],
                       HDF5_file["network/neurons/position"][()],
                       HDF5_file["network/neurons/rotation"][()],
                       HDF5_file["network/neurons/maxDendRadius"][:],
                       HDF5_file["network/neurons/maxAxonRadius"][:],
                       HDF5_file["network/neurons/virtualNeuron"][:],
                       HDF5_file["network/neurons/volumeID"][:],
                       HDF5_file["network/neurons/axonDensityType"][:],
                       HDF5_file["network/neurons/axonDensity"][:],
                       HDF5_file["network/neurons/axonDensityRadius"][:],
                       HDF5_file["network/neurons/axonDensityBoundsXYZ"][:],
                       HDF5_file["network/neurons/morphology"][:],
                       HDF5_file["network/neurons/parameterID"][:],
                       HDF5_file["network/neurons/modulationID"][:]):

            n = dict([])

            if type(name) == np.ndarray:
                # Old version of savefiles give different output
                name = name[0]
                neuron_id = neuron_id[0]
                hoc = hoc[0]
                dend_radius = dend_radius[0]
                axon_radius = axon_radius[0]

            if type(name) in [bytes, np.bytes_]:
                n["name"] = name.decode()
            else:
                n["name"] = name

            if morph is not None:
                if type(morph) in [bytes, np.bytes_]:
                    n["morphology"] = morph.decode()
                else:
                    n["morphology"] = morph

            # Naming convention is TYPE_X, where XX is a number starting from 0
            n["type"] = n["name"].split("_")[0]

            n["neuronID"] = neuron_id

            if type(vID) in [bytes, np.bytes_]:
                n["volumeID"] = vID.decode()
            else:
                n["volumeID"] = vID

            if type(hoc) in [bytes, np.bytes_]:
                n["hoc"] = hoc.decode()
            else:
                n["hoc"] = hoc

            n["position"] = pos
            n["rotation"] = rot.reshape(3, 3)
            n["maxDendRadius"] = dend_radius
            n["maxAxonRadius"] = axon_radius
            n["virtualNeuron"] = virtual

            if len(axon_density_type) == 0:
                n["axonDensityType"] = None
            elif type(axon_density_type) in [bytes, np.bytes_]:
                n["axonDensityType"] = axon_density_type.decode()
            else:
                n["axonDensityType"] = axon_density_type

            if len(axon_density) > 0:
                if type(axon_density) in [bytes, np.bytes_]:
                    n["axonDensity"] = axon_density.decode()
                else:
                    n["axonDensity"] = axon_density
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

    # OLD version does not include parameterID and modulationID

    def extract_neurons_OLD(self, HDF5file):

        neurons = []

        for name, neuron_id, hoc, pos, rot, dend_radius, axon_radius, virtual, vID, \
            axon_density_type, axon_density, axon_density_radius, \
            axon_density_bounds_xyz, \
            morph \
                in zip(HDF5file["network/neurons/name"][:],
                       HDF5file["network/neurons/neuronID"][:],
                       HDF5file["network/neurons/hoc"][:],
                       HDF5file["network/neurons/position"][()],
                       HDF5file["network/neurons/rotation"][()],
                       HDF5file["network/neurons/maxDendRadius"][:],
                       HDF5file["network/neurons/maxAxonRadius"][:],
                       HDF5file["network/neurons/virtualNeuron"][:],
                       HDF5file["network/neurons/volumeID"][:],
                       HDF5file["network/neurons/axonDensityType"][:],
                       HDF5file["network/neurons/axonDensity"][:],
                       HDF5file["network/neurons/axonDensityRadius"][:],
                       HDF5file["network/neurons/axonDensityBoundsXYZ"][:],
                       HDF5file["network/neurons/morphology"][:]):

            n = dict([])

            if type(name) == np.ndarray:
                # Old version of savefiles give different output
                name = name[0]
                neuron_id = neuron_id[0]
                hoc = hoc[0]
                dend_radius = dend_radius[0]
                axon_radius = axon_radius[0]

            if type(name) in [bytes, np.bytes_]:
                n["name"] = name.decode()
            else:
                n["name"] = name

            if morph is not None:
                if type(morph) in [bytes, np.bytes_]:
                    n["morphology"] = morph.decode()
                else:
                    n["morphology"] = morph

            # Naming convention is TYPE_X, where XX is a number starting from 0
            n["type"] = n["name"].split("_")[0]

            n["neuronID"] = neuron_id

            if type(vID) in [bytes, np.bytes_]:
                n["volumeID"] = vID.decode()
            else:
                n["volumeID"] = vID

            if type(hoc) in [bytes, np.bytes_]:
                n["hoc"] = hoc.decode()
            else:
                n["hoc"] = hoc

            n["position"] = pos
            n["rotation"] = rot.reshape(3, 3)
            n["maxDendRadius"] = dend_radius
            n["maxAxonRadius"] = axon_radius
            n["virtualNeuron"] = virtual

            if len(axon_density_type) == 0:
                n["axonDensityType"] = None
            elif type(axon_density_type) in [bytes, np.bytes_]:
                n["axonDensityType"] = axon_density_type.decode()
            else:
                n["axonDensityType"] = axon_density_type

            if len(axon_density) > 0:
                if type(axon_density) in [bytes, np.bytes_]:
                    n["axonDensity"] = axon_density.decode()
                else:
                    n["axonDensity"] = axon_density
            else:
                n["axonDensity"] = None

            if n["axonDensityType"] == "xyz":
                n["axonDensityBoundsXYZ"] = axon_density_bounds_xyz
            else:
                n["axonDensityBoundsXYZ"] = None

            n["axonDensityRadius"] = axon_density_radius

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

        prototype_info = self.config[neuron_info["name"]]

        from snudda.Neuron_morphology import NeuronMorphology
        neuron = NeuronMorphology(name=neuron_info["name"],
                                  position=neuron_info["position"],
                                  rotation=neuron_info["rotation"],
                                  swc_filename=prototype_info["morphology"],
                                  mech_filename=prototype_info["mechanisms"],
                                  load_morphology=True)

        return neuron

    ############################################################################

    def synapse_iterator(self, chunk_size=1000000, dataType="synapses"):

        type_str_dict = {"synapses": "network/synapses",
                       "gapJunctions": "network/gapJunctions"}
        data_str = type_str_dict[dataType]

        with h5py.File(self.network_file, 'r') as f:
            num_rows = f[data_str].shape[0]
            if num_rows == 0:
                # No synapses
                return

            chunk_size = min(num_rows, chunk_size)
            num_steps = int(np.ceil(num_rows / chunk_size))

            row_start = 0

            for row_end in np.linspace(chunk_size, num_rows, num_steps, dtype=int):
                synapses = f[data_str][row_start:row_end, :]
                row_start = row_end

                yield synapses

    ############################################################################

    def gap_junction_iterator(self, chunk_size=1000000):

        with h5py.File(self.network_file, 'r') as f:
            num_rows = f["network/gapJunctions"].shape[0]
            chunk_size = min(num_rows, chunk_size)

            num_steps = int(np.ceil(num_rows / chunk_size))
            row_start = 0

            for rowEnd in np.linspace(chunk_size, num_rows, num_steps, dtype=int):
                gap_junctions = f["network/gapJunctions"][row_start:rowEnd, :]
                row_start = rowEnd

                yield gap_junctions

    ############################################################################

    def _row_eval_post_pre(self, row, num_neurons):
        return row[1] * num_neurons + row[0]

    def _row_eval_post(self, row, num_neurons):
        return row[1]

    ############################################################################

    def find_synapses_SLOW(self, pre_id, n_max=1000000):

        print("Finding synapses originating from " + str(pre_id) + ", this is slow")

        synapses = np.zeros((n_max, 13), dtype=np.int32)
        syn_ctr = 0

        if np.issubdtype(type(pre_id), np.integer):
            for synList in self.synapse_iterator():
                for syn in synList:
                    if syn[0] == pre_id:
                        synapses[syn_ctr, :] = syn
                        syn_ctr += 1

        else:
            for synList in self.synapse_iterator():
                for syn in synList:
                    if syn[0] in pre_id:
                        synapses[syn_ctr, :] = syn
                        syn_ctr += 1

        with h5py.File(self.network_file, 'r') as f:
            synapseCoords = synapses[:, 2:5][:syn_ctr, :] \
                            * f["meta/voxelSize"][()] \
                            + f["meta/simulationOrigo"][()]

        return synapses[:syn_ctr, :], synapseCoords

    ############################################################################

    # Either give preID and postID, or just postID

    def findSynapses(self, preID=None, postID=None, silent=True):

        if postID is None:
            return self.find_synapses_SLOW(pre_id=preID)

        with h5py.File(self.network_file, 'r') as f:

            assert postID is not None, "Must specify at least postID"

            nRows = f["network/synapses"].shape[0]
            nNeurons = f["network/neurons/neuronID"].shape[0]

            if preID is None:
                rowEval = self._row_eval_post
                valTarget = postID
            else:
                rowEval = self._row_eval_post_pre
                valTarget = postID * nNeurons + preID

            idxA1 = 0
            idxA2 = nRows - 1

            idxFound = None

            # We use idxA1 and idxA2 as upper and lower range within which we
            # hope to find one of the synapses. Once we found a synapse row
            # we go up and down in matrix to find the range of the synapses
            # matching the requested condition. This works because the matrix is
            # sorted on postID, and then preID if postID matches

            if rowEval(f["network/synapses"][idxA1, :], nNeurons) == valTarget:
                idxFound = idxA1

            if rowEval(f["network/synapses"][idxA2, :], nNeurons) == valTarget:
                idxFound = idxA2

            while idxA1 < idxA2 and idxFound is None:

                idxNext = int(np.round((idxA1 + idxA2) / 2))
                valNext = rowEval(f["network/synapses"][idxNext, :], nNeurons)

                # print("synRow = " + str( f["network/synapses"][idxNext,:]))
                # print("valTarget = " + str(valTarget) + " , valNext = " + str(valNext))
                # print("idxNext = " + str(idxNext), " - " + str(idxA1) + " " + str(idxA2))

                if valNext < valTarget:
                    idxA1 = idxNext
                elif valNext > valTarget:
                    idxA2 = idxNext
                else:
                    # We found a hit
                    idxFound = idxNext
                    break

            if idxFound is None:
                # No synapses found
                print("No synapses found")
                return None, None

            # Find start of synapse range
            idxB1 = idxFound
            valB1 = rowEval(f["network/synapses"][idxB1 - 1, :], nNeurons)

            while valB1 == valTarget and idxB1 > 0:
                idxB1 -= 1
                valB1 = rowEval(f["network/synapses"][idxB1 - 1, :], nNeurons)

            # Find end of synapse range
            idxB2 = idxFound

            if idxB2 + 1 < f["network/synapses"].shape[0]:
                valB2 = rowEval(f["network/synapses"][idxB2 + 1, :], nNeurons)

                while valB2 == valTarget and idxB2 + 1 < f["network/synapses"].shape[0]:
                    idxB2 += 1
                    valB2 = rowEval(f["network/synapses"][idxB2 + 1, :], nNeurons)

            synapses = f["network/synapses"][idxB1:idxB2 + 1, :].copy()

            if not silent:
                print("Synapse range, first " + str(idxB1) + ", last " + str(idxB2))
                print(str(synapses))

            # Calculate coordinates
            synapse_coords = synapses[:, 2:5] \
                            * f["meta/voxelSize"][()] \
                            + f["meta/simulationOrigo"][()]

            return synapses, synapse_coords

    ############################################################################

    def find_latest_file(self):

        files = glob('save/network-connect-voxel-pruned-synapse-file-*.hdf5')

        mod_time = [os.path.getmtime(f) for f in files]
        idx = np.argsort(mod_time)

        print("Using the newest file: " + files[idx[-1]])

        return files[idx[-1]]

    ############################################################################

    # Returns cellID of all neurons of neuronType

    def get_cell_id_of_type(self, neuron_type, num_neurons=None, random_permute=False):

        cell_id = [x["neuronID"] for x in self.data["neurons"] \
                   if x["type"] == neuron_type]

        if num_neurons is not None:
            if random_permute:
                # Do not use this if you have a simulation with multiple
                # workers... they might randomize differently, and you might
                # get more or less neurons in total than you wanted
                keep_idx = np.random.permutation(len(cell_id))[:num_neurons]
                cell_id = np.array([cell_id[x] for x in keep_idx])
            else:
                cell_id = np.array([cell_id[x] for x in range(num_neurons)])

            if len(cell_id) < num_neurons:
                print("getCellIDofType: wanted " + str(num_neurons) \
                      + " only got " + str(len(cell_id)) \
                      + " neurons of type " + str(neuron_type))

        # Double check that all of the same type
        assert np.array([self.data["neurons"][x]["type"] == neuron_type \
                         for x in cell_id]).all()

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

    nl = Snuddaload(args.networkFile, loadSynapses=loadSynapses)

    if args.listN:
        print("Neurons in network: ")

        for nid, name, pos in [(x["neuronID"], x["name"], x["position"]) \
                               for x in nl.data["neurons"]]:
            print("%d : %s  (x: %f, y: %f, z: %f)" % (nid, name, pos[0], pos[1], pos[2]))

    if args.listT is not None:
        if args.listT == "?":
            print("List neuron types in network:")

            nTypes = np.unique([x["type"] for x in nl.data["neurons"]])
            for nt in nTypes:
                num = len([x["type"] for x in nl.data["neurons"] if x["type"] == nt])
                print(nt + " (" + str(num) + " total)")

        else:
            print("Neurons of type " + args.listT + ":")
            nOfType = [(x["neuronID"], x["name"]) for x in nl.data["neurons"] \
                       if x["type"] == args.listT]
            for nid, name in nOfType:
                print("%d : %s" % (nid, name))

    if args.listPre:
        print("List neurons pre-synaptic to neuronID = " + str(args.listPre) \
              + " (" + str(nl.data["neurons"][args.listPre]["name"]) + ")")
        synapses = nl.findSynapses(postID=args.listPre)
        preID = np.unique(synapses[0][:, 0])

        for nid, name in [(x["neuronID"], x["name"]) for x in nl.data["neurons"] \
                          if x["neuronID"] in preID]:
            nSyn = np.sum(synapses[0][:, 0] == nid)
            print("%d : %s (%d synapses)" % (nid, name, nSyn))

    if args.listPost:
        print("List neurons post-synaptic to neuronID = " + str(args.listPost) \
              + " (" + str(nl.data["neurons"][args.listPost]["name"]) + ")")
        synapses = nl.findSynapses(preID=args.listPost)
        postID = np.unique(synapses[0][:, 1])

        for nid, name in [(x["neuronID"], x["name"]) for x in nl.data["neurons"] \
                          if x["neuronID"] in postID]:
            nSyn = np.sum(synapses[0][:, 1] == nid)
            print("%d : %s (%d synapses)" % (nid, name, nSyn))

        # List neurons of network

    # syn = nl.findSynapses(22,5)

    # syn2 = nl.findSynapses(postID=5)

    # cellID = nl.getCellIDofType(neuronType="FSN")
