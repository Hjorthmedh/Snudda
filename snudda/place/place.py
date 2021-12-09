# snudda_place.py
#
# Johannes Hjorth, Royal Institute of Technology (KTH)
# Human Brain Project 2019
#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907, No 945539
# (Human Brain Project SGA1, SGA2, SGA3).

#
import sys

import numexpr
import numpy as np
import scipy.cluster
import os
from collections import OrderedDict
import h5py
import json

from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.snudda_path import snudda_parse_path, snudda_path_exists
from snudda.place.region_mesh import RegionMesh
from snudda.place.rotation import SnuddaRotate

''' This code places all neurons in space, but does not setup their
    connectivity. That is done by detect.py and prune.py '''


class SnuddaPlace(object):

    """ Places neurons in 3D space. Use detect to add connections, and prune to remove redundant connections. """

    def __init__(self,
                 config_file=None,
                 network_path=None,
                 verbose=False,
                 log_file=None,
                 rc=None,
                 d_view=None,
                 h5libver=None,
                 raytrace_borders=False,
                 random_seed=None,
                 griddata_interpolation=False):  # Setting this to true is 5x slower

        """
        Constructor.

        Args:
            config_file (str) : Path to config file, e.g. network-config.json in network_path
            network_path (str) : Path to network directory
            verbose (bool) : Print extra information on screen
            log_file (str) : Log file for place
            rc : ipyparallel remote client
            d_view : ipyparallel direct view object
            h5libver : Version of h5py library
            raytrace_borders (bool) : Should positions in border voxels be raytraces, slower but more accurate
            random_seed (int) : Numpy random seed
            griddata_interpolation (bool) : Should we interpolate density data (5x slower)

        """

        if not config_file and network_path:
            config_file = os.path.join(network_path, "network-config.json")

        if not network_path and config_file:
            network_path = os.path.dirname(config_file)

        if not log_file and network_path:
            log_dir = os.path.join(network_path, "log")
            os.makedirs(log_dir, exist_ok=True)
            log_file = open(os.path.join(log_dir, "place-neurons.txt"), "w")

        self.network_path = network_path
        self.config_file = config_file

        self.verbose = verbose
        self.log_file = log_file

        if self.network_path:
            self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        else:
            self.write_log("No network_path given, not setting position_file. Remember to pass it to write_data.")
            self.position_file = None

        self.rc = rc
        self.d_view = d_view

        if self.rc and not self.d_view:
            self.d_view = self.rc.direct_view(targets='all')

        if h5libver is None:
            self.h5libver = "latest"
        else:
            self.h5libver = h5libver

        self.write_log("Using hdf5 version: " + str(self.h5libver))

        self.griddata_interpolation = griddata_interpolation

        # List of all neurons
        self.neurons = []
        self.neuron_prototypes = {}
        self.random_seed = random_seed
        self.random_generator = None
        self.rotate_helper = None

        self.raytrace_borders = raytrace_borders

        # This defines the neuron units/channels. The dictionary lists all the
        # members of each unit, the neuronChannel gives the individual neurons
        # channel membership
        self.nPopulationUnits = 1
        self.population_unit_placement_method = "random"
        self.population_units = dict([])
        self.population_unit = None

        # These are the dimensions of our space, dMin also needs a "padding"
        # region outside the space where fake neurons are placed. This is to
        # avoid boundary effects, without padding we would get too high density
        # at the edges
        self.volume = dict([])

        # self.read_config()  # -- Now called from core.py

    def __del__(self):

        if self.rc:
            # Cleanup memory on workers
            from snudda.utils import cleanup
            cleanup(self.rc, "place")

    ############################################################################

    def place(self):

        """ Place neurons in 3D space. """

        self.parse_config()
        self.write_data()

    ############################################################################

    def write_log(self, text):

        """ Write text to log file. """

        if self.log_file is not None:
            self.log_file.write(text + "\n")

        if self.verbose:
            print(text, flush=True)

    ############################################################################

    def add_neurons(self,
                    swc_path,
                    num_neurons,
                    param_filename=None,
                    mech_filename=None,
                    modulation=None,
                    name="Unnamed",
                    hoc=None,
                    volume_id=None,
                    virtual_neuron=False,
                    axon_density=None,
                    parameter_key=None,
                    morphology_key=None,
                    modulation_key=None):

        """
        Add neurons to volume specified.

        Args:
            swc_path (str): Path to morphology directory (or single morphology)
            num_neurons (int): Number of neurons to add
            param_filename (str): Path to parameter file
            mech_filename (str): Path to mechanism file
            modulation (str): Path to neuromodulation file
            name (str): Name of neuron population, e.g. DSPN (which will become DSPN_0, DSPN_1, etc...)
            hoc (str): Path to hoc file (currently disabled)
            volume_id (str): ID of the volume to place neurons in
            virtual_neuron (bool): Real or virtual neuron, the latter can be used to model axons giving input to network
            axon_density (str): Axon density

            parameter_key (str, optional): Parameter Key to use, default=None (Randomise between all available)
            morphology_key (str, optional): Morphology Key to use, default=None (Randomise between all available)
            modulation_key (str, optional): Key for neuromodulation, default=None (Randomise between all available)
        """

        assert volume_id is not None, f"You must specify a volume for neuron {name}"
        assert hoc is None, "Currently only support hoc=None, since we can have multiple parameter, morph combos now"

        neuron_prototype = NeuronPrototype(neuron_name=name,
                                           neuron_path=None,
                                           morphology_path=swc_path,
                                           parameter_path=param_filename,
                                           mechanism_path=mech_filename,
                                           modulation_path=modulation,
                                           load_morphology=False,
                                           virtual_neuron=virtual_neuron)

        neuron_type = name.split("_")[0]
        neuron_positions = self.volume[volume_id]["mesh"].place_neurons(num_neurons, neuron_type)

        first_added = True

        neuron_rotations = self.rotate_helper.get_rotations(volume_name=volume_id, neuron_type=neuron_type,
                                                            neuron_positions=neuron_positions,
                                                            rng=self.random_generator)

        for coords, rotation in zip(neuron_positions, neuron_rotations):
            # We set loadMorphology = False, to preserve memory
            # Only morphology loaded for nm then, to get axon and dend
            # radius needed for connectivity

            # Pick a random parameter set
            # parameter.json can be a list of lists, this allows you to select the
            # parameter set randomly
            # modulation.json is similarly formatted, pick a parameter set here
            if parameter_key is None:
                parameter_id = self.random_generator.integers(1000000)
            else:
                parameter_id = None

            if modulation_key is None:
                modulation_id = self.random_generator.integers(1000000)
            else:
                modulation_id = None

            if morphology_key is None:
                morphology_id = self.random_generator.integers(1000000)
            else:
                morphology_id = None

            n = neuron_prototype.clone(position=coords,
                                       rotation=rotation,
                                       morphology_id=morphology_id,
                                       parameter_id=parameter_id,
                                       modulation_id=modulation_id,
                                       parameter_key=parameter_key,
                                       morphology_key=morphology_key,
                                       modulation_key=modulation_key)

            # self.writeLog("Place " + str(self.cellPos[i,:]))

            n.neuron_id = len(self.neurons)
            n.volume_id = volume_id

            assert axon_density is None or len(n.axon) == 0, \
                "!!! ERROR: Neuron: " + str(n.name) + " has both axon and axon density."

            n.axon_density = axon_density
            self.neurons.append(n)

            # This info is used by workers to speed things up
            if first_added:
                first_added = False
                self.neuron_prototypes[n.name] = n

    ############################################################################

    def parse_config(self, config_file=None, resort_neurons=True):

        """ Parse network config_file """

        if config_file is None:
            config_file = self.config_file

        if config_file is None:
            self.write_log("No configuration file specified")
            sys.exit(-1)

        if not os.path.exists(config_file):
            self.write_log(f"Config file does not exist: {config_file}")
            self.write_log("Run snudda init <your directory> first")
            sys.exit(-1)

        self.write_log(f"Parsing configuration file {config_file}")

        cfg_file = open(config_file, 'r')

        try:
            config = json.load(cfg_file, object_pairs_hook=OrderedDict)
        finally:
            cfg_file.close()

        if not config:
            self.write_log("Warning, empty network config.")

        if self.random_seed is None:
            if "RandomSeed" in config and "place" in config["RandomSeed"]:
                self.random_seed = config["RandomSeed"]["place"]
                self.write_log(f"Reading random seed from config file: {self.random_seed}")
            else:
                # No random seed given, invent one
                self.random_seed = 1001
                self.write_log(f"No random seed provided, using: {self.random_seed}")
        else:
            self.write_log(f"Using random seed provided by command line: {self.random_seed}")

        self.random_generator = np.random.default_rng(self.random_seed + 115)

        if self.log_file is None:
            mesh_log_filename = "mesh-log.txt"
        else:
            mesh_log_filename = self.log_file.name + "-mesh"
        mesh_logfile = open(mesh_log_filename, 'wt')

        # First handle volume definitions
        volume_def = config["Volume"]

        # Setup random seeds for all volumes
        ss = np.random.SeedSequence(self.random_seed)
        all_seeds = ss.generate_state(len(volume_def))
        all_vd = sorted(volume_def.keys())

        vol_seed = dict()
        for vd, seed in zip(all_vd, all_seeds):
            vol_seed[vd] = seed

        for volume_id, vol_def in volume_def.items():

            self.volume[volume_id] = vol_def

            if "meshFile" in vol_def:

                assert "dMin" in vol_def, "You must specify dMin if using a mesh" \
                                          + " for volume " + str(volume_id)

                if "meshBinWidth" not in vol_def or not vol_def["meshBinWidth"]:
                    self.write_log("No meshBinWidth specified, using 1e-4")
                    mesh_bin_width = 1e-4
                else:
                    mesh_bin_width = vol_def["meshBinWidth"]

                self.write_log(f"Using mesh_bin_width {mesh_bin_width}")

                if "-cube-mesh-" in vol_def["meshFile"] or "slice.obj" in vol_def["meshFile"]:
                    self.write_log("Cube or slice mesh, switching to serial processing.")
                    d_view = None
                else:
                    d_view = self.d_view

                if snudda_path_exists(vol_def["meshFile"]):
                    mesh_file = snudda_parse_path(vol_def["meshFile"])
                elif os.path.exists(os.path.join(self.network_path, vol_def["meshFile"])):
                    mesh_file = os.path.join(self.network_path, vol_def["meshFile"])
                else:
                    self.write_log(f"Unable to find mesh file {vol_def['meshFile']}")
                    sys.exit(-1)

                self.volume[volume_id]["mesh"] \
                    = RegionMesh(mesh_file,
                                 d_view=d_view,
                                 raytrace_borders=self.raytrace_borders,
                                 d_min=vol_def["dMin"],
                                 bin_width=mesh_bin_width,
                                 log_file=mesh_logfile,
                                 random_seed=vol_seed[volume_id])

                if "density" in self.volume[volume_id]:
                    # We need to set up the neuron density functions also
                    # TODO: Here density for each neuron type in the volume is defined
                    #       as a numexpr evaluated string of x,y,z stored in a dictionary
                    #       with the neuron type as key.
                    #       Add ability to also specify a density file.
                    for neuron_type in self.volume[volume_id]["density"]:

                        density_func = None

                        if "densityFunction" in self.volume[volume_id]["density"][neuron_type]:
                            density_str = self.volume[volume_id]["density"][neuron_type]["densityFunction"]
                            density_func = lambda x, y, z: numexpr.evaluate(density_str)

                        if "densityFile" in self.volume[volume_id]["density"][neuron_type]:
                            density_file = self.volume[volume_id]["density"][neuron_type]["densityFile"]

                            # We need to load the data from the file
                            from scipy.interpolate import griddata
                            with open(snudda_parse_path(density_file), "r") as f:
                                density_data = json.load(f, object_pairs_hook=OrderedDict)

                                assert volume_id in density_data and neuron_type in density_data[volume_id], \
                                    f"Volume {volume_id} does not contain data for neuron type {neuron_type}"

                                assert "Coordinates" in density_data[volume_id][neuron_type] \
                                       and "Density" in density_data[volume_id][neuron_type], \
                                    (f"Missing Coordinates and/or Density data for "
                                     f"volume {volume_id}, neuron type {neuron_type}")

                                coord = np.array(density_data[volume_id][neuron_type]["Coordinates"]) * 1e-6  # Convert to SI
                                density = np.array(density_data[volume_id][neuron_type]["Density"])

                                if self.griddata_interpolation:
                                    density_func_helper = lambda pos: griddata(points=coord, values=density,
                                                                               xi=pos, method="linear",
                                                                               fill_value=0)
                                else:
                                    density_func_helper = lambda pos: griddata(points=coord, values=density,
                                                                               xi=pos, method="nearest",
                                                                               fill_value=0)

                                density_func = lambda x, y, z: density_func_helper(np.array([x, y, z]).transpose())

                        self.volume[volume_id]["mesh"].define_density(neuron_type, density_func)

            self.write_log("Using dimensions from config file")

        # Setup for rotations
        self.rotate_helper = SnuddaRotate(self.config_file)

        assert "Neurons" in config, \
            "No neurons defined. Is this config file old format?"

        # Read in the neurons
        for name, definition in config["Neurons"].items():

            neuron_name = name
            morph = definition["morphology"]
            param = definition["parameters"]
            mech = definition["mechanisms"]

            if "modulation" in definition:
                modulation = definition["modulation"]
            else:
                # Modulation optional
                modulation = None

            num = definition["num"]
            volume_id = definition["volumeID"]

            if "neuronType" in definition:
                # type is "neuron" or "virtual" (provides input only)
                model_type = definition["neuronType"]
            else:
                model_type = "neuron"

            if 'hoc' in definition:
                hoc = definition["hoc"]
            else:
                hoc = None

            if model_type == "virtual":
                # Virtual neurons gets spikes from a file
                mech = ""
                hoc = ""
                virtual_neuron = True
            else:
                virtual_neuron = False

            if "axonDensity" in definition:
                axon_density = definition["axonDensity"]
            else:
                axon_density = None

            if "parameterKey" in definition:
                parameter_key = definition["parameterKey"]
            else:
                parameter_key = None

            if "morphologyKey" in definition:
                morphology_key = definition["morphologyKey"]
            else:
                morphology_key = None

            if "modulationKey" in definition:
                modulation_key = definition["modulationKey"]
            else:
                modulation_key = None

            self.write_log(f"Adding: {num} {neuron_name}")
            self.add_neurons(name=neuron_name,
                             swc_path=morph,
                             param_filename=param,
                             mech_filename=mech,
                             modulation=modulation,
                             num_neurons=num,
                             hoc=hoc,
                             volume_id=volume_id,
                             virtual_neuron=virtual_neuron,
                             axon_density=axon_density,
                             parameter_key=parameter_key,
                             morphology_key=morphology_key,
                             modulation_key=modulation_key)

        self.config_file = config_file

        # We reorder neurons, sorting their IDs after position
        # -- UPDATE: Now we spatial cluster neurons depending on number of workers
        if resort_neurons:
            self.sort_neurons(sort_idx=self.cluster_neurons())

        if False:  # Debug purposes, make sure neuron ranges are ok
            self.plot_ranges()

        if "PopulationUnits" in config:
            self.define_population_units(config["PopulationUnits"])

        mesh_logfile.close()

    ############################################################################

    def all_neuron_positions(self):

        """ Returns all neuron positions as a n x 3 matrix. """

        n_neurons = len(self.neurons)
        pos = np.zeros((n_neurons, 3))

        for i in range(0, n_neurons):
            pos[i, :] = self.neurons[i].position

        return pos

    ############################################################################

    def all_neuron_rotations(self):

        """ Returns all neuron rotations as a n x 3 x 3 matrix. """

        n_neurons = len(self.neurons)
        rot = np.zeros((n_neurons, 3, 3))

        for i in range(0, n_neurons):
            rot[i, :, :] = self.neurons[i].rotation

        return rot

    ############################################################################

    def all_neuron_names(self):

        """ Returns all neuron names as a list. """

        return map(lambda x: x.name, self.neurons)

    ############################################################################

    def write_data(self, file_name=None):

        """ Writes positition data to HDF5 file file_name. """

        if not file_name:
            file_name = self.position_file

        assert len(self.neurons) > 0, "No neurons to save!"

        self.write_log(f"Writing data to HDF5 file: {file_name}")

        pos_file = h5py.File(file_name, "w", libver=self.h5libver)

        with open(self.config_file, 'r') as cfg_file:
            config = json.load(cfg_file, object_pairs_hook=OrderedDict)

        # Meta data
        save_meta_data = [(self.config_file, "configFile"),
                          (json.dumps(config), "config")]

        meta_group = pos_file.create_group("meta")

        for data, dataName in save_meta_data:
            meta_group.create_dataset(dataName, data=data)

        network_group = pos_file.create_group("network")

        # Neuron information
        neuron_group = network_group.create_group("neurons")

        # If the name list is longer than 20 chars, increase S20
        name_list = [n.name.encode("ascii", "ignore") for n in self.neurons]
        str_type = 'S' + str(max(1, max([len(x) for x in name_list])))
        neuron_group.create_dataset("name", (len(name_list),), str_type, name_list,
                                    compression="gzip")

        neuron_id_list = np.arange(len(self.neurons))
        neuron_group.create_dataset("neuronID", (len(neuron_id_list),),
                                    'int', neuron_id_list)

        volume_id_list = [n.volume_id.encode("ascii", "ignore")
                          for n in self.neurons]
        str_type_vid = 'S' + str(max(1, max([len(x) for x in volume_id_list])))

        neuron_group.create_dataset("volumeID",
                                    (len(volume_id_list),), str_type_vid, volume_id_list,
                                    compression="gzip")

        hoc_list = [n.hoc.encode("ascii", "ignore") for n in self.neurons]
        max_hoc_len = max([len(x) for x in hoc_list])
        max_hoc_len = max(max_hoc_len, 10)  # In case there are none
        neuron_group.create_dataset("hoc", (len(hoc_list),), f"S{max_hoc_len}", hoc_list,
                                    compression="gzip")

        swc_list = [n.swc_filename.encode("ascii", "ignore") for n in self.neurons]
        max_swc_len = max([len(x) for x in swc_list])
        neuron_group.create_dataset("morphology", (len(swc_list),), f"S{max_swc_len}", swc_list,
                                    compression="gzip")

        neuron_path = [n.neuron_path.encode("ascii", "ignore") for n in self.neurons]
        max_np_len = max([len(x) for x in neuron_path])
        neuron_group.create_dataset("neuronPath", (len(neuron_path),), f"S{max_np_len}", neuron_path)

        virtual_neuron_list = np.array([n.virtual_neuron for n in self.neurons], dtype=bool)
        virtual_neuron = neuron_group.create_dataset("virtualNeuron",
                                                     data=virtual_neuron_list)

        # Create dataset, filled further down
        neuron_position = neuron_group.create_dataset("position",
                                                      (len(self.neurons), 3),
                                                      "float",
                                                      compression="gzip")
        neuron_rotation = neuron_group.create_dataset("rotation",
                                                      (len(self.neurons), 9),
                                                      "float",
                                                      compression="gzip")

        neuron_dend_radius = neuron_group.create_dataset("maxDendRadius",
                                                         (len(self.neurons),),
                                                         "float",
                                                         compression="gzip")

        neuron_axon_radius = neuron_group.create_dataset("maxAxonRadius",
                                                         (len(self.neurons),),
                                                         "float",
                                                         compression="gzip")

        neuron_param_id = neuron_group.create_dataset("parameterID",
                                                      (len(self.neurons),),
                                                      "int",
                                                      compression="gzip")

        neuron_morph_id = neuron_group.create_dataset("morphologyID",
                                                      (len(self.neurons),),
                                                      "int",
                                                      compression="gzip")

        neuron_modulation_id = neuron_group.create_dataset("modulationID",
                                                           (len(self.neurons),),
                                                           "int",
                                                           compression="gzip")

        pk_list = [n.parameter_key.encode("ascii", "ignore")
                   if n.parameter_key is not None else ""
                   for n in self.neurons]
        pk_str_type = 'S' + str(max(1, max([len(x) for x in pk_list])))

        mk_list = [n.morphology_key.encode("ascii", "ignore")
                   if n.morphology_key is not None else ""
                   for n in self.neurons]
        mk_str_type = 'S' + str(max(1, max([len(x) for x in mk_list])))

        mok_list = [n.modulation_key.encode("ascii", "ignore")
                    if n.modulation_key is not None else ""
                    for n in self.neurons]
        mok_str_type = 'S' + str(max(1, max([len(x) for x in mok_list])))

        neuron_param_key = neuron_group.create_dataset("parameterKey",
                                                       (len(self.neurons),),
                                                       pk_str_type,
                                                       compression="gzip")

        neuron_morph_key = neuron_group.create_dataset("morphologyKey",
                                                       (len(self.neurons),),
                                                       mk_str_type,
                                                       compression="gzip")

        neuron_modulation_key = neuron_group.create_dataset("modulationKey",
                                                            (len(self.neurons),),
                                                            mok_str_type,
                                                            compression="gzip")

        for (i, n) in enumerate(self.neurons):
            neuron_position[i] = n.position
            neuron_rotation[i] = n.rotation.reshape(1, 9)
            neuron_dend_radius[i] = n.max_dend_radius
            neuron_axon_radius[i] = n.max_axon_radius
            neuron_param_id[i] = -1 if n.parameter_id is None else n.parameter_id
            neuron_morph_id[i] = -1 if n.morphology_id is None else n.morphology_id
            neuron_modulation_id[i] = -1 if n.modulation_id is None else n.modulation_id

            if n.parameter_key:
                neuron_param_key[i] = n.parameter_key

            if n.morphology_key:
                neuron_morph_key[i] = n.morphology_key

            if n.modulation_key:
                neuron_modulation_key[i] = n.modulation_key

        # Store input information
        if self.population_unit is None:
            # If no population units were defined, then set them all to 0 (= no population unit)
            self.population_unit = np.zeros((len(self.neurons),), dtype=int)

        neuron_group.create_dataset("populationUnitID", data=self.population_unit, dtype=int)
        # neuron_group.create_dataset("nPopulationUnits", data=self.nPopulationUnits, dtype=int)

        # Variable for axon density "r", "xyz" or "" (No axon density)
        axon_density_type = [n.axon_density[0].encode("ascii", "ignore")
                             if n.axon_density is not None
                             else b""
                             for n in self.neurons]

        ad_str_type2 = "S" + str(max(1, max([len(x) if x is not None else 1
                                             for x in axon_density_type])))
        neuron_group.create_dataset("axonDensityType", (len(axon_density_type),),
                                    ad_str_type2, data=axon_density_type,
                                    compression="gzip")

        axon_density = [n.axon_density[1].encode("ascii", "ignore")
                        if n.axon_density is not None
                        else b""
                        for n in self.neurons]
        ad_str_type = "S" + str(max(1, max([len(x) if x is not None else 1
                                            for x in axon_density])))

        neuron_group.create_dataset("axonDensity", (len(axon_density),),
                                    ad_str_type, data=axon_density, compression="gzip")

        axon_density_radius = [n.axon_density[2]
                               if n.axon_density is not None and n.axon_density[0] == "r"
                               else np.nan for n in self.neurons]

        neuron_group.create_dataset("axonDensityRadius", data=axon_density_radius)

        # We also need to save axonDensityBoundsXYZ, and nAxon points for the
        # non-spherical axon density option

        axon_density_bounds_xyz = np.nan * np.zeros((len(self.neurons), 6))

        for ni, n in enumerate(self.neurons):

            if n.axon_density is None:
                # No axon density specified, skip
                continue

            if n.axon_density[0] == "xyz":

                try:
                    axon_density_bounds_xyz[ni, :] = np.array(n.axon_density[2])
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    self.write_log(f"Incorrect density string: {n.axon_density}")
                    sys.exit(-1)

        neuron_group.create_dataset("axonDensityBoundsXYZ", data=axon_density_bounds_xyz)

        pos_file.close()

    ############################################################################

    def define_population_units(self, population_unit_info):

        """ Defines population units.

        Args:
            population_unit_info (dict): Has keys "AllUnitID" with a list of all Unit IDs, and <VolumeID> which points
                                         to a dictionary. This dictionary has keys:
                                         "method" : "random" or "radialDensity"                    "
                                         "unitID" : ID of population unit
                                         "fractionOfNeurons" : How large fraction of neurons belong to this unit (used by "random" method)
                                         "neuronTypes" : List of Neuron types that belong to this population unit
                                         "numNeurons" : Number of neurons in each population unit, only used with radialDensity method
                                         "structure" : Name of structure population unit is located in (VolumeID)
                                         "centres" : Centre of radial density
                                         "ProbabilityFunctions" : Probability function defining unit membership, function of radius

        """

        method_lookup = {"random": self.random_labeling,
                         "radialDensity": self.population_unit_density_labeling}

        for unit_id in population_unit_info["AllUnitID"]:
            self.population_units[unit_id] = []

        for volume_id in population_unit_info:
            if volume_id in ["AllUnitID"]:
                continue  # Not a population unit, metadata.

            neuron_id = self.volume_neurons(volume_id)
            method_name = population_unit_info[volume_id]["method"]

            assert method_name in method_lookup, \
                (f"Unknown population placement method {method_name}. "
                 f"Valid options are {', '.join([x for x in method_lookup])}")

            method_lookup[method_name](population_unit_info[volume_id], neuron_id)

    ############################################################################

    def random_labeling(self, population_unit_info, neuron_id):

        """
        Creates random labeling.

        Args:
             neuron_id (list) : All potential neuron ID
             population_unit_info (dict): Dictionary with "method" = "random"
                                          "unitID" = ID of population unit
                                          "fractionOfNeurons" = Fraction of neurons in this unit
                                          "neuronTypes" = List of neuron types that are in this unit
                                          "structure" = Name of the structure the unit is located in

        """

        self.init_population_units()  # This initialises population unit labelling if not already allocated

        unit_id = population_unit_info["unitID"]
        fraction_of_neurons = population_unit_info["fractionOfNeurons"]
        neuron_types = population_unit_info["neuronTypes"]  # list of neuron types that belong to this population unit
        structure_name = population_unit_info["structure"]

        # First we need to generate unit lists (with fractions) for each neuron type
        units_available = dict()
        for uid, fon, neuron_type_list in zip(unit_id, fraction_of_neurons, neuron_types):
            for nt in neuron_type_list:
                if nt not in units_available:
                    units_available[nt] = dict()
                    units_available[nt]["unit"] = []
                    units_available[nt]["fraction"] = []

                units_available[nt]["unit"].append(uid)
                units_available[nt]["fraction"].append(fon)

        all_neuron_types = [n.name.split("_")[0] for n in self.neurons]
        # Next we check that no fraction sum is larger than 1

        for neuron_type in units_available:
            assert np.sum(units_available[neuron_type]["fraction"]) <= 1, \
                (f"Population unit fraction sum for Neuron type {neuron_type} "
                 f"in structure {structure_name} sums to more than 1.")

            assert (np.array(units_available[neuron_type]["fraction"]) >= 0).all(), \
                f"Population unit fractions must be >= 0. Please check {neuron_type} in {structure_name}"

            cum_fraction = np.cumsum(units_available[neuron_type]["fraction"])

            neurons_of_type = [self.neurons[nid].neuron_id
                               for (nid, n_type) in zip(neuron_id, all_neuron_types)
                               if n_type == neuron_type]

            rand_num = self.random_generator.uniform(size=len(neurons_of_type))

            for nid, rn in zip(neurons_of_type, rand_num):
                # If our randum number is smaller than the first fraction, then neuron in first pop unit
                # if random number is between first and second cumulative fraction, then second pop unit
                # If larger than last cum_fraction, then no pop unit was picked (and we get -1)
                idx = len(cum_fraction) - np.sum(rn <= cum_fraction)

                if idx < len(cum_fraction):
                    unit_id = units_available[nt]["unit"][idx]
                    self.population_unit[nid] = unit_id
                    self.population_units[unit_id].append(nid)
                else:
                    self.population_unit[nid] = 0

    ############################################################################

    def population_unit_density_labeling(self, population_unit_info, neuron_id):

        """
        Creates population units based on radial density functions.

        Args:
            neuron_id (list) : All potential neuron ID
            population_unit_info (dict): "method" must be "radialDensity"
                                         "neuronTypes" list of neuron types
                                         "centres" of radial probabilityes, one per neuron type
                                         "numNeurons" : Number of neurons in each population unit, only used with radialDensity method
                                         "probabilityFunctions" list of probability functions of r (as str)
                                         "unitID" ID of population unit
        """

        assert population_unit_info["method"] == "radialDensity"
        self.init_population_units()  # This initialises population unit labelling if not alraedy allocated

        neuron_types = population_unit_info["neuronTypes"]
        centres = np.array(population_unit_info["centres"])
        probability_functions = population_unit_info["ProbabilityFunctions"]
        unit_id = population_unit_info["unitID"]

        if "numNeurons" in population_unit_info and population_unit_info["numNeurons"] is not None:
            num_neurons = population_unit_info["numNeurons"]

            if np.isscalar(num_neurons):
                num_neurons = np.full((len(centres),), num_neurons)
        else:
            num_neurons = [None for x in centres]

        assert len(neuron_types) == len(centres) == len(probability_functions) == len(unit_id) == len(num_neurons)

        # xyz = self.all_neuron_positions()
        unit_probability = np.zeros(centres.shape[0])

        for nid in neuron_id:

            pos = self.neurons[nid].position
            neuron_type = self.neurons[nid].name.split("_")[0]

            for idx, (centre_pos, neuron_type_list, p_func) \
                    in enumerate(zip(centres, neuron_types, probability_functions)):

                if neuron_type in neuron_type_list:
                    d = np.linalg.norm(pos-centre_pos)
                    unit_probability[idx] = numexpr.evaluate(p_func)
                else:
                    unit_probability[idx] = 0  # That unit does not contain this neuron type

            # Next we randomise membership
            rand_num = self.random_generator.uniform(size=len(unit_probability))
            member_flag = rand_num < unit_probability

            # Currently we only allow a neuron to be member of one population unit
            n_flags = np.sum(member_flag)
            if n_flags == 0:
                self.population_unit[nid] = 0
            elif n_flags == 1:
                self.population_unit[nid] = unit_id[np.where(member_flag)[0][0]]
            else:
                # More than one unit, pick the one that had smallest relative randnum
                idx = np.argmax(np.multiply(np.divide(unit_probability, rand_num), member_flag))
                self.population_unit[nid] = unit_id[idx]

        # Also update dictionary with lists of neurons of that unit
        for uid in unit_id:
            # Channel 0 is unassigned, no channel, poor homeless neurons!
            self.population_units[uid] = np.where(self.population_unit == uid)[0]

        # Finally if numNeurons is specified, we need to reduce the number of neurons belonging to that unit
        for u_id, n_neurons in zip(unit_id, num_neurons):
            if n_neurons is not None:
                assert len(self.population_units[u_id]) >= n_neurons, \
                    f"Unable to pick {n_neurons} for population unit {u_id}, only {len(self.population_units[u_id])} available."

                if n_neurons < len(self.population_units[u_id]):

                    perm_nid = np.random.permutation(self.population_units[u_id])
                    keep_nid = perm_nid[:n_neurons]
                    self.population_units[u_id] = keep_nid

                    remove_nid = perm_nid[n_neurons:]
                    for rid in remove_nid:
                        self.population_unit[rid] = 0

                    if 0 in self.population_units:
                        self.population_units[0] = np.array(list(set(self.population_units[0]).union(remove_nid)))
                    else:
                        self.population_units[0] = remove_nid


    ############################################################################

    def init_population_units(self):

        """ Initialise population units. If none are given they are all set to 0."""

        if not self.population_unit:
            # If no population units were defined, then set them all to 0 (= no population unit)
            self.population_unit = np.zeros((len(self.neurons),), dtype=int)

    # TODO: In prune we later need to gather all synapses belonging to each neuron, which means opening
    #       the hypervoxel files that contain the worker's neurons' synapses. Therefore it is good to have
    #       only nearby neurons on the same worker. Come up with a better scheme for sorting neurons.
    #

    def plot_ranges(self):

        from matplotlib import pyplot as plt

        n_workers = len(self.d_view) if self.d_view is not None else 1
        range_borders = np.linspace(0, len(self.neurons), n_workers + 1).astype(int)

        colours = np.zeros((len(self.neurons),))
        r_start = 0
        for idx, r_end in enumerate(range_borders[1:]):
            colours[r_start:r_end] = idx+1
            r_start = r_end

        xyz = self.all_neuron_positions()

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=colours, alpha=0.5)
        plt.show()

    def cluster_neurons(self, n_trials=5):

        """
        Cluster neurons, so that nearby neurons are grouped on same worker, to speed up simulations.

        Args:
            n_trials (int) : Number of trials for k-means clustering (default 5)
        """

        n_workers = len(self.d_view) if self.d_view is not None else 1
        n_clusters = np.maximum(n_workers*5, 100)
        n_clusters = np.minimum(n_clusters, len(self.neurons))

        xyz = self.all_neuron_positions()
        centroids, labels = scipy.cluster.vq.kmeans2(xyz, n_clusters, minit="points")

        n_centroids = centroids.shape[0]
        assert n_centroids == n_clusters
        cluster_member_list = [[] for x in range(n_centroids)]

        # Find the members of each cluster
        for idx, val in enumerate(labels):
            cluster_member_list[val].append(idx)

        num_neurons = xyz.shape[0]
        range_borders = np.linspace(0, num_neurons, n_workers + 1).astype(int)

        global_centroid_order = list(np.argsort(np.sum(centroids, axis=1)))

        neuron_order = -np.ones((num_neurons,), dtype=int)
        neuron_order_ctr = 0
        range_start = 0

        for range_end in range_borders[1:]:

            # While within a range, we want the clusters closest to the current primary cluster
            current_cluster = global_centroid_order[0]
            d = np.linalg.norm(centroids - centroids[current_cluster, :], axis=1)
            local_centroid_order = list(np.argsort(d))

            while range_start < range_end and len(local_centroid_order) > 0:
                while len(cluster_member_list[local_centroid_order[0]]) == 0:
                    del local_centroid_order[0]
                current_cluster = local_centroid_order[0]

                take_n = np.minimum(len(cluster_member_list[current_cluster]), range_end - range_start)
                neuron_order[neuron_order_ctr:neuron_order_ctr+take_n] = cluster_member_list[current_cluster][:take_n]
                neuron_order_ctr += take_n

                del cluster_member_list[current_cluster][:take_n]
                range_start += take_n

            while len(global_centroid_order) > 0 and len(cluster_member_list[global_centroid_order[0]]) == 0:
                del global_centroid_order[0]

            if len(global_centroid_order) == 0:
                break

        # Sometimes the original cluster is bad? Try again...
        if np.count_nonzero(neuron_order < 0) > 0 and n_trials > 1:
            self.write_log(f"Redoing place:neuron_clustering, {np.count_nonzero(neuron_order < 0)} "
                           f"neurons unaccounted for",
                           is_error=True)
            self.write_log(f"incorrect neuron_order={neuron_order} (printed for debugging)")
            neuron_order = self.cluster_neurons(n_trials=n_trials-1)

        # TODO: This occured once on Tegner, why did it happen?
        assert np.count_nonzero(neuron_order < 0) == 0, \
            "cluster_neurons: Not all neurons accounted for. Please rerun place."

        # Just some check that all is ok
        assert (np.diff(np.sort(neuron_order)) == 1).all(), "cluster_neurons: There are gaps in the sorting, error"

        # TODO: Verify that sort order is ok

        return neuron_order

    def sort_neurons(self, sort_idx=None):

        """ Sorting neurons. If no argument is given they will be sorted along x,y,z axis.

            To use cluster sorting, use:
                sp.sort_neurons(sort_idx=sp.cluster_neurons())

        Args:
            sort_idx (list, optional) : Sort order

            """

        if sort_idx is None:
            # This changes the neuron IDs so the neurons are sorted along x,y or z
            xyz = self.all_neuron_positions()
            sort_idx = np.lexsort(xyz[:, [2, 1, 0]].transpose())  # x, y, z sort order

        self.write_log("Re-sorting the neuron IDs")

        for newIdx, oldIdx in enumerate(sort_idx):
            self.neurons[oldIdx].neuron_id = newIdx

        self.neurons = [self.neurons[x] for x in sort_idx]

        for idx, n in enumerate(self.neurons):
            assert idx == self.neurons[idx].neuron_id, \
                "Something went wrong with sorting"

    def volume_neurons(self, volume_id):

        return [n.neuron_id for n in self.neurons if n.volume_id == volume_id]

    ############################################################################


if __name__ == "__main__":

    assert False, "Please use snudda.py place networks/yournetwork"
