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

import numpy as np
import os
from collections import OrderedDict
import h5py
import json

from .neuron_morphology import NeuronMorphology
from .region_mesh import RegionMesh

''' This code places all neurons in space, but does not setup their
  connectivity. That is done in another script. '''


class SnuddaPlace(object):

    def __init__(self,
                 config_file=None,
                 network_path=None,
                 verbose=True,
                 log_file=None,
                 d_view=None,
                 lb_view=None,
                 h5libver="latest",
                 raytrace_borders=False,
                 random_seed=None):

        if not config_file and network_path:
            config_file = os.path.join(network_path, "network-config.json")

        if not network_path and config_file:
            network_path = os.path.dirname(config_file)

        self.network_path = network_path
        self.config_file = config_file

        self.verbose = verbose
        self.log_file = log_file

        self.d_view = d_view
        self.lb_view = lb_view

        self.h5libver = h5libver
        self.write_log("Using hdf5 version: " + str(self.h5libver))

        # List of all neurons
        self.neurons = []
        self.neuronPrototypes = {}
        self.random_seed = random_seed
        self.random_generator = None

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

    ############################################################################

    def write_log(self, text):
        if self.log_file is not None:
            self.log_file.write(text + "\n")
            print(text)
        else:
            if self.verbose:
                print(text)

    ############################################################################

    def add_neurons(self,
                    swc_filename,
                    num_neurons,
                    param_data=None,
                    mech_filename=None,
                    modulation=None,
                    name="Unnamed",
                    hoc=None,
                    volume_id=None,
                    rotation_mode="random",
                    virtual_neuron=False,
                    axon_density=None):

        assert volume_id is not None, "You must specify a volume for neuron " + name

        nm = NeuronMorphology(swc_filename=swc_filename,
                              param_data=param_data,
                              mech_filename=mech_filename,
                              name=name,
                              hoc=hoc,
                              virtual_neuron=virtual_neuron)

        neuron_type = name.split("_")[0]
        neuron_coords = self.volume[volume_id]["mesh"].place_neurons(num_neurons, neuron_type)

        first_added = True

        for coords in neuron_coords:
            # We set loadMorphology = False, to preserve memory
            # Only morphology loaded for nm then, to get axon and dend
            # radius needed for connectivity

            # Pick a random parameterset
            # parameter.json can be a list of lists, this allows you to select the
            # parameterset randomly
            # modulation.json is similarly formatted, pick a parameter set here
            parameter_id = self.random_generator.integers(1000000)
            modulation_id = self.random_generator.integers(1000000)

            if rotation_mode == "random":
                rotation = nm.rand_rotation_matrix(rand_nums=self.random_generator.random(size=(3,)))
            elif rotation_mode is None or rotation_mode == "":
                self.write_log("Rotation mode: None (disabled) for " + name)
                rotation = np.eye(3)
            else:
                self.write_log("Unknown rotation mode: " + str(rotation_mode)
                               + ", valid modes '' or 'random'.")
                assert False, "Unknown rotation mode: " + str(rotation_mode)

            n = nm.clone(position=coords,
                         rotation=rotation,
                         load_morphology=False,
                         parameter_id=parameter_id,
                         modulation_id=modulation_id)

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
                self.neuronPrototypes[n.name] = n

    ############################################################################

    def read_config(self, config_file=None):

        if config_file is None:
            config_file = self.config_file

        if config_file is None:
            self.write_log("No configuration file specified")
            os.sys.exit(-1)

        if not os.path.exists(config_file):
            self.write_log("Config file does not exist: " + str(config_file))
            self.write_log("Run snudda init <your directory> first")
            os.sys.exit(-1)

        self.write_log("Parsing configuration file " + config_file)

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
                self.write_log(f"Reading random see from config file: {self.random_seed}")
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

                if "meshBinWidth" not in vol_def:
                    self.write_log("No meshBinWidth specified, using 1e-4")
                    mesh_bin_width = 1e-4
                else:
                    mesh_bin_width = vol_def["meshBinWidth"]

                self.write_log(f"Using mesh_bin_width {mesh_bin_width}")

                if "-cube-mesh-" in vol_def["meshFile"]:
                    self.write_log("Cube mesh, switching to serial processing.")
                    d_view = None
                    lb_view = None
                else:
                    d_view = self.d_view
                    lb_view = self.lb_view

                if os.path.exists(vol_def["meshFile"]):
                    mesh_file = vol_def["meshFile"]
                elif os.path.exists(os.path.join(self.network_path, vol_def["meshFile"])):
                    mesh_file = os.path.join(self.network_path, vol_def["meshFile"])
                else:
                    self.write_log(f"Unable to find mesh file {vol_def['meshFile']}")
                    os.sys.exit(-1)

                self.volume[volume_id]["mesh"] \
                    = RegionMesh(mesh_file,
                                 d_view=d_view,
                                 lb_view=lb_view,
                                 raytrace_borders=self.raytrace_borders,
                                 d_min=vol_def["dMin"],
                                 bin_width=mesh_bin_width,
                                 log_file=mesh_logfile,
                                 random_seed=vol_seed[volume_id])

            self.write_log("Using dimensions from config file")

        if "PopulationUnits" in config:
            self.population_unit_placement_method = config["PopulationUnits"]["method"]
            self.nPopulationUnits = config["PopulationUnits"]["nPopulationUnits"]

            if self.population_unit_placement_method == "populationUnitSpheres":
                self.population_unit_radius = config["PopulationUnits"]["radius"]
                self.population_unit_centres = config["PopulationUnits"]["centres"]

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

            rotation_mode = definition["rotationMode"]

            if "axonDensity" in definition:
                axon_density = definition["axonDensity"]
            else:
                axon_density = None

            self.write_log(f"Adding: {num} {neuron_name}")
            self.add_neurons(name=neuron_name,
                             swc_filename=morph,
                             param_data=param,
                             mech_filename=mech,
                             modulation=modulation,
                             num_neurons=num,
                             hoc=hoc,
                             volume_id=volume_id,
                             virtual_neuron=virtual_neuron,
                             rotation_mode=rotation_mode,
                             axon_density=axon_density)

        self.config_file = config_file

        # We reorder neurons, sorting their IDs after position
        self.sort_neurons()

        if self.population_unit_placement_method is not None:
            self.define_population_units(method=self.population_unit_placement_method)

        mesh_logfile.close()

    ############################################################################

    def all_neuron_positions(self):
        n_neurons = len(self.neurons)
        pos = np.zeros((n_neurons, 3))

        for i in range(0, n_neurons):
            pos[i, :] = self.neurons[i].position

        return pos

    ############################################################################

    def all_neuron_rotations(self):

        n_neurons = len(self.neurons)
        rot = np.zeros((n_neurons, 3, 3))

        for i in range(0, n_neurons):
            rot[i, :, :] = self.neurons[i].rotation

        return rot

    ############################################################################

    def all_neuron_names(self):
        return map(lambda x: x.name, self.neurons)

    ############################################################################

    def write_data(self, file_name):

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
        neuron_group.create_dataset("hoc", (len(hoc_list),), 'S' + str(max_hoc_len), hoc_list,
                                    compression="gzip")

        swc_list = [n.swc_filename.encode("ascii", "ignore") for n in self.neurons]
        max_swc_len = max([len(x) for x in swc_list])
        neuron_group.create_dataset("morphology", (len(swc_list),), 'S' + str(max_swc_len), swc_list,
                                    compression="gzip")

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
        neuron_modulation_id = neuron_group.create_dataset("modulationID",
                                                           (len(self.neurons),),
                                                           "int",
                                                           compression="gzip")

        for (i, n) in enumerate(self.neurons):
            neuron_position[i] = n.position
            neuron_rotation[i] = n.rotation.reshape(1, 9)
            neuron_dend_radius[i] = n.max_dend_radius
            neuron_axon_radius[i] = n.max_axon_radius
            neuron_param_id[i] = n.parameter_id
            neuron_modulation_id[i] = n.modulation_id

        # Store input information
        neuron_group.create_dataset("populationUnitID", data=self.population_unit, dtype=int)
        neuron_group.create_dataset("nPopulationUnits", data=self.nPopulationUnits, dtype=int)

        if self.population_unit_placement_method is not None:
            neuron_group.create_dataset("populationUnitPlacementMethod", data=self.population_unit_placement_method)
        else:
            neuron_group.create_dataset("populationUnitPlacementMethod", data="")

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
                    os.sys.exit(-1)

        neuron_group.create_dataset("axonDensityBoundsXYZ", data=axon_density_bounds_xyz)

        pos_file.close()

    ############################################################################

    def define_population_units(self, method="random", num_population_units=None):

        if num_population_units is None:
            num_population_units = self.nPopulationUnits

        if method == "random":
            self.random_labeling()
        elif method == "populationUnitSpheres":
            self.population_unit_spheres_labeling(self.population_unit_centres, self.population_unit_radius)
        else:
            self.population_unit = np.zeros((len(self.neurons),), dtype=int)
            self.population_units = dict([])

    ############################################################################

    def random_labeling(self, num_population_units=None):

        if num_population_units is None:
            num_population_units = self.nPopulationUnits

        self.population_unit = self.random_generator.integers(num_population_units, size=len(self.neurons))

        self.population_units = dict([])

        for i in range(0, num_population_units):
            self.population_units[i] = np.where(self.population_unit == i)[0]

    ############################################################################

    def population_unit_spheres_labeling(self, population_unit_centres, population_unit_radius):

        xyz = self.all_neuron_positions()

        centres = np.array(population_unit_centres)
        self.population_unit = np.zeros((xyz.shape[0],), dtype=int)

        if centres.shape[1] == 0:
            print("No population centres specified.")
            return

        for (ctr, pos) in enumerate(xyz):
            d = [np.linalg.norm(pos - c) for c in centres]
            idx = np.argsort(d)

            if d[idx[0]] <= population_unit_radius:
                self.population_unit[ctr] = idx[0] + 1  # We reserve 0 for no channel

        num_population_units = np.max(self.population_unit) + 1

        for i in range(0, num_population_units):
            # Channel 0 is unassigned, no channel, poor homeless neurons!
            self.population_units[i] = np.where(self.population_unit == i)[0]

    ############################################################################

    def sort_neurons(self):

        # This changes the neuron IDs so the neurons are sorted along x,y or z
        xyz = self.all_neuron_positions()

        sort_idx = np.lexsort(xyz[:, [2, 1, 0]].transpose())  # x, y, z sort order

        self.write_log("Re-sorting the neuron IDs after location")

        for newIdx, oldIdx in enumerate(sort_idx):
            self.neurons[oldIdx].neuron_id = newIdx

        self.neurons = [self.neurons[x] for x in sort_idx]

        for idx, n in enumerate(self.neurons):
            assert idx == self.neurons[idx].neuron_id, \
                "Something went wrong with sorting"

    ############################################################################


if __name__ == "__main__":

    assert False, "Please use snudda.py place networks/yournetwork"
