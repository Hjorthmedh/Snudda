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

# TODO: Let place also decide where separate axonal morphologies should be located, save swc morph, and position + rotation
#       Then detect can just load that info directly. Simplifies load_neuron
#       HDF5 filen, ha en matrix med info för axon morph, axon pos, axon rotation, axon parent

import json
import os
import sys

import h5py
import numexpr
import numpy as np
import scipy.cluster
from scipy.interpolate import griddata

from snudda.utils.snudda_path import get_snudda_data
from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.place.region_mesh_redux import NeuronPlacer, RegionMeshRedux

from snudda.place.rotation import SnuddaRotate
from snudda.utils.snudda_path import snudda_parse_path, snudda_path_exists, snudda_simplify_path

from snudda.neurons.morphology_data import MorphologyData

''' This code places all neurons in space, but does not setup their
    connectivity. That is done by detect.py and prune.py '''


class SnuddaPlace(object):
    """ Places neurons in 3D space. Use detect to add connections, and prune to remove redundant connections. """

    def __init__(self,
                 config_file=None,
                 network_path=None,
                 snudda_data=None,
                 verbose=False,
                 log_file=None,
                 rc=None,
                 d_view=None,
                 h5libver=None,
                 random_seed=None,
                 griddata_interpolation=False,
                 morphologies_stay_inside=True):

        """
        Constructor.

        Args:
            config_file (str) : Path to config file, e.g. network-config.json in network_path
            network_path (str) : Path to network directory
            snudda_data (str): Path to snudda data
            verbose (bool) : Print extra information on screen
            log_file (str) : Log file for place
            rc : ipyparallel remote client
            d_view : ipyparallel direct view object
            h5libver : Version of h5py library
            random_seed (int) : Numpy random seed
            griddata_interpolation (bool) : Should we interpolate density data (5x slower)

        """

        self.rc = rc
        self.d_view = d_view

        if self.rc and not self.d_view:
            self.d_view = self.rc.direct_view(targets='all')

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
        self.config = None
        self.morphologies_stay_inside = morphologies_stay_inside

        self.axon_config_cache = None

        self.snudda_data = get_snudda_data(snudda_data=snudda_data,
                                           config_file=self.config_file,
                                           network_path=self.network_path)

        self.verbose = verbose
        self.log_file = log_file

        if self.snudda_data is None:
            self.write_log(f"Warning, snudda_data is not set!")

        if self.network_path:
            self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        else:
            self.write_log("No network_path given, not setting position_file. Remember to pass it to write_data.")
            self.position_file = None

        if h5libver is None:
            self.h5libver = "latest"
        else:
            self.h5libver = h5libver

        self.write_log(f"Using hdf5 version: {self.h5libver}")

        self.griddata_interpolation = griddata_interpolation

        # List of all neurons
        self.neurons = []
        self.neuron_prototypes = {}
        self.random_seed = random_seed
        self.random_generator = None
        self.rotate_helper = None

        # This defines the neuron units/channels. The dictionary lists all the
        # members of each unit, the neuronChannel gives the individual neurons
        # channel membership
        self.num_population_units = 1
        self.population_unit_placement_method = "random"
        self.population_units = dict([])
        self.population_unit = None

        # These are the dimensions of our space, d_min also needs a "padding"
        # region outside the space where fake neurons are placed. This is to
        # avoid boundary effects, without padding we would get too high density
        # at the edges
        self.volume = dict([])

    def __del__(self):

        if self.rc:
            # Cleanup memory on workers
            from snudda.utils import cleanup
            cleanup(self.rc, "place")

    ############################################################################

    def place(self):

        """ Place neurons in 3D space. """

        self.parse_config()

        if self.morphologies_stay_inside:
            self.avoid_edges_parallel()

        self.write_data()

    ############################################################################

    def write_log(self, text, flush=True, is_error=False, force_print=False):  # Change flush to False in future, debug

        """
        Writes to log file. Use setup_log first. Text is only written to screen if self.verbose=True,
        or is_error = True, or force_print = True.

        test (str) : Text to write
        flush (bool) : Should all writes be flushed to disk directly?
        is_error (bool) : Is this an error, always written.
        force_print (bool) : Force printing, even if self.verbose=False.
        """

        if self.log_file is not None:
            self.log_file.write(f"{text}\n")
            if flush:
                self.log_file.flush()

        if self.verbose or is_error or force_print:
            print(text, flush=True)

    ############################################################################

    def add_neurons(self,
                    swc_path,
                    num_neurons,
                    param_filename=None,
                    mech_filename=None,
                    modulation=None,
                    reaction_diffusion=None,
                    name="Unnamed",
                    hoc=None,
                    volume_id=None,
                    virtual_neuron=False,
                    axon_density=None,
                    parameter_key=None,
                    morphology_key=None,
                    modulation_key=None,
                    config=None,
                    rng=None):

        """
        Add neurons to volume specified.

        Args:
            swc_path (str): Path to morphology directory (or single morphology)
            num_neurons (int): Number of neurons to add
            param_filename (str): Path to parameter file
            mech_filename (str): Path to mechanism file
            modulation (str): Path to neuromodulation file
            reaction_diffusion (str): Path to RxD reaction diffusion file
            name (str): Name of neuron population, e.g. DSPN (which will become DSPN_0, DSPN_1, etc...)
            hoc (str): Path to hoc file (currently disabled)
            volume_id (str): ID of the volume to place neurons in
            virtual_neuron (bool): Real or virtual neuron, the latter can be used to model axons giving input to network
            axon_density (str): Axon density

            parameter_key (str, optional): Parameter Key to use, default=None (Randomise between all available)
            morphology_key (str, optional): Morphology Key to use, default=None (Randomise between all available)
            modulation_key (str, optional): Key for neuromodulation, default=None (Randomise between all available)
            config (dict): dict representation of network-config.json
        """

        assert volume_id is not None, f"You must specify a volume for neuron {name}"
        assert hoc is None, "Currently only support hoc=None, since we can have multiple parameter, morph combos now"

        if rng is None:
            rng = self.random_generator

        neuron_prototype = NeuronPrototype(neuron_name=name,
                                           neuron_path=None,
                                           snudda_data=self.snudda_data,
                                           morphology_path=swc_path,
                                           parameter_path=param_filename,
                                           mechanism_path=mech_filename,
                                           modulation_path=modulation,
                                           reaction_diffusion_path=reaction_diffusion,
                                           load_morphology=False,
                                           virtual_neuron=virtual_neuron,
                                           verbose=self.verbose)

        neuron_type = name.split("_")[0]
        neuron_positions = self.volume[volume_id]["mesh"].place_neurons(num_neurons, neuron_type)

        first_added = True

        neuron_rotations = self.rotate_helper.get_rotations(volume_name=volume_id, neuron_type=neuron_type,
                                                            neuron_positions=neuron_positions,
                                                            rng=rng)

        # Are there any projections from this neuron type in the config file?
        axon_info = self.generate_extra_axon_info(source_neuron=name,
                                                  position=neuron_positions,
                                                  config=config,
                                                  rng=rng)

        for idx, (coords, rotation) in enumerate(zip(neuron_positions, neuron_rotations)):

            # We set loadMorphology = False, to preserve memory
            # Only morphology loaded for nm then, to get axon and dend
            # radius needed for connectivity

            # Pick a random parameter set
            # parameter.json can be a list of lists, this allows you to select the
            # parameter set randomly
            # modulation.json is similarly formatted, pick a parameter set here
            if parameter_key is None:
                parameter_id = rng.integers(1000000)
            else:
                parameter_id = None

            if modulation_key is None:
                modulation_id = rng.integers(1000000)
            else:
                modulation_id = None

            if morphology_key is None:
                morphology_id = rng.integers(1000000)
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

            n.neuron_id = len(self.neurons)
            n.volume_id = volume_id

            n.axon_density = axon_density

            if axon_info is not None and len(axon_info) > 0:
                for axon_name, axon_position, axon_rotation, axon_swc in axon_info:
                    n.add_morphology(swc_file=axon_swc[idx],
                                     name=axon_name,
                                     position=axon_position[idx, :],
                                     rotation=axon_rotation[idx, :].reshape((3, 3)),
                                     lazy_loading=True)

            self.neurons.append(n)

            # This info is used by workers to speed things up
            if first_added:
                first_added = False
                self.neuron_prototypes[n.name] = n

    ############################################################################

    def parse_config(self, config_file=None, resort_neurons=True):

        """ Pase network config_file"""

        if config_file is None:
            config_file = self.config_file

        if config_file is None:
            self.write_log("No config file specified.", is_error=True)
            raise ValueError("No config file specified.")

        if not os.path.isfile(config_file):
            self.write_log(f"Missing config file {config_file}", is_error=True)
            raise ValueError(f"Missing config file {config_file}")

        self.write_log(f"Parsing place config file {config_file}")

        with open(config_file, "r") as f:
            config = json.load(f)

        if not config:
            self.write_log("Warning, empty network config.")

        if self.random_seed is None:
            if "random_seed" in config and "place" in config["random_seed"]:
                self.random_seed = config["random_seed"]["place"]
                self.write_log(f"Reading random seed from config file: {self.random_seed}")
        else:
            self.write_log(f"Using random seed provided by command line: {self.random_seed}")

        # Setup for rotations
        self.rotate_helper = SnuddaRotate(self.config_file)

        if self.random_seed:
            self.random_generator = np.random.default_rng(self.random_seed + 115)
        else:
            self.random_generator = np.random.default_rng()

        # Setup random seeds for all regions
        ss = np.random.SeedSequence(self.random_seed)
        all_seeds = ss.generate_state(len(config["regions"]))

        for (region_name, region_data), region_seed in zip(config["regions"].items(), all_seeds):
            total_num_neurons = region_data.get("num_neurons")

            volume_data = region_data["volume"]
            self.volume[region_name] = volume_data

            if "random_seed" not in volume_data:
                self.volume[region_name]["random_seed"] = region_seed
            else:
                region_seed = self.volume[region_name]["random_seed"]

            region_rnd = np.random.default_rng(region_seed + 123)

            if snudda_path_exists(volume_data["mesh_file"], self.snudda_data):
                mesh_file = snudda_parse_path(volume_data["mesh_file"], self.snudda_data)
            elif os.path.isfile(os.path.join(self.network_path, volume_data["mesh_file"])):
                mesh_file = os.path.join(self.network_path, volume_data["mesh_file"])
            else:
                raise ValueError(f"Unable to find mesh_file {volume_data['mesh_file']}")

            if "num_putative_points" in volume_data:
                num_putative_points = int(volume_data["num_putative_points"])
            else:
                num_putative_points = None

            self.volume[region_name]["mesh"] \
                = NeuronPlacer(mesh_path=mesh_file,
                               d_min=self.volume[region_name]["d_min"],
                               random_seed=region_seed,
                               n_putative_points=num_putative_points)

            if "density" in self.volume[region_name]:
                for neuron_type in self.volume[region_name]["density"]:
                    density_func = None

                    if "density_function" in self.volume[region_name]["density"][neuron_type]:
                        density_str = self.volume[region_name]["density"][neuron_type]["density_function"]
                        density_func = lambda x, y, z, d_str=density_str: numexpr.evaluate(d_str)

                    if "density_file" in self.volume[region_name]["density"][neuron_type]:
                        density_file = self.volume[region_name]["density"][neuron_type]["density_file"]

                        # We need to load the data from the file
                        from scipy.interpolate import griddata
                        with open(snudda_parse_path(density_file, self.snudda_data), "r") as f:
                            density_data = json.load(f)

                            assert region_name in density_data and neuron_type in density_data[region_name], \
                                f"Volume {region_name} does not contain data for neuron type {neuron_type}"

                            assert "coordinates" in density_data[region_name][neuron_type] \
                                   and "density" in density_data[region_name][neuron_type], \
                                (f"Missing coordinates and/or Density data for "
                                 f"volume {region_name}, neuron type {neuron_type}")

                            # Convert to SI (* 1e-6)
                            coord = np.array(density_data[region_name][neuron_type]["coordinates"]) * 1e-6
                            density = np.array(density_data[region_name][neuron_type]["density"])

                            if self.griddata_interpolation:
                                density_func = lambda x, y, z, c=coord, d=density: \
                                    griddata(points=c, values=d,
                                             xi=np.array([x, y, z]), method="linear",
                                             fill_value=0).transpose()
                            else:
                                density_func = lambda x, y, z, c=coord, d=density: \
                                    griddata(points=c, values=d,
                                             xi=np.array([x, y, z]), method="nearest",
                                             fill_value=0).transpose()

                    self.volume[region_name]["mesh"].define_density(neuron_type, density_func)

            if "neurons" not in region_data:
                self.write_log(f"No neurons specified for volume {region_name}")

            number_of_added_neurons = 0

            for neuron_type, neuron_data in region_data["neurons"].items():

                model_type = neuron_data.get("neuron_type", "neuron")  # default "neuron"
                # rotation_mode currently not used?!

                axon_density = neuron_data.get("axon_density")

                if "num_neurons" in neuron_data:
                    num_neurons = neuron_data["num_neurons"]
                elif "fraction" in neuron_data:
                    if total_num_neurons is None:
                        raise ValueError(f"If fraction is specified, then total_num_neurons for the region {region_name} must be set")

                    if neuron_data["fraction"] < 0 or neuron_data["fraction"] > 1:
                        raise ValueError(f"{neuron_type}: Neuron 'fraction' must be between 0 and 1.")
                    num_neurons = int(neuron_data["fraction"] * total_num_neurons)
                else:
                    raise ValueError(f"You need to specify 'fraction' or 'num_neurons' for {neuron_type}")

                # print(f"{neuron_type = }, {num_neurons = }")

                n_types = len(neuron_data["neuron_path"])

                if isinstance(num_neurons, int):
                    n_neurons = np.full((n_types, ), int(num_neurons/n_types))

                    extra_n = region_rnd.integers(low=0, high=n_types, size=num_neurons-np.sum(n_neurons))

                    for en in extra_n:
                        n_neurons[en] += 1
                elif len(num_neurons) != n_types:
                    raise ValueError(f"num_neurons can be a scalar or a vector, if it is a vector "
                                     f"({num_neurons = }) then its length must equal number of "
                                     f"neuron_path:s given ({n_types})")
                else:
                    n_neurons = num_neurons

                # RxD reaction diffusion config file
                default_reaction_diffusion = neuron_data.get("reaction_diffusion")

                parameter_key_list = SnuddaPlace.replicate_str(neuron_data.get("parameter_key"),
                                                               n_neurons, f"{neuron_type} parameter_key")
                morphology_key_list = SnuddaPlace.replicate_str(neuron_data.get("morphology_key"),
                                                                n_neurons, f"{neuron_type} morphology_key")
                modulation_key_list = SnuddaPlace.replicate_str(neuron_data.get("modulation_key"),
                                                                n_neurons, f"{neuron_type} modulation_key")

                for (neuron_name, neuron_path), num, parameter_key, morphology_key, modulation_key \
                        in zip(neuron_data["neuron_path"].items(), n_neurons,
                               parameter_key_list, morphology_key_list, modulation_key_list):

                    print(f"{neuron_name = }, {num = }, {neuron_path = }")

                    if neuron_name.split("_")[0] != neuron_type and neuron_name != neuron_type:
                        raise ValueError(f"The keys in neuron_path must be {neuron_name}_X where X is usually a number")

                    morph = os.path.join(neuron_path, "morphology")
                    param = os.path.join(neuron_path, "parameters.json")
                    mech = os.path.join(neuron_path, "mechanisms.json")

                    modulation = os.path.join(neuron_path, "modulation.json")
                    if not snudda_path_exists(modulation, snudda_data=self.snudda_data):
                        modulation = None

                    if default_reaction_diffusion is None:
                        reaction_diffusion = os.path.join(neuron_path, "reaction_diffusion.json")
                        if not snudda_path_exists(reaction_diffusion, snudda_data=self.snudda_data):
                            reaction_diffusion = None
                    else:
                        reaction_diffusion = default_reaction_diffusion

                    if model_type == "virtual":
                        param = None
                        mech = None
                        modulation = None
                        hoc = None
                        virtual_neuron = True
                    else:
                        virtual_neuron = False

                    self.add_neurons(name=neuron_name,
                                     swc_path=morph,
                                     param_filename=param,
                                     mech_filename=mech,
                                     modulation=modulation,
                                     reaction_diffusion=reaction_diffusion,
                                     num_neurons=num,
                                     hoc=None,
                                     volume_id=region_name,
                                     virtual_neuron=virtual_neuron,
                                     axon_density=axon_density,
                                     parameter_key=parameter_key,
                                     morphology_key=morphology_key,
                                     modulation_key=modulation_key,
                                     config=config,
                                     rng=region_rnd)

                    number_of_added_neurons += num

            if total_num_neurons is not None and number_of_added_neurons > total_num_neurons*1.01:
                raise ValueError(f"{region_name} should have {total_num_neurons} but {number_of_added_neurons} were added.")

        self.config_file = config_file
        self.config = config

        # We reorder neurons, sorting their IDs after position
        # -- UPDATE: Now we spatial cluster neurons depending on number of workers
        if resort_neurons:
            self.sort_neurons(sort_idx=self.cluster_neurons(rng=region_rnd))

        if False:  # Debug purposes, make sure neuron ranges are ok
            self.plot_ranges()

        self.define_population_units(config)

    @staticmethod
    def get_var_helper(data, key, default_value=None):
        return data[key] if key in data else default_value

    ############################################################################

    def avoid_edges_parallel(self):

        ss = np.random.SeedSequence(self.random_seed + 100)
        neuron_random_seed = ss.generate_state(len(self.neurons))

        bend_neuron_info = []

        for neuron in self.neurons:
            neuron_type = neuron.name.split("_")[0]
            neuron_config = self.config["regions"][neuron.volume_id]["neurons"][neuron_type]

            if "stay_inside_mesh" in neuron_config and neuron_config["stay_inside_mesh"]:
                volume_id = neuron_config["volume_id"]
                mesh_file = self.config["regions"][volume_id]["volume"]["mesh_file"]

                if isinstance(neuron_config["stay_inside_mesh"], dict):
                    if "k_dist" in neuron_config["stay_inside_mesh"]:
                        k_dist = neuron_config["stay_inside_mesh"]["k_dist"]
                    else:
                        k_dist = 30e-6

                    if "n_random" in neuron_config["stay_inside_mesh"]:
                        n_random = neuron_config["stay_inside_mesh"]["n_random"]
                    else:
                        n_random = 5

                    if "max_angle" in neuron_config["stay_inside_mesh"]:
                        max_angle = neuron_config["stay_inside_mesh"]["max_angle"]
                    else:
                        max_angle = 0.1  # radians
                else:
                    k_dist = 30e-6
                    n_random = 5
                    max_angle = 0.1  # radians

                bend_neuron_info.append((neuron.neuron_id, neuron.name, neuron.swc_filename,
                                         neuron.position, neuron.rotation,
                                         neuron_random_seed[neuron.neuron_id],
                                         volume_id, mesh_file, k_dist, n_random, max_angle))

        bend_morph_path = os.path.join(self.network_path, "modified_morphologies")

        if not os.path.isdir(bend_morph_path):
            os.mkdir(bend_morph_path)

        if self.d_view is None:
            # Make sure we use the same random seeds if we run in serial, as would have been used in parallel

            modified_neurons = self.avoid_edges_helper(bend_neuron_info=bend_neuron_info, network_path=self.network_path)

        else:

            # Make random permutation of neurons, to spread out the edge neurons
            unsorted_neuron_id = self.random_generator.permutation(len(bend_neuron_info))
            bend_neuron_info = [bend_neuron_info[idx] for idx in unsorted_neuron_id]

            with self.d_view.sync_imports():
                from snudda.place import SnuddaPlace

            self.d_view.scatter("bend_neuron_info", bend_neuron_info, block=True)
            self.d_view.push({"config_file": self.config_file,
                              "network_path": self.network_path,
                              "snudda_data": self.snudda_data},
                             block=True)

            cmd_str = f"sp = SnuddaPlace(config_file=config_file,network_path=network_path,snudda_data=snudda_data)"
            self.d_view.execute(cmd_str, block=True)

            cmd_str3 = f"modified_neurons = SnuddaPlace.avoid_edges_helper(bend_neuron_info=bend_neuron_info, network_path=network_path)"
            self.d_view.execute(cmd_str3, block=True)

            modified_neurons = self.d_view.gather("modified_neurons", block=True)

        for neuron_id, new_morphology in modified_neurons:
            # Replace the original morphology with the warped morphology, morphology includes rotation
            self.neurons[neuron_id].swc_filename = new_morphology
            self.neurons[neuron_id].rotation = np.eye(3)

    @staticmethod
    def avoid_edges_helper(bend_neuron_info, network_path):

        # TODO: We need name, swc_file, position, rotation
        # This needs to be passed, since self.neurons is not pickleable...

        from snudda.place.bend_morphologies import BendMorphologies

        bend_morph = dict()
        bend_morph_path = os.path.join(network_path, "modified_morphologies")

        modified_morphologies = []

        for neuron_id, neuron_name, swc_filename, position, rotation, random_seed, volume_id, mesh_file,\
            k_dist, n_random, max_angle in bend_neuron_info:

            if volume_id not in bend_morph:
                bend_morph[volume_id] = BendMorphologies(region_mesh=mesh_file, rng=None)

            # Returns None if unchanged
            new_morph_name = os.path.join(bend_morph_path, f"{neuron_name}-{neuron_id}.swc")
            new_morphology = bend_morph[volume_id].edge_avoiding_morphology(swc_file=swc_filename,
                                                                            new_file=new_morph_name,
                                                                            original_position=position,
                                                                            original_rotation=rotation,
                                                                            random_seed=random_seed,
                                                                            k_dist=k_dist,
                                                                            n_random=n_random,
                                                                            max_angle=max_angle)

            if new_morphology:
                modified_morphologies.append((neuron_id, new_morphology))

        return modified_morphologies

    ############################################################################

    def generate_extra_axon_info(self, source_neuron, position, config, rng):

        axon_info = []

        if self.axon_config_cache is None:

            self.axon_config_cache = dict()

            for region in config["regions"]:
                for neuron_key, neuron_info in config["regions"][region]["neurons"].items():
                    if "axon_config" in neuron_info:
                        axon_config = snudda_parse_path(neuron_info["axon_config"], self.snudda_data)
                        with open(axon_config, "r") as f:
                            self.axon_config_cache[neuron_key] = json.load(f)

        # Each neuron type has a set of axons to choose from
        source_neuron_type = source_neuron.split("_")[0]

        if source_neuron_type in self.axon_config_cache:
            axon_config = self.axon_config_cache[source_neuron_type]

            for axon_name, axon_data in axon_config.items():

                axon_position, axon_rotation, axon_swc \
                    = self.get_projection_axon_location(source_position=position,
                                                        proj_info=axon_data,
                                                        rng=rng)
                axon_info.append([axon_name, axon_position, axon_rotation, axon_swc])

            return axon_info

        else:
            # No extra axons for neuron
            return []

    #########################################################################

    def generate_extra_axon_info_old(self, source_neuron, position, config, rng):

        axon_info = []

        if self.axon_config_cache is None:

            self.axon_config_cache = dict()

            for neuron_key, neuron_info in config["neurons"].items():
                if "axon_config" in neuron_info:
                    axon_config = snudda_parse_path(neuron_info["axon_config"], self.snudda_data)
                    with open(axon_config, "r") as f:
                        self.axon_config_cache[neuron_key] = json.load(f)

        if source_neuron in self.axon_config_cache:
            axon_config = self.axon_config_cache[source_neuron]

            for axon_name, axon_data in axon_config.items():

                axon_position, axon_rotation, axon_swc \
                    = self.get_projection_axon_location(source_position=position,
                                                        proj_info=axon_data,
                                                        rng=rng)
                axon_info.append([axon_name, axon_position, axon_rotation, axon_swc])

            return axon_info

        else:
            # No extra axons for neuron
            return []

    def get_projection_axon_location(self, source_position, proj_info, rng, patch_hull=True):
        if 'projection' not in proj_info:
            raise KeyError("No 'projection' entry in the projection config!")
        
        proj_cfg = proj_info["projection"]
        if "projection_file" in proj_cfg and \
            ("source" in proj_cfg or \
           "destination" in proj_cfg):
           raise NotImplementedError("Projections should specify either a file or a mapping!")

        elif "source" in proj_cfg \
            and "destination" in proj_cfg:
            source = np.array(proj_cfg["source"])*1e-6
            destination = np.array(proj_cfg["destination"])*1e-6

        elif "projection_file" in proj_cfg :
            with open(proj_cfg["projection_file"], 'r') as f:
                proj_file_data = json.load(f)
            source = np.array(proj_file_data["source"])*1e-6
            destination = np.array(proj_file_data["destination"])*1e-6
            
        else: 
            raise NotImplementedError("Unknown projection configuration!")

        # specify the rotations of the termination zones
        rotation_cfg = proj_info.get("rotation", {})
        if "rotation" in rotation_cfg:
            rotation = np.array(rotation_cfg["rotation"])
            rot_position = destination
        
        elif "rotation_file" in rotation_cfg:
            with open(rotation_cfg["rotation_file"], 'r') as f:
                rotation_data = json.load(f)
            rotation = np.array(rotation_data["rotation"])
            rot_position = np.array(rotation_data["position"])*1e-6
        else:
            rotation = None
            rot_position = None

        target_centres = griddata(points=source,
                                  values=destination,
                                  xi=source_position,
                                  method="linear")

        # this checks if there are values outside of the convex hull
        to_patch = np.where(np.isnan(np.sum(target_centres, axis=1)))[0]

        # and patch missing entries
        if patch_hull and len(to_patch)>0:
            self.write_log(f"Patched {len(to_patch)}/{len(target_centres)}")
            target_centres_patched = griddata(points=source,
                                              values=destination,
                                              xi=source_position,
                                              method="nearest")
            target_centres[to_patch] = target_centres_patched[to_patch]

        # which coordinates to use for selecting rotation
        mapping = rotation_cfg.get('mapping', 'target')
        if mapping == "target":
            xi = target_centres
        elif mapping == "source":
            xi = source_position
        else:
            raise NotImplementedError(f"Unknown mapping '{mapping}'!")

        if rotation is not None:
            target_rotation = griddata(points=rot_position,
                                       values=rotation,
                                       xi=xi, method="linear")
            # if the rotation is specified as a field of rotation vectors, 
            # then these need to be converted to matrices.
            if target_rotation[0].shape == (3,):
                rotation_matrices = \
                [SnuddaRotate.rotation_matrix_from_vectors(np.array([0, 0, 1]), rv).flatten()\
                 for rv in target_rotation] 
                target_rotation = np.array(rotation_matrices)
                
        else:
            target_rotation = [None for x in range(source_position.shape[0])]

        num_axons = len(proj_info["morphologies"])
        axon_id = rng.choice(num_axons, source_position.shape[0])
        axon_swc = [proj_info["morphologies"][x] for x in axon_id]

        return target_centres, target_rotation, axon_swc

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

    def gather_extra_axons(self):

        ax_neuron = []
        ax_name = []
        ax_swc = []
        ax_position = []
        ax_rotation = []

        for n_idx, n in enumerate(self.neurons):
            for md_key, md in n.morphology_data.items():
                if md_key != "neuron":
                    ax_neuron.append(n_idx)  # which neuron does axon belong to
                    ax_name.append(md_key)
                    ax_swc.append(md.swc_file)
                    ax_position.append(md.position)
                    ax_rotation.append(md.rotation.reshape((1, 9)))

        if len(ax_neuron) > 0:
            ax_neuron = np.array(ax_neuron)
            ax_position = np.vstack(ax_position)
            ax_rotation = np.vstack(ax_rotation)

        return ax_neuron, ax_name, ax_position, ax_rotation, ax_swc

    ############################################################################

    def write_data(self, file_name=None):

        """ Writes position data to HDF5 file file_name. """

        if not file_name:
            file_name = self.position_file

        assert len(self.neurons) > 0, "No neurons to save!"

        self.write_log(f"Writing data to HDF5 file: {file_name}")

        pos_file = h5py.File(file_name, "w", libver=self.h5libver)

        with open(self.config_file, 'r') as cfg_file:
            config = json.load(cfg_file)

        # Meta data
        save_meta_data = [(self.config_file, "config_file"),
                          (json.dumps(config), "config"),
                          (snudda_parse_path("$DATA", self.snudda_data), "snudda_data")]

        meta_group = pos_file.create_group("meta")

        for data, data_name in save_meta_data:
            meta_group.create_dataset(data_name, data=data)

        network_group = pos_file.create_group("network")

        # Neuron information
        neuron_group = network_group.create_group("neurons")

        # If the name list is longer than 20 chars, increase S20
        name_list = [n.name.encode("ascii", "ignore") for n in self.neurons]
        str_type = 'S' + str(max(1, max([len(x) for x in name_list])))
        neuron_group.create_dataset("name", (len(name_list),), str_type, name_list,
                                    compression="gzip")

        neuron_id_list = np.arange(len(self.neurons))
        neuron_group.create_dataset("neuron_id", (len(neuron_id_list),),
                                    'int', neuron_id_list)

        volume_id_list = [n.volume_id.encode("ascii", "ignore")
                          for n in self.neurons]
        str_type_vid = 'S' + str(max(1, max([len(x) for x in volume_id_list])))

        neuron_group.create_dataset("volume_id",
                                    (len(volume_id_list),), str_type_vid, volume_id_list,
                                    compression="gzip")

        hoc_list = [snudda_simplify_path(n.hoc, self.snudda_data).encode("ascii", "ignore") if hasattr(n, "hoc") else "" for n in self.neurons]
        max_hoc_len = max([len(x) for x in hoc_list])
        max_hoc_len = max(max_hoc_len, 10)  # In case there are none
        neuron_group.create_dataset("hoc", (len(hoc_list),), f"S{max_hoc_len}", hoc_list,
                                    compression="gzip")

        swc_list = [snudda_simplify_path(n.swc_filename, self.snudda_data).encode("ascii", "ignore") for n in self.neurons]
        max_swc_len = max([len(x) for x in swc_list])

        neuron_group.create_dataset("morphology", (len(swc_list),), data=swc_list,
                                    dtype=h5py.special_dtype(vlen=bytes),
                                    compression="gzip")

        neuron_path = [snudda_simplify_path(n.neuron_path, self.snudda_data).encode("ascii", "ignore") for n in self.neurons]
        max_np_len = max([len(x) for x in neuron_path])
        # neuron_group.create_dataset("neuron_path", (len(neuron_path),), f"S{max_np_len}", neuron_path)
        neuron_group.create_dataset("neuron_path", (len(neuron_path),), data=neuron_path,
                                    dtype=h5py.special_dtype(vlen=str))

        virtual_neuron_list = np.array([n.virtual_neuron for n in self.neurons], dtype=bool)
        virtual_neuron = neuron_group.create_dataset("virtual_neuron",
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

        # Write axons to hdf5 file
        ax_neuron, ax_name, ax_position, ax_rotation, ax_swc = self.gather_extra_axons()

        axon_group = neuron_group.create_group("extra_axons")
        axon_group.create_dataset("parent_neuron", data=ax_neuron)
        axon_group.create_dataset("name", (len(ax_name), ), data=ax_name,
                                  dtype=h5py.special_dtype(vlen=bytes), compression="gzip")
        axon_group.create_dataset("position", data=ax_position)
        axon_group.create_dataset("rotation", data=ax_rotation)
        axon_group.create_dataset("morphology", (len(ax_swc), ), data=ax_swc,
                                  dtype=h5py.special_dtype(vlen=bytes), compression="gzip")

        # TODO: Parent tree info, eller motsvarande, måste sparas också!

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

        rd_list = [n.reaction_diffusion.encode("ascii", "ignore")
                   if n.reaction_diffusion is not None else ""
                   for n in self.neurons]
        rd_str_type = 'S' + str(max(1, max([len(x) for x in rd_list])))

        neuron_param_key = neuron_group.create_dataset("parameter_key",
                                                       (len(self.neurons),),
                                                       pk_str_type,
                                                       compression="gzip")

        neuron_morph_key = neuron_group.create_dataset("morphology_key",
                                                       (len(self.neurons),),
                                                       mk_str_type,
                                                       compression="gzip")

        neuron_modulation_key = neuron_group.create_dataset("modulation_key",
                                                            (len(self.neurons),),
                                                            mok_str_type,
                                                            compression="gzip")

        reaction_diffusion = neuron_group.create_dataset("reaction_diffusion_file",
                                                         (len(self.neurons),),
                                                         rd_str_type)

        neuron_pos_all = np.zeros((len(self.neurons), 3))
        neuron_rot_all = np.zeros((len(self.neurons), 9))

        neuron_param_key_list = []
        neuron_param_key_idx = []
        neuron_morph_key_list = []
        neuron_morph_key_idx = []
        neuron_mod_key_list = []
        neuron_mod_key_idx = []

        reacdiff_list = []
        reacdiff_key = []

        for (i, n) in enumerate(self.neurons):
            neuron_pos_all[i, :] = n.position
            neuron_rot_all[i, :] = n.rotation.reshape(1, 9)

            if n.parameter_key:
                neuron_param_key_list.append(n.parameter_key)
                neuron_param_key_idx.append(i)

            if n.morphology_key:
                neuron_morph_key_list.append(n.morphology_key)
                neuron_morph_key_idx.append(i)

            if n.modulation_key:
                neuron_mod_key_list.append(n.modulation_key)
                neuron_mod_key_idx.append(i)

            if n.reaction_diffusion:
                reacdiff_list.append(n.reaction_diffusion)
                reacdiff_key.append(i)

        neuron_position[:, :] = neuron_pos_all
        neuron_rotation[:, :] = neuron_rot_all

        if len(neuron_param_key_list) > 0:
            neuron_param_key[neuron_param_key_idx] = neuron_param_key_list

        if len(neuron_morph_key_list) > 0:
            neuron_morph_key[neuron_morph_key_idx] = neuron_morph_key_list

        if len(neuron_mod_key_list) > 0:
            neuron_modulation_key[neuron_mod_key_idx] = neuron_mod_key_list

        if len(reacdiff_list) > 0:
            reaction_diffusion[reacdiff_key] = reacdiff_list

        # Store input information
        if self.population_unit is None:
            # If no population units were defined, then set them all to 0 (= no population unit)
            self.population_unit = np.zeros((len(self.neurons),), dtype=int)

        neuron_group.create_dataset("population_unit_id", data=self.population_unit, dtype=int)

        # Variable for axon density "r", "xyz" or "" (No axon density)
        axon_density_type = [n.axon_density[0].encode("ascii", "ignore")
                             if n.axon_density is not None
                             else b""
                             for n in self.neurons]

        ad_str_type2 = "S" + str(max(1, max([len(x) if x is not None else 1
                                             for x in axon_density_type])))
        neuron_group.create_dataset("axon_density_type", (len(axon_density_type),),
                                    ad_str_type2, data=axon_density_type,
                                    compression="gzip")

        axon_density = [n.axon_density[1].encode("ascii", "ignore")
                        if n.axon_density is not None
                        else b""
                        for n in self.neurons]
        ad_str_type = "S" + str(max(1, max([len(x) if x is not None else 1
                                            for x in axon_density])))

        neuron_group.create_dataset("axon_density", (len(axon_density),),
                                    ad_str_type, data=axon_density, compression="gzip")

        axon_density_radius = [n.axon_density[2]
                               if n.axon_density is not None and n.axon_density[0] == "r"
                               else np.nan for n in self.neurons]

        neuron_group.create_dataset("axon_density_radius", data=axon_density_radius)

        # We also need to save axon_density_bounds_xyz, and num_axon points for the
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

        neuron_group.create_dataset("axon_density_bounds_xyz", data=axon_density_bounds_xyz)

        pos_file.close()

        self.write_log(f"Write done.")

    ############################################################################

    def define_population_units(self, config):

        """ Defines population units.

        """

        method_lookup = {"random": self.random_labeling,
                         "radial_density": self.population_unit_density_labeling,
                         "mesh": self.population_unit_mesh}

        for region_name in self.config["regions"]:
            if "population_units" in self.config["regions"][region_name]:
                for unit_id in self.config["regions"][region_name]["population_units"]["unit_id"]:
                    self.population_units[unit_id] = []

                neuron_id = self.volume_neurons(region_name)
                method_name = self.config["regions"][region_name]["population_units"]["method"]

                assert method_name in method_lookup, \
                    (f"Unknown population placement method {method_name}. "
                     f"Valid options are {', '.join([x for x in method_lookup])}")

                method_lookup[method_name](self.config["regions"][region_name]["population_units"], neuron_id)

    ############################################################################

    def random_labeling(self, population_unit_info, neuron_id):

        """
        Creates random labeling.

        Args:
             neuron_id (list) : All potential neuron ID
             population_unit_info (dict): Dictionary with "method" = "random"
                                          "unit_id" = ID of population unit
                                          "fraction_of_neurons" = Fraction of neurons in this unit
                                          "neuron_types" = List of neuron types that are in this unit
                                          "structure" = Name of the structure the unit is located in

        """

        self.init_population_units()  # This initialises population unit labelling if not already allocated

        unit_id = population_unit_info["unit_id"]
        fraction_of_neurons = population_unit_info["fraction_of_neurons"]
        neuron_types = population_unit_info["neuron_types"]  # list of neuron types that belong to this population unit
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
                # If larger than last cum_fraction, then no pop unit was picked (and we set pop unit 0)
                idx = len(cum_fraction) - np.sum(rn <= cum_fraction)

                if idx < len(cum_fraction):
                    uid = units_available[neuron_type]["unit"][idx]
                    self.population_unit[nid] = uid
                    self.population_units[uid].append(nid)
                else:
                    self.population_unit[nid] = 0

    ############################################################################

    def population_unit_density_labeling(self, population_unit_info, neuron_id):

        """
        Creates population units based on radial density functions.

        Args:
            neuron_id (list) : All potential neuron ID
            population_unit_info (dict): "method" must be "radial_density"
                                         "neuron_types" list of neuron types
                                         "centres" of radial probabilityes, one per neuron type
                                         "num_neurons" : Number of neurons in each population unit, only used with radial_density method
                                         "probability_functions" list of probability functions of r (as str)
                                         "unit_id" ID of population unit
        """

        assert population_unit_info["method"] == "radial_density"
        self.init_population_units()  # This initialises population unit labelling if not alraedy allocated

        neuron_types = population_unit_info["neuron_types"]
        centres = np.array(population_unit_info["centres"])
        probability_functions = population_unit_info["probability_functions"]
        unit_id = np.array(population_unit_info["unit_id"])

        if "num_neurons" in population_unit_info and population_unit_info["num_neurons"] is not None:
            num_neurons = population_unit_info["num_neurons"]

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
                    d = np.linalg.norm(pos - centre_pos)
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

        # Finally if num_neurons is specified, we need to reduce the number of neurons belonging to that unit
        for u_id, n_neurons in zip(unit_id, num_neurons):
            if n_neurons is not None:

                assert len(self.population_units[u_id]) >= n_neurons, \
                    f"Unable to pick {n_neurons} for population unit {u_id}, only {len(self.population_units[u_id])} ({neuron_types[np.where(unit_id == u_id)[0][0]]}) available."

                if n_neurons < len(self.population_units[u_id]):

                    perm_nid = np.random.permutation(self.population_units[u_id])
                    keep_nid = perm_nid[:n_neurons]
                    self.population_units[u_id] = keep_nid

                    remove_nid = perm_nid[n_neurons:]
                    for rid in remove_nid:
                        self.population_unit[rid] = 0

                    if 0 in self.population_units:
                        self.population_units[0] = np.sort(np.array(list(set(self.population_units[0]).union(remove_nid))))
                    else:
                        self.population_units[0] = remove_nid

    def population_unit_mesh(self, population_unit_info, neuron_id):

        self.init_population_units()  # This initialises population unit labelling if not already allocated

        unit_id = population_unit_info["unit_id"]
        mesh_file = population_unit_info["mesh_file"]
        fraction_of_neurons = population_unit_info["fraction_of_neurons"]
        neuron_types = population_unit_info["neuron_types"]  # list of neuron types that belong to this population unit
        structure_name = population_unit_info["structure"]

        pos = np.vstack([self.neurons[nid].position for nid in neuron_id])
        model_neuron_types = [self.neurons[nid].name.split("_")[0] for nid in neuron_id]

        member_probability = np.zeros(shape=(pos.shape[0], len(mesh_file)))

        for idx, (mf, frac, nts) in enumerate(zip(mesh_file, fraction_of_neurons, neuron_types)):

            # This checks if neurons are of the types that are included in population unit
            has_nt = np.array([n in nts for n in model_neuron_types], dtype=bool)

            rm = RegionMeshRedux(mf, verbose=self.verbose)
            member_probability[:, idx] = np.logical_and(rm.check_inside(pos), has_nt) * frac

        # If the probability sums to more than 1, then normalise it, otherwise keep smaller
        member_probability = np.divide(member_probability, np.maximum(1, np.sum(member_probability, axis =1).reshape(len(member_probability),1)))
        
        # Also, we need to add population unit 0 as an option, since choice needs P_sum = 1
        full_member_probability = np.zeros(shape=(member_probability.shape[0], member_probability.shape[1]+1))
        full_member_probability[:, :-1] = member_probability
        full_member_probability[:, -1] = np.maximum(1 - np.sum(member_probability, axis=1), 0)
        all_unit_id = unit_id + [0]

        # Normalise to 1
        row_sums = np.sum(full_member_probability, axis=1)
        row_sums = row_sums[:, np.newaxis]
        full_member_probability = full_member_probability / row_sums

        for idx, (nid, P) in enumerate(zip(neuron_id, full_member_probability)):

            uid = self.random_generator.choice(all_unit_id, p=P)
            self.population_unit[nid] = uid

            if uid > 0:
                self.population_units[uid].append(nid)

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
            colours[r_start:r_end] = idx + 1
            r_start = r_end

        xyz = self.all_neuron_positions()

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=colours, alpha=0.5)
        plt.show()

    def cluster_neurons(self, n_trials=10, rng=None):

        """
        Cluster neurons, so that nearby neurons are grouped on same worker, to speed up simulations.

        Args:
            n_trials (int) : Number of trials for k-means clustering (default 5)
        """

        if rng is None:
            rng = self.random_generator

        n_workers = len(self.d_view) if self.d_view is not None else 1
        n_clusters = np.maximum(n_workers * 5, 100)
        n_clusters = np.minimum(n_clusters, len(self.neurons))

        if n_workers > 1:
            self.write_log(f"Neurons order is optimised for {n_workers} workers. "
                           f"For reproducibility always use the same number of workers.")

        xyz = self.all_neuron_positions()

        centroids, labels = scipy.cluster.vq.kmeans2(xyz, n_clusters, minit="points", seed=rng)

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
                neuron_order[neuron_order_ctr:neuron_order_ctr + take_n] = cluster_member_list[current_cluster][:take_n]
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
            neuron_order = self.cluster_neurons(n_trials=n_trials - 1, rng=rng)

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

        for new_idx, old_idx in enumerate(sort_idx):
            self.neurons[old_idx].neuron_id = new_idx

        self.neurons = [self.neurons[x] for x in sort_idx]

        for idx, n in enumerate(self.neurons):
            assert idx == self.neurons[idx].neuron_id, \
                "Something went wrong with sorting"

    def volume_neurons(self, volume_id):

        return [n.neuron_id for n in self.neurons if n.volume_id == volume_id]

    @staticmethod
    def replicate_str(string, n_replicas, variable_name=None):

        # variable_name is just used to help with error message if incorrect input

        if string is None or isinstance(string, str):
            rep_str = [string for x in range(len(n_replicas))]
        elif len(string) == n_replicas:
            rep_str = string
        else:
            raise ValueError(f"Expected a str or a list of str of length {len(n_replicas)} for {variable_name}")

        return rep_str

    ############################################################################


if __name__ == "__main__":
    assert False, "Please use snudda.py place networks/yournetwork"
