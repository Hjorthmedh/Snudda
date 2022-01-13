# Rewriting the create network config file to make it more general

# !!! Currently writing full path for the files in the snudda data directory
#     this makes moving between computers difficult. Fix...

#
# Add a function so that $SNUDDADATA refers to the base datapath for snudda
#
import numpy as np
import numexpr
import sys
import os.path
import glob
import collections
from snudda.place.create_cube_mesh import create_cube_mesh
from snudda.place.create_slice_mesh import create_slice_mesh

import json
from snudda.utils.snudda_path import snudda_parse_path, snudda_simplify_path
from snudda.utils.snudda_path import snudda_path_exists
from snudda.utils.snudda_path import snudda_isdir
from snudda.utils.snudda_path import snudda_isfile
from snudda.utils.numpy_encoder import NumpyEncoder


class SnuddaInit(object):
    """ Creates network-config.json in network_path. """

    def __init__(self,
                 network_path=None,
                 struct_def=None,
                 neurons_dir=None,
                 config_file=None,
                 random_seed=None,
                 connection_override_file=None):

        """Constructor

           Args:
           network_path (str): location of network files
           struct_def (dict, optional): definition of struct to create
           neurons_dir (str, optional): path to neurons, default is $SNUDDA_DATA/neurons
           config_file (str, optional): name of network config file, default network-config.json
           random_seed (int, optional): random seed"""

        # print("CreateConfig")

        self.network_data = collections.OrderedDict([])

        self.network_data["RandomSeed"], self.init_rng = SnuddaInit.setup_random_seeds(random_seed)
        self.network_data["Volume"] = collections.OrderedDict([])
        self.num_neurons_total = 0

        if config_file:
            self.config_file = config_file
        elif network_path:
            self.config_file = os.path.join(network_path, "network-config.json")
        else:
            self.config_file = None

        if network_path:
            self.network_path = network_path
        elif config_file:
            self.network_path = os.path.dirname(config_file)
        else:
            self.network_path = ""

        if self.config_file and self.network_path:
            assert self.network_path == os.path.dirname(self.config_file), \
               f"network_path {self.network_path} and config_file path {self.config_file} must match"

        self.network_data["Connectivity"] = dict([])
        self.network_data["Neurons"] = dict([])

        struct_func = {"Striatum": self.define_striatum,
                       "GPe": self.define_GPe,
                       "GPi": self.define_GPi,
                       "STN": self.define_STN,
                       "SNr": self.define_SNr,
                       "Cortex": self.define_cortex,
                       "Thalamus": self.define_thalamus}

        if struct_def:
            for sn in struct_def:
                print("Adding " + sn + " with " + str(struct_def[sn]) + " neurons")
                struct_func[sn](num_neurons=struct_def[sn], neurons_dir=neurons_dir)

            if connection_override_file:
                self.replace_connectivity(connection_file=connection_override_file)

            # Only write JSON file if the structDef was not empty
            self.write_json(self.config_file)
        else:

            if connection_override_file:
                self.replace_connectivity(connection_file=connection_override_file)

            # print("No structDef defined, not writing JSON file in init")

    ############################################################################

    # meshBinWidth is used when voxelising the mesh to determine which part of
    # space is inside the mesh. For smaller structures we might need to use
    # a smaller meshBinWidth than the default 1e-4

    def define_structure(self,
                         struct_name,
                         struct_mesh,
                         d_min=None,
                         struct_centre=None,
                         side_len=None,
                         slice_depth=None,
                         mesh_bin_width=None):
        """
        Sets up definition for a brain structure (e.g. Cortex, Striatum, ...).

        Args:
            struct_name (str): Name of brain structure
            struct_mesh (str): Path to wavefront obj file with 3D mesh of structure or 'cube' or 'slice'
            d_min (float): Minimum distance between somas (puts upper limit on neuron density)
            struct_centre ((float, float, float)): Location of brain structure (centre)
            side_len (float, optional): side of cube, or slice
            slice_depth (float, optional): depth of slice
            mesh_bin_width (float): discretisation of 3D mesh during cell placement
        """

        if d_min is None:
            d_min = 15e-6

        if struct_mesh == "cube":
            assert slice_depth is None, "define_structure: sliceDepth is not used for cubes, please set to None"
            assert side_len is not None, "define_structure: cube needs sideLen specified"
            assert struct_centre is not None, "define_structure: cube needs a structCentre"

            struct_mesh = os.path.join(self.network_path, "mesh", f"{struct_name}-cube-mesh-{side_len}.obj")

            if mesh_bin_width is None:
                mesh_bin_width = side_len / 3.0
                print("Setting mesh_bin_width to " + str(mesh_bin_width))

            create_cube_mesh(file_name=struct_mesh,
                             centre_point=struct_centre,
                             side_len=side_len,
                             description=f"{struct_name} cube mesh, centre: {struct_centre}, side: {side_len}")

        elif struct_mesh == "slice":

            struct_mesh = os.path.join(self.network_path, "mesh", f"{struct_name}-slice-mesh-150mum-depth.obj")

            # 2019-11-26 : Anya said that her sagital striatal slices
            # were 2.36 x 2.36 mm. So that can be an upper limit

            if side_len is None:
                side_len = 200e-6

            if slice_depth is None:
                slice_depth = 150e-6

            print("Using slice depth: " + str(slice_depth))

            if mesh_bin_width is None:
                mesh_bin_width = np.minimum(side_len, slice_depth) / 3.0
                print("Setting meshBinWidth to " + str(mesh_bin_width))

            create_slice_mesh(file_name=struct_mesh,
                              centre_point=np.array([0, 0, 0]),
                              x_len=side_len,
                              y_len=side_len,
                              z_len=slice_depth,
                              description=struct_name + " slice mesh")

        if not snudda_path_exists(struct_mesh):
            print(f"Warning struct mesh {struct_mesh} is missing!")

        assert struct_name not in self.network_data["Volume"], \
            "defineStruct: Volume " + struct_name + " is already defined."

        self.network_data["Volume"][struct_name] = \
            self.define_volume(d_min=d_min,
                               mesh_file=struct_mesh,
                               mesh_bin_width=mesh_bin_width)

    ############################################################################

    @staticmethod
    def define_volume(mesh_file=None, d_min=15e-6, mesh_bin_width=1e-4):

        vol = dict([])
        vol["type"] = "mesh"
        vol["dMin"] = d_min
        vol["meshFile"] = mesh_file
        vol["meshBinWidth"] = mesh_bin_width

        return vol

    ############################################################################

    # This allows the user to specify a rotation field for neurons,
    # see examples/notebooks/example_of_neuronrotations.ipynb
    def define_rotation(self, volume_id, neuron_type, rotation_mode, rotation_field_file=None):

        if "neuronOrientation" not in self.network_data["Volume"][volume_id]:
            self.network_data["Volume"][volume_id]["neuronOrientation"] = collections.OrderedDict()

        self.network_data["Volume"][volume_id]["neuronOrientation"][neuron_type] = collections.OrderedDict()
        self.network_data["Volume"][volume_id]["neuronOrientation"][neuron_type]["rotationMode"] = rotation_mode

        if rotation_field_file:
            self.network_data["Volume"][volume_id]["neuronOrientation"][neuron_type]["rotationFieldFile"] \
                = rotation_field_file

    ############################################################################

    # conductance and conductanceStd -- allow variation of conductances
    #
    # channelParamDictionary = dictionary specifying other parameters, such as
    #                   fascilitation and depression of AMPA/NMDA channels etc

    def add_neuron_target(self, neuron_name, target_name, connection_type,
                          dist_pruning,
                          f1,
                          soft_max,
                          mu2,
                          a3,
                          dist_pruning_other=None,
                          f1_other=None,
                          soft_max_other=None,
                          mu2_other=None,
                          a3_other=None,
                          conductance=None,
                          cluster_synapses=False,
                          mod_file=None,
                          parameter_file=None,
                          channel_param_dictionary=None):

        if conductance is None:
            conductance = [1.0e-9, 0]

        # if(connectionType == "GapJunction"):
        #  assert f1 is None and softMax is None and mu2 is None and a3 is None\
        #    and f1_other is None and softMax_other is None and mu2_other is None \
        #    and a3_other is None, \
        #    "addNeuronTarget: " + str(neuronName) \
        #    + ", pruning not currently available for gap junctions"

        if parameter_file is not None:
            if channel_param_dictionary is None:
                channel_param_dictionary = dict([])

            channel_param_dictionary["parameterFile"] = parameter_file

        if mod_file is not None:
            if channel_param_dictionary is None:
                channel_param_dictionary = dict([])

            channel_param_dictionary["modFile"] = mod_file

        if type(conductance) == list:
            cond = conductance[0]
            cond_std = conductance[1]
        else:
            cond = conductance
            cond_std = 0

        con_info = dict([])
        con_info["conductance"] = [cond, cond_std]  # Mean, Std
        con_info["channelParameters"] = channel_param_dictionary
        pruning_info = dict([])
        pruning_info["f1"] = f1
        pruning_info["softMax"] = soft_max
        pruning_info["mu2"] = mu2
        pruning_info["a3"] = a3
        pruning_info["distPruning"] = dist_pruning
        pruning_info["cluster"] = cluster_synapses
        con_info["pruning"] = pruning_info

        # pruneInfo = (distPruning,f1,softMax,mu2,a3)

        if (dist_pruning_other is not None
                or f1_other is not None
                or soft_max_other is not None
                or mu2_other is not None
                or a3_other is not None):

            # If any of the other varibles is set,
            # then all other "other" variables that are
            # not set assume default values

            if dist_pruning_other is None:
                dist_pruning_other = dist_pruning

            if f1_other is None:
                f1_other = f1

            if soft_max_other is None:
                soft_max_other = soft_max

            if mu2_other is None:
                mu2_other = mu2

            if a3_other is None:
                a3_other = a3

            pruning_info_other = dict([])
            pruning_info_other["f1"] = f1_other
            pruning_info_other["softMax"] = soft_max_other
            pruning_info_other["mu2"] = mu2_other
            pruning_info_other["a3"] = a3_other
            pruning_info_other["distPruning"] = dist_pruning_other
            pruning_info_other["cluster"] = cluster_synapses

            # Different pruning rules for within and between neuron units
            con_info["pruningOther"] = pruning_info_other

        # Json did not like tuples in keys, so we separate by comma
        nt_key = f"{neuron_name},{target_name}"
        if nt_key not in self.network_data["Connectivity"]:
            self.network_data["Connectivity"][nt_key] = dict([])

        self.network_data["Connectivity"][nt_key][connection_type] = con_info

    ############################################################################

    def replace_connectivity(self, connection_file=None, connection_dict=None):

        """ Replaces the default connectivity.
        
        Args: 
            connection_file : Path to JSON file with connection block
            connection_dict : OrderedDict (or dict) with connection block
        
        """

        assert (connection_file is not None) + (connection_dict is not None) == 1, \
            f"replace_connectivity: One of connection_file and connection_dict should be given"

        if connection_file:
            connection_file = snudda_parse_path(connection_file)
            print(f"Reading connectivity from {connection_file}")
            assert os.path.isfile(connection_file), f"Connection JSON file {connection_file} does not exist."

            with open(connection_file, "r") as f:
                connection_dict = json.load(f, object_pairs_hook=collections.OrderedDict)

        assert type(connection_dict) in [collections.OrderedDict, dict]

        assert "Connectivity" in connection_dict, \
            f"replace_connectivity: Missing the Connectivity key in dictionary:\n{connection_dict}"

        self.network_data["Connectivity"] = connection_dict["Connectivity"]

    ############################################################################

    # modelType is "neuron" or "virtual" (= just provides input to network)
    # For axonDensity when it is "xyz" we assume that soma is at 0,0,0

    # neuronDir contains all the neurons in separate directories
    # Each of those directories have a config and morphology subdirectory

    def add_neurons(self, name,
                    neuron_dir,
                    num_neurons,
                    axon_density=None,
                    model_type="neuron",
                    volume_id=None,
                    rotation_mode="random"):

        if num_neurons <= 0:
            return

        if axon_density is not None:
            if axon_density[0] == "r":
                # Verify axon density function
                r = np.linspace(0, axon_density[2], 10)  # r is used in eval
                try:
                    numexpr.evaluate(axon_density[1])
                except:
                    print("!!! Axon density failed test: " + str(axon_density))
                    print("Inparameter: r = 1-D array of radius in meter")
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    sys.exit(-1)
            elif axon_density[0] == "xyz":
                x = np.linspace(axon_density[2][0], axon_density[2][1], 10)  # x,y,z used in eval below
                y = np.linspace(axon_density[2][2], axon_density[2][3], 10)
                z = np.linspace(axon_density[2][4], axon_density[2][5], 10)
                try:
                    numexpr.evaluate(axon_density[1])
                except:
                    print("!!! Axon density failed test: " + str(axon_density))
                    print("Inparameters: x,y,z three 1-D arrays (units in meter)")
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    sys.exit(-1)

                print("Checking boundaries, to make sure P is not too high")
                x = np.zeros((8, 1))
                y = np.zeros((8, 1))
                z = np.zeros((8, 1))
                ctr = 0

                for xx in axon_density[2][0:2]:
                    for yy in axon_density[2][2:4]:
                        for zz in axon_density[2][4:6]:
                            x[ctr] = xx
                            y[ctr] = yy
                            z[ctr] = zz
                            ctr += 1

                # axon_density function of x,y,z defined above
                p_corner = numexpr.evaluate(axon_density[1]) * (3e-6 ** 3)

                for P, xx, yy, zz in zip(p_corner, x, y, z):
                    print(f"{name} axon density P({xx}, {yy}, {zz}) = {P}")

                if (p_corner > 0.01).any():
                    print("Axon density too high at boundary!!")
                    print("Please increase bounding box")
                    sys.exit(-1)

                # print(str(axonDensity[3]) + " " + str(name) \
                #      + " axon points to place")

        print(f"Adding neurons: {name} from dir {snudda_parse_path(neuron_dir)}")
        # TODO: We should force users to use same name as the directory name
        # ie, fs/FS_0 directory should be named FS_0

        # Find which neurons are available in neuronDir
        dir_list = glob.glob(snudda_parse_path(neuron_dir) + "/*")
        neuron_file_list = []

        assert len(dir_list) > 0, f"Neuron dir {snudda_parse_path(neuron_dir)} is empty!"

        for fd in dir_list:

            d = snudda_simplify_path(fd)

            if snudda_isdir(d):
                # We want to maintain the $SNUDDA_DATA keyword in the path so that the user can move
                # the config file between systems and still run it.
                par_file = os.path.join(d, "parameters.json")
                mech_file = os.path.join(d, "mechanisms.json")
                modulation_file = os.path.join(d, "modulation.json")
                if not snudda_path_exists(modulation_file):
                    modulation_file = None

                sd = snudda_parse_path(d)
                swc_file = self.get_morphologies(sd)
                hoc_file = glob.glob(os.path.join(sd, "*hoc"))
                hoc_file = [snudda_simplify_path(hf) for hf in hoc_file]

                assert len(hoc_file) <= 1, f"Morph dir {sd} contains more than one hoc file"

                if len(hoc_file) == 0:
                    hoc_file = [None]

                neuron_file_list.append((d,
                                        swc_file,
                                        par_file,
                                        mech_file,
                                        modulation_file,
                                        hoc_file[0]))

        # First check how many unique cells we hava available, then we
        # calculate how many of each to use in simulation
        n_ind = len(neuron_file_list)
        assert n_ind > 0, \
            f"No swc morphologies found in {neuron_dir}.\nObs, each morphology should have its own subdirectory."

        n_of_each_ind = np.zeros((n_ind,))
        n_of_each_ind[:] = int(num_neurons / n_ind)
        still_to_add = int(num_neurons - np.sum(n_of_each_ind))
        add_idx = self.init_rng.permutation(n_ind)[0:still_to_add]
        n_of_each_ind[add_idx] += 1

        # Add the neurons to config

        for ctr, ((nrnDir, swc_file, par_file, mech_file, modulation_file, hoc_file), num) \
                in enumerate(zip(neuron_file_list, n_of_each_ind)):

            if int(num) == 0:
                continue

            unique_name = name + "_" + str(ctr)
            cell_data = dict([])

            if not snudda_isfile(par_file) and model_type != "virtual":
                print(f"Parameter file not found: {snudda_parse_path(par_file)}")

            if not snudda_isfile(mech_file) and model_type != "virtual":
                print(f"Mechanism file not found: {snudda_parse_path(mech_file)}")

            if hoc_file is not None and not snudda_isfile(hoc_file):
                print(f"Hoc file not found: {snudda_parse_path(hoc_file)}")

            cell_data["morphology"] = swc_file
            cell_data["parameters"] = par_file
            cell_data["mechanisms"] = mech_file

            if modulation_file is not None:
                # Modulation is optional
                cell_data["modulation"] = modulation_file

            cell_data["num"] = int(num)
            cell_data["hoc"] = hoc_file

            cell_data["neuronType"] = model_type
            cell_data["rotationMode"] = rotation_mode
            cell_data["volumeID"] = volume_id

            if axon_density is not None:
                cell_data["axonDensity"] = axon_density

            self.network_data["Neurons"][unique_name] = cell_data

    @staticmethod
    def get_morphologies(neuron_dir):
        """
        Returns SWC morphology(s) path or file, depending on 'morphology' if specified in parameters.json or not.

        If 'morphology' in parameters.json exists then the path where these morphologies are stored is returned.
        If it does not exist it is assumed that there is exactly one SWC file present in the neuron_dir,
        and that is then returned.

        Args:
            neuron_dir (str): Path to neuron directory, may contain $SNUDDA_DATA, shorthand for SNUDDA_DATA directory
        """

        parameter_file = os.path.join(neuron_dir, "parameters.json")

        if os.path.isfile(parameter_file):

            # First check if the morphologies are listed in the parameter file
            with open(parameter_file, "r") as f:
                par_data = json.load(f, object_pairs_hook=collections.OrderedDict)

            # We now expect a dictionary of parameter sets. If it is a list, we convert it to a dictionary
            if type(par_data) == list:
                par_data = {"default": par_data}

            meta_file = os.path.join(neuron_dir, "meta.json")
            if os.path.isfile(meta_file):
                with open(meta_file, "r") as mf:
                    meta_data = json.load(mf, object_pairs_hook=collections.OrderedDict)

                has_meta = True
            else:
                has_meta = False

        else:
            print("No parameter.json file.")
            has_meta = False

        if has_meta:

            morph_dir = os.path.join(neuron_dir, "morphology")
            morph_dir_full = snudda_parse_path(morph_dir)
            assert os.path.exists(morph_dir_full), \
                f"Morphology directory missing: {morph_dir_full}"

            # Also check that all morphologies listed exists
            missing_par_key = []
            missing_morphology_tag = []
            missing_morph = []
            for par_key in par_data.keys():
                if par_key not in meta_data:
                    missing_par_key.append(par_key)
                else:
                    for morph_key in meta_data[par_key].keys():
                        if "morphology" not in meta_data[par_key][morph_key]:
                            missing_morphology_tag.append((par_key, morph_key))
                        elif not os.path.isfile(os.path.join(morph_dir_full,
                                                             meta_data[par_key][morph_key]["morphology"])):
                            missing_morph.append(meta_data[par_key][morph_key]["morphology"])

            assert len(missing_par_key) == 0, \
                f"Missing parameter key(s) {', '.join(missing_par_key)} in {meta_file}"

            assert len(missing_morphology_tag) == 0, \
                f"Missing morphology tag(s) for {', '.join(missing_morphology_tag)} in {meta_file}"

            assert len(missing_morph) == 0, \
                f"The following morphologies in {meta_file} are missing: {', '.join(missing_morph)}"

            return snudda_simplify_path(morph_dir)

        else:
            swc_file = glob.glob(os.path.join(snudda_parse_path(neuron_dir), "*swc"))
            assert len(swc_file) == 1, \
                (f"If no morphology is given in parameter.json then " 
                 f"{snudda_parse_path(neuron_dir)} should contain exactly one swc file")

            swc_file = snudda_simplify_path(swc_file[0])
            return swc_file

    def add_neuron_density(self, volume_id, neuron_type, density_func=None, density_file=None):

        assert volume_id in self.network_data["Volume"], f"Volume {volume_id} not defined"

        assert (density_func is None) + (density_file is None) == 1, \
            f"Volume {volume_id}, neuron type {neuron_type}: Only one of density_func and density_file should be set"

        if "density" not in self.network_data["Volume"][volume_id]:
            self.network_data["Volume"][volume_id]["density"] = dict()

        self.network_data["Volume"][volume_id]["density"][neuron_type] = dict()

        if density_func:
            self.network_data["Volume"][volume_id]["density"][neuron_type]["densityFunction"] = density_func

        if density_file:
            self.network_data["Volume"][volume_id]["density"][neuron_type]["densityFile"] = density_file

        ############################################################################

    def write_json(self, filename=None):

        if not filename:
            filename = self.config_file

        # Create directory if it does not already exist
        dir_name = os.path.dirname(filename)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        print(f"Writing {filename}")

        import json
        with open(filename, 'w') as f:
            json.dump(self.network_data, f, indent=4, cls=NumpyEncoder)

    ###########################################################################################################

    def neuron_projection(self, neuron_name, target_name,
                          projection_name,
                          projection_file,
                          source_volume,
                          dest_volume,
                          projection_radius,
                          number_of_targets,
                          number_of_synapses,
                          dendrite_synapse_density,
                          connection_type,
                          dist_pruning=None,
                          f1=None,
                          soft_max=None,
                          mu2=None,
                          a3=None,
                          cluster_synapses=False,
                          dist_pruning_other=None,
                          f1_other=None,
                          soft_max_other=None,
                          mu2_other=None,
                          a3_other=None,
                          conductance=None,
                          mod_file=None,
                          parameter_file=None,
                          channel_param_dictionary=None):

        self.add_neuron_target(neuron_name, target_name,
                               connection_type,
                               dist_pruning,
                               f1,
                               soft_max,
                               mu2,
                               a3,
                               cluster_synapses=cluster_synapses,
                               dist_pruning_other=None,
                               f1_other=None,
                               soft_max_other=None,
                               mu2_other=None,
                               a3_other=None,
                               conductance=None,
                               mod_file=None,
                               parameter_file=None,
                               channel_param_dictionary=channel_param_dictionary)

        # Next we need to add the connection mapping specific parameters

        nt_key = f"{neuron_name},{target_name}"
        con_info = self.network_data["Connectivity"][nt_key][connection_type]

        con_info["projectionFile"] = projection_file
        con_info["projectionName"] = projection_name
        con_info["source"] = source_volume
        con_info["destination"] = dest_volume
        con_info["projectionRadius"] = projection_radius
        con_info["numberOfTargets"] = number_of_targets
        con_info["numberOfSynapses"] = number_of_synapses
        con_info["dendriteSynapseDensity"] = dendrite_synapse_density

    ############################################################################

    # Population Units here refer to processing units, where the neurons within a Population Unit
    # might have different connectivity than neurons belonging to different population Units

    # Centre of striatum mesh is [3540e-6,4645e-6,5081e-6]
    # Radius is now in SI units also (meters)

    def add_population_unit_density(self,
                                    structure_name,
                                    neuron_types,
                                    unit_centre,
                                    probability_function,  # Function of d (distance to centre) as string
                                    num_neurons=None,
                                    unit_id=None):

        if type(neuron_types) != list:
            neuron_types = list(neuron_types)

        unit_id = self.setup_population_unit(unit_id)

        if structure_name not in self.network_data["PopulationUnits"]:
            self.network_data["PopulationUnits"][structure_name] = collections.OrderedDict()
            self.network_data["PopulationUnits"][structure_name]["method"] = "radialDensity"

            self.network_data["PopulationUnits"][structure_name]["centres"] = [unit_centre]
            self.network_data["PopulationUnits"][structure_name]["ProbabilityFunctions"] = [probability_function]
            self.network_data["PopulationUnits"][structure_name]["unitID"] = [unit_id]
            self.network_data["PopulationUnits"][structure_name]["neuronTypes"] = [neuron_types]
            self.network_data["PopulationUnits"][structure_name]["numNeurons"] = [num_neurons]

        else:
            self.network_data["PopulationUnits"][structure_name]["centres"].append(unit_centre)
            self.network_data["PopulationUnits"][structure_name]["ProbabilityFunctions"].append(probability_function)
            self.network_data["PopulationUnits"][structure_name]["unitID"].append(unit_id)
            self.network_data["PopulationUnits"][structure_name]["neuronTypes"].append(neuron_types)
            self.network_data["PopulationUnits"][structure_name]["numNeurons"].append(num_neurons)


    def add_population_unit_random(self, structure_name, neuron_types, fraction_of_neurons, unit_id=None):

        if type(neuron_types) != list:
            neuron_types = list(neuron_types)

        unit_id = self.setup_population_unit(unit_id)

        if structure_name not in self.network_data["PopulationUnits"]:
            self.network_data["PopulationUnits"][structure_name] = collections.OrderedDict()
            self.network_data["PopulationUnits"][structure_name]["method"] = "random"

            self.network_data["PopulationUnits"][structure_name]["fractionOfNeurons"] = [fraction_of_neurons]
            self.network_data["PopulationUnits"][structure_name]["unitID"] = [unit_id]
            self.network_data["PopulationUnits"][structure_name]["neuronTypes"] = [neuron_types]
            self.network_data["PopulationUnits"][structure_name]["structure"] = structure_name
        else:
            self.network_data["PopulationUnits"][structure_name]["fractionOfNeurons"].append(fraction_of_neurons)
            self.network_data["PopulationUnits"][structure_name]["unitID"].append(unit_id)
            self.network_data["PopulationUnits"][structure_name]["neuronTypes"].append(neuron_types)

    # Helper function, returns next free unit_id and sets up data structures (user can also choose own unit_id
    # but it must be unique and not already used
    def setup_population_unit(self, unit_id=None):

        if "PopulationUnits" not in self.network_data:
            self.network_data["PopulationUnits"] = collections.OrderedDict([])
            self.network_data["PopulationUnits"]["AllUnitID"] = []

        if not unit_id:
            if self.network_data["PopulationUnits"]["AllUnitID"]:
                unit_id = np.max(self.network_data["PopulationUnits"]["AllUnitID"]) + 1
            else:
                unit_id = 1

        assert unit_id not in self.network_data["PopulationUnits"]["AllUnitID"], f"Unit id {unit_id} already in use"
        self.network_data["PopulationUnits"]["AllUnitID"].append(unit_id)

        return unit_id

    ############################################################################

    # Normally: nNeurons = number of neurons set, then the fractions specified
    # fMSD1, fMSD2, fFS, fChIN, fLTS are used to  calculate the number of neurons
    # of each type.
    #
    # If nNeurons is set to None, then we use nMSD1, nMSD2, nFS, nChIN, nLTS
    # to set the number of neurons. This is useful if you want to create a small
    # test network.

    # Divide by fTot since we are not including all neurons and we want the
    # proportions to sum to 1.0 (f means fraction)

    # mesh_file can be used to override default mesh file

    def define_striatum(self,
                        num_neurons=None,
                        f_dSPN=0.475,
                        f_iSPN=0.475,
                        f_FS=0.013,
                        f_ChIN=0.011,
                        f_LTS=0.007,
                        num_dSPN=None,
                        num_iSPN=None,
                        num_FS=None,
                        num_ChIN=None,
                        num_LTS=None,
                        volume_type=None,
                        side_len=None,
                        slice_depth=None,
                        neurons_dir=None,
                        neuron_density=80500,
                        within_population_unit_SPN_modifier=1,
                        between_population_unit_SPN_modifier=1,
                        mesh_file=None,
                        mesh_bin_width=None,
                        d_min=None,
                        cluster_FS_synapses=False,
                        cluster_SPN_synapses=False):

        get_val = lambda x: 0 if x is None else x
        if num_neurons is None:
            self.num_dSPN = get_val(num_dSPN)
            self.num_iSPN = get_val(num_iSPN)
            self.num_FS = get_val(num_FS)
            self.num_ChIN = get_val(num_ChIN)
            self.num_LTS = get_val(num_LTS)

            self.num_neurons_total += self.num_FS + self.num_dSPN + self.num_iSPN + self.num_ChIN + self.num_LTS
            num_neurons = self.num_neurons_total

            if self.num_neurons_total <= 0:
                # No neurons specified, skipping structure
                return
        else:
            if num_neurons <= 0:
                # No neurons specified, skipping structure
                return

            f_tot = f_dSPN + f_iSPN + f_FS + f_ChIN + f_LTS

            self.num_FS = np.round(f_FS * num_neurons / f_tot)
            self.num_dSPN = np.round(f_dSPN * num_neurons / f_tot)
            self.num_iSPN = np.round(f_iSPN * num_neurons / f_tot)
            self.num_ChIN = np.round(f_ChIN * num_neurons / f_tot)
            self.num_LTS = np.round(f_LTS * num_neurons / f_tot)

            self.num_neurons_total += self.num_FS + self.num_dSPN + self.num_iSPN + self.num_ChIN + self.num_LTS

            if abs(num_neurons - self.num_neurons_total) > 5:
                print("Striatum should have " + str(num_neurons) + " but " + str(self.num_neurons_total) \
                      + " are being requested, check fractions set for defineStriatum.")

        assert volume_type is None or mesh_file is None, "You should not specify both volume_type and mesh_file"

        if mesh_file:

            assert mesh_bin_width, "If you specify mesh_file you need to specify mesh_bin_width (e.g 1e-4)"

            self.define_structure(struct_name="Striatum",
                                  struct_mesh=mesh_file,
                                  mesh_bin_width=mesh_bin_width,
                                  d_min=d_min)

        elif volume_type == "mouseStriatum":
            if mesh_bin_width is None:
                mesh_bin_width = 1e-4

            self.define_structure(struct_name="Striatum",
                                  struct_mesh=os.path.join("$SNUDDA_DATA", "mesh", "Striatum-d.obj"),
                                  mesh_bin_width=mesh_bin_width,
                                  d_min=d_min)

            density_file = os.path.join("$SNUDDA_DATA", "density", "dorsal_striatum_density.json")

            self.add_neuron_density(volume_id="Striatum", neuron_type="dSPN", density_file=density_file)
            self.add_neuron_density(volume_id="Striatum", neuron_type="iSPN", density_file=density_file)
            self.add_neuron_density(volume_id="Striatum", neuron_type="FS", density_file=density_file)
            self.add_neuron_density(volume_id="Striatum", neuron_type="LTS", density_file=density_file)
            self.add_neuron_density(volume_id="Striatum", neuron_type="ChIN", density_file=density_file)

        elif volume_type == "slice":
            self.define_structure(struct_name="Striatum",
                                  struct_mesh="slice",
                                  side_len=side_len,
                                  slice_depth=slice_depth,
                                  d_min=d_min)

        elif num_neurons <= 1e6:  # 1e6
            print("Using cube for striatum")
            # 1.73 million neurons, volume of allen striatal mesh is 21.5mm3
            striatum_volume = 1e-9 * num_neurons / neuron_density  # 80.5e3
            striatum_side_len = striatum_volume ** (1. / 3)
            striatum_centre = np.array([4750e-6,4000e-6, 7750e-6])

            if num_neurons < 500:
                mesh_bin_width = striatum_side_len
            elif num_neurons < 5000:
                mesh_bin_width = striatum_side_len / 5
            else:
                mesh_bin_width = striatum_side_len / 10

            # Reduced striatum, due to few neurons
            self.define_structure(struct_name="Striatum",
                                  struct_mesh="cube",
                                  struct_centre=striatum_centre,
                                  side_len=striatum_side_len,
                                  mesh_bin_width=mesh_bin_width,
                                  d_min=d_min)

        else:
            # Default, full size striatum
            self.define_structure(struct_name="Striatum",
                                  struct_mesh=os.path.join("$SNUDDA_DATA", "mesh", "Striatum-d.obj"),
                                  mesh_bin_width=1e-4,
                                  d_min=d_min)

            density_file = os.path.join("$SNUDDA_DATA", "density", "dorsal_striatum_density.json")

            self.add_neuron_density(volume_id="Striatum", neuron_type="dSPN", density_file=density_file)
            self.add_neuron_density(volume_id="Striatum", neuron_type="iSPN", density_file=density_file)
            self.add_neuron_density(volume_id="Striatum", neuron_type="FS", density_file=density_file)
            self.add_neuron_density(volume_id="Striatum", neuron_type="LTS", density_file=density_file)
            self.add_neuron_density(volume_id="Striatum", neuron_type="ChIN", density_file=density_file)

        if neurons_dir is None:
            neurons_dir = os.path.join("$SNUDDA_DATA", "neurons")

        print(f"Neurons for striatum read from {snudda_parse_path(neurons_dir)}/striatum")

        FS_dir = os.path.join(neurons_dir, "striatum", "fs")
        dSPN_dir = os.path.join(neurons_dir, "striatum", "dspn")
        iSPN_dir = os.path.join(neurons_dir, "striatum", "ispn")
        ChIN_dir = os.path.join(neurons_dir, "striatum", "chin")
        LTS_dir = os.path.join(neurons_dir, "striatum", "lts")

        # Add the neurons

        self.add_neurons(name="FS", neuron_dir=FS_dir,
                         num_neurons=self.num_FS,
                         volume_id="Striatum")

        self.add_neurons(name="dSPN", neuron_dir=dSPN_dir,
                         num_neurons=self.num_dSPN,
                         volume_id="Striatum")

        self.add_neurons(name="iSPN", neuron_dir=iSPN_dir,
                         num_neurons=self.num_iSPN,
                         volume_id="Striatum")

        # ChIN axon density,
        # We start with the axon length per unit volume, then we scale it
        # to synapses per unit volume
        # This will then be used to precompute a lookup table
        # Guestimated density from Suzuki 2001, J Neurosci, figure 1bb
        # see directory morphology/ChINdensityEstimate

        # "I Janickova et al. 2017 så har de 2018 varicosities i en area på 655 um²,
        # deras slices är 70 um tjocka och om man antar att det inte är några
        # varicositites som täcker varandra så är volym-densiteten/mm³: 4.4*10⁷/mm3"
        # 1.7e6/24*0.01 = 708 ChIN per mm3
        # 4.4e7 / 708 = 62000 varicosities per ChIN
        #
        # 325 ChIN synapser per MS
        # 2-5 ChIN per MS
        # --> 65-160 synapser between a ChIN-MS pair
        # --> Each ChIN connect to 400 - 950 MS
        #
        # Number of MS within 350 micrometer radius 4*pi*(350e-6)^3/3*1.76e6/24e-9
        # --> 13100 MS reachable by ChIN at most (or rather number of MS somas
        # within radius of axonal arbour)
        # -->  3-7% connectivity probability??

        # ChINaxonDensity = ("6*5000*1e12/3*np.exp(-d/60e-6)",350e-6)

        # func type, density function, max axon radius
        # OLD: ChIN_axon_density = ("r", "5000*1e12/3*np.exp(-r/120e-6)", 350e-6)
        ChIN_axon_density = ("r", "5000*1e12/3*exp(-r/120e-6)", 350e-6)

        self.add_neurons(name="ChIN", neuron_dir=ChIN_dir,
                         num_neurons=self.num_ChIN,
                         axon_density=ChIN_axon_density,
                         volume_id="Striatum")

        ############################################################################

        # Add LTS neuron

        # OBS, the SWC coordinates assume that the soma is centred at 0,0,0
        # Func type, Density function, [[xmin,xmax,ymin,ymax,zmin,zmax]], nAxonPoints

        # See plotLTSdensity.py

        # LTS_density_str = "12*3000*1e12*( 0.25*np.exp(-(((x-200e-6)/100e-6)**2 + ((y-0)/50e-6)**2 + ((z-0)/30e-6)**2)) + 1*np.exp(-(((x-300e-6)/300e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/10e-6)**2)) + 1*np.exp(-(((x-700e-6)/100e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/15e-6)**2)) )",
        LTS_density_str = ("12*3000*1e12*( 0.25*exp(-(((x-200e-6)/100e-6)**2 "
                           "+ ((y-0)/50e-6)**2 + ((z-0)/30e-6)**2)) "
                           "+ 1*exp(-(((x-300e-6)/300e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/10e-6)**2)) "
                           "+ 1*exp(-(((x-700e-6)/100e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/15e-6)**2)) )")

        LTS_axon_density = ("xyz",
                            LTS_density_str,
                            [-200e-6, 900e-6, -100e-6, 100e-6, -30e-6, 30e-6])

        # !!! Remember to update bounding box

        self.add_neurons(name="LTS", neuron_dir=LTS_dir,
                         num_neurons=self.num_LTS,
                         axon_density=LTS_axon_density,
                         volume_id="Striatum")

        # Define FS targets

        # Szydlowski SN, Pollak Dorocic I, Planert H, Carlen M, Meletis K,
        # Silberberg G (2013) Target selectivity of feedforward inhibition
        # by striatal fast-spiking interneurons. J Neurosci
        # --> FS does not target ChIN

        # FS_dist_dep_pruning = "np.exp(-(0.5*d/60e-6)**2)"  # updated 2019-10-31
        FS_dist_dep_pruning = "exp(-(0.5*d/60e-6)**2)"  # Using numexpr.evaluate now, so no np. needed
        # Temp disable dist dep pruning
        # FSDistDepPruning = None
        FS_gGABA = [1.1e-9, 1.5e-9]  # cond (1nS Gittis et al 2010), condStd
        FS_to_LTS_gGABA = [1.1e-10, 1.5e-10]  # cond (1nS Gittis et al 2010), condStd
        FS_gGapJunction = [0.5e-9, 0.1e-9]
        # (gap junctions: 0.5nS, P=0.3 -- Galarreta Hestrin 2002, Koos Tepper 1999)
        # total 8.4nS ?? Gittis et al 2010??

        # File with FS->FS parameters (dont have one yet)
        pfFSFS = None  # Gittis 2010?
        pfFSLTS = None

        # pfFSdSPN = "synapses/v1/trace_table.txt-FD-model-parameters.json"
        # pfFSiSPN = "synapses/v1/trace_table.txt-FI-model-parameters.json"
        pfFSdSPN = os.path.join("$SNUDDA_DATA", "synapses", "striatum", "PlanertFitting-FD-tmgaba-fit.json")
        pfFSiSPN = os.path.join("$SNUDDA_DATA", "synapses", "striatum", "PlanertFitting-FI-tmgaba-fit.json")

        # Increased from a3=0.1 to a3=0.7 to match FS-FS connectivity from Gittis
        self.add_neuron_target(neuron_name="FS",
                               target_name="FS",
                               connection_type="GABA",
                               dist_pruning=None,
                               f1=0.15, soft_max=5, mu2=2, a3=1,
                               conductance=FS_gGABA,
                               cluster_synapses=cluster_FS_synapses,
                               parameter_file=pfFSFS,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.33e-3, 1e3),
                                                         "tau2": (5.7e-3, 1e3)})
        # !!! Double check that channelParamDictionary works, and SI units gets
        # converted to natural units

        self.add_neuron_target(neuron_name="FS",
                               target_name="dSPN",
                               connection_type="GABA",
                               dist_pruning=FS_dist_dep_pruning,
                               f1=0.5, soft_max=5, mu2=2, a3=1.0,
                               conductance=FS_gGABA,
                               cluster_synapses=cluster_FS_synapses,
                               parameter_file=pfFSdSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.2e-3, 1e3),
                                                         "tau2": (8e-3, 1e3)})

        self.add_neuron_target(neuron_name="FS",
                               target_name="iSPN",
                               connection_type="GABA",
                               dist_pruning=FS_dist_dep_pruning,
                               f1=0.5, soft_max=5, mu2=2, a3=0.9,
                               conductance=FS_gGABA,
                               cluster_synapses=cluster_FS_synapses,
                               parameter_file=pfFSiSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.2e-3, 1e3),
                                                         "tau2": (8e-3, 1e3)})

        self.add_neuron_target(neuron_name="FS",
                               target_name="LTS",
                               connection_type="GABA",
                               dist_pruning=None,
                               f1=0.15, soft_max=3, mu2=2, a3=1.0,
                               conductance=FS_to_LTS_gGABA,
                               cluster_synapses=cluster_FS_synapses,
                               parameter_file=pfFSLTS,
                               mod_file="tmGabaA",
                               channel_param_dictionary=None)

        # FS-FS gap junction, currently without pruning
        if True:
            self.add_neuron_target(neuron_name="FS",
                                   target_name="FS",
                                   connection_type="GapJunction",
                                   dist_pruning=None,
                                   f1=0.7, soft_max=8, mu2=2, a3=1.0,
                                   conductance=FS_gGapJunction,
                                   cluster_synapses=False,
                                   channel_param_dictionary=None)

        ## Define MSD1 targets

        # 3e-6 voxel method
        MSP11 = 1.0  # 0.55
        MSP12 = 1.0  # 0.20

        # Taverna 2008, fig 3E&F:
        # D1D1 22.6+/-3pS per synapse, 37+/-15 synapses (approx)
        # D2D1 24.6+/-6pS per synapse, 75+/-30 synapses (approx)
        # D2D2 24+/-1.5pS per synapse, 78+/-11 synapses (approx)

        # !!! But Taverna 2008 analyse aggregates all synapses into a conductance
        # measure?? if so, we need to divide the values by 3 or 4.
        #

        # !!! UPDATE: Assume 24pS per channel, and 10 channels per synapse

        MSD1gGABA = [0.24e-9, 0.1e-9]
        # Koos, Tepper 1999 says max 0.75nS?
        MSD1GABAfailRate = 0.7  # Taverna 2008, figure 2

        # OLD: Previously: 23pA * 50 receptors = 1.15e-9 -- Taverna 2008, fig3
        # OLD: std ~ +/- 8 receptors, we used before:  [1.15e-9, 0.18e-9]

        # !!! TODO: When this runs we do not know how many population units will be added...

        P11withinUnit = MSP11 * within_population_unit_SPN_modifier
        P11betweenUnit = MSP11 * between_population_unit_SPN_modifier
        P12withinUnit = MSP12 * within_population_unit_SPN_modifier
        P12betweenUnit = MSP12 * between_population_unit_SPN_modifier

        # pfdSPNdSPN = "synapses/v1/trace_table.txt-DD-model-parameters.json"
        # pfdSPNiSPN = "synapses/v1/trace_table.txt-DI-model-parameters.json"
        pfdSPNdSPN = os.path.join("$SNUDDA_DATA", "synapses", "striatum", "PlanertFitting-DD-tmgaba-fit.json")
        pfdSPNiSPN = os.path.join("$SNUDDA_DATA", "synapses", "striatum", "PlanertFitting-DI-tmgaba-fit.json")
        pfdSPNChIN = None

        # Argument for distance dependent SPN-SPN synapses:
        # Koos, Tepper, Wilson 2004 -- SPN-SPN more distally

        # From this paper, https://www.frontiersin.org/articles/10.3389/fnana.2010.00150/full,
        #
        # This is in contrast to the axon collateral synapses between SPNs
        # (Tunstall et al., 2002), which typically evoke significantly
        # smaller IPSPs/IPSCs than FSI-evoked synaptic responses when
        # recorded somatically (Koós et al., 2004; Tepper et al., 2004,
        # 2008; Gustafson et al., 2006) due to a combination of
        # predominantly distal synaptic locations (88%; Wilson and Groves,
        # 1980) and relatively few synaptic (2–3) connections made by each
        # SPN on each postsynaptic SPN (Koós et al., 2004)
        #
        # Also, In Kai's Thesis on the first page, He used this reference,
        # https://www.sciencedirect.com/science/article/pii/S0166223612001191?via%3Dihub,
        #

        # With Taverna conductances, we see that the response is much stronger than Planert 2010.
        # We try to introduce distance dependent pruning to see if removing strong proximal synapses
        # will give a better match to experimental data.

        # SPN2SPNdistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)"
        SPN2SPNdistDepPruning = "1-exp(-(0.4*d/60e-6)**2)"

        # Chuhma about 20pA response from 10% SPN, we need to reduce activity, try dist dep pruning
        # (already so few synapses and connectivity)
        # SPN2ChINDistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)"
        SPN2ChINDistDepPruning = "1-exp(-(0.4*d/60e-6)**2)"

        # old f1 = 0.15
        self.add_neuron_target(neuron_name="dSPN",
                               target_name="dSPN",
                               connection_type="GABA",
                               dist_pruning=SPN2SPNdistDepPruning,
                               f1=0.38*0.75, soft_max=3, mu2=2.4,
                               a3=P11withinUnit,
                               a3_other=P11betweenUnit,
                               conductance=MSD1gGABA,
                               cluster_synapses=cluster_SPN_synapses,
                               parameter_file=pfdSPNdSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.3e-3, 1e3),
                                                         "tau2": (12.4e-3, 1e3),
                                                         "failRate": MSD1GABAfailRate})

        # old f1 = 0.15
        self.add_neuron_target(neuron_name="dSPN",
                               target_name="iSPN",
                               connection_type="GABA",
                               dist_pruning=SPN2SPNdistDepPruning,
                               f1=0.20*0.82, soft_max=3, mu2=2.4,
                               a3=P12withinUnit,
                               a3_other=P12betweenUnit,
                               conductance=MSD1gGABA,
                               cluster_synapses=cluster_SPN_synapses,
                               parameter_file=pfdSPNiSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.3e-3, 1e3),
                                                         "tau2": (12.4e-3, 1e3),
                                                         "failRate": MSD1GABAfailRate})

        # Doig, Magill, Apicella, Bolam, Sharott 2014:
        # 5166 +/- 285 GABA synapses on ChIN (antag att 95% av dem är från MS?)
        # 2859 +/- Assymetrical (Glut) synapses on ChIN

        # Set a3 pruning to 0.1, to remove 90% of connected pairs
        # removed softMax = 3 (want to get 5000 MSD1+D2 synapses on ChIN)

        self.add_neuron_target(neuron_name="dSPN",
                               target_name="ChIN",
                               connection_type="GABA",
                               dist_pruning=SPN2ChINDistDepPruning,
                               f1=0.1, soft_max=3, mu2=2.4, a3=0.1,
                               conductance=MSD1gGABA,
                               cluster_synapses=cluster_SPN_synapses,
                               parameter_file=pfdSPNChIN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"failRate": MSD1GABAfailRate})

        ## Define MSD2 targets

        # 3e-6 voxel method
        MSP21 = 1.0  # 0.50
        MSP22 = 1.0  # 0.95

        # OLD: 24pA * 51 receptors = 1.15e-9 -- Taverna 2008, fig3
        # OLD: std ~ +/- 10 receptors [1.24e-9, 0.24e-9]

        # Taverna 2008, fig 3E&F:
        # D1D1 22.6+/-3pS per synapse, 37+/-15 synapses (approx)
        # D2D1 24.6+/-6pS per synapse, 75+/-30 synapses (approx)
        # D2D2 24+/-1.5pS per synapse, 78+/-11 synapses (approx)

        # !!! But Taverna 2008 analyse aggregates all synapses into a conductance
        # measure?? if so, we need to divide the values by 3 or 4.
        #

        # !!! UPDATE: Assume 24pS per channel, and 10 channels per synapse
        # Because in Taverna 2008 iSPN has more receptors in total, we increase
        # softMax from 3 to 4

        MSD2gGABA = [0.24e-9, 0.1e-9]
        MSD2GABAfailRate = 0.4  # Taverna 2008, 2mM

        # Voxel method
        P21withinUnit = MSP21 * within_population_unit_SPN_modifier
        P21betweenUnit = MSP21 * between_population_unit_SPN_modifier
        P22withinUnit = MSP22 * within_population_unit_SPN_modifier
        P22betweenUnit = MSP22 * between_population_unit_SPN_modifier

        pfiSPNdSPN = os.path.join("$SNUDDA_DATA", "synapses", "striatum", "PlanertFitting-ID-tmgaba-fit.json")
        pfiSPNiSPN = os.path.join("$SNUDDA_DATA", "synapses", "striatum", "PlanertFitting-II-tmgaba-fit.json")
        pfiSPNChIN = None

        # GABA decay från Taverna 2008

        # old f1 = 0.15
        self.add_neuron_target(neuron_name="iSPN",
                               target_name="dSPN",
                               connection_type="GABA",
                               dist_pruning=SPN2SPNdistDepPruning,
                               f1=0.3*0.93, soft_max=4, mu2=2.4,
                               a3=P21withinUnit,
                               a3_other=P21betweenUnit,
                               conductance=MSD2gGABA,
                               cluster_synapses=cluster_SPN_synapses,
                               parameter_file=pfiSPNdSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.3e-3, 1e3),
                                                         "tau2": (12.4e-3, 1e3),
                                                         "failRate": MSD2GABAfailRate})

        # old f1 = 0.15
        self.add_neuron_target(neuron_name="iSPN",
                               target_name="iSPN",
                               connection_type="GABA",
                               dist_pruning=SPN2SPNdistDepPruning,
                               f1=0.55, soft_max=4, mu2=2.4,
                               a3=P22withinUnit,
                               a3_other=P22betweenUnit,
                               conductance=MSD2gGABA,
                               cluster_synapses=cluster_SPN_synapses,
                               parameter_file=pfiSPNiSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.3e-3, 1e3),
                                                         "tau2": (12.4e-3, 1e3),
                                                         "failRate": MSD2GABAfailRate})

        # See comment for dSPN to ChIN
        self.add_neuron_target(neuron_name="iSPN",
                               target_name="ChIN",
                               connection_type="GABA",
                               dist_pruning=SPN2ChINDistDepPruning,
                               f1=0.1, soft_max=3, mu2=2.4, a3=0.1,
                               conductance=MSD2gGABA,
                               cluster_synapses=cluster_SPN_synapses,
                               parameter_file=pfiSPNChIN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"failRate": MSD2GABAfailRate})

        ## Define ChIN targets

        # Nelson AB, Hammack N, Yang CF, Shah NM, Seal RP, Kreitzer AC
        # (2014) Striatal choliner- gic interneurons Drive GABA release
        # from dopamine terminals. Neuron
        # Mamaligas, Ford 2016 -- connectivity, 2-5ChIN per MS (in slice)

        ChINgGABA = 1e-9  # If just one value given, then gSTD = 0
        ChINgACh = 1e-9  # FIXME

        # Run 1142 -- No mu2
        # Run 1150 -- Mu2 2.4
        # Run 1153 -- Mu2 D1: 5, D2: 10 (för att testa fler värden)

        # Guzman et al 2003 "Dopaminergic Modulation of Axon Collaterals Interconnecting Spiny Neurons of the Rat Striatum"
        # 325 ChIN inputs per MS (2500 * 0.13)

        # Do ChIN co-release GABA?!! otherwise should be ACh

        pfChINdSPN = None
        pfChINiSPN = None
        pfChINLTS = None

        # !!! SET RELEASE TO GABA FOR NOW

        # ================================================================
        # commenting gabaergic ChIN -> SPN connections Feb. 25th 2020 (RL)

        if False:
            self.add_neuron_target(neuron_name="ChIN",
                                   target_name="dSPN",
                                   connection_type="GABA",
                                   dist_pruning=None,
                                   f1=0.5, soft_max=10, mu2=15, a3=0.1,  # SM 15
                                   conductance=ChINgGABA,
                                   cluster_synapses=False,
                                   parameter_file=pfChINdSPN,
                                   mod_file="tmGabaA",
                                   channel_param_dictionary=None)

            # TEST SETTING THIS TO ACh (SHOULD BE GABA), will this change?
            # !!!

            self.add_neuron_target(neuron_name="ChIN",
                                   target_name="iSPN",
                                   connection_type="GABA",
                                   dist_pruning=None,
                                   f1=0.5, soft_max=10, mu2=10, a3=0.1,  # SM 12
                                   conductance=ChINgGABA,
                                   cluster_synapses=False,
                                   parameter_file=pfChINiSPN,
                                   mod_file="tmGabaA",
                                   channel_param_dictionary=None)
        # ================================================================

        # We got an increasing connection distribution with distance, looks fishy
        # !!! Should be ACh, lets try set it to GABA and see if that changes things
        # --- trying same pruning as for ChIN to MSD2
        if False:
            self.add_neuron_target(neuron_name="ChIN",
                                   target_name="LTS",
                                   connection_type="ACh",
                                   dist_pruning=None,
                                   f1=0.5, soft_max=None, mu2=10, a3=None,  # SM 12
                                   conductance=ChINgACh,
                                   cluster_synapses=False,
                                   parameter_file=pfChINLTS,
                                   mod_file="ACh",  # !!! DOES NOT YET EXIST --- FIXME
                                   channel_param_dictionary=None)

        # !!! USE SAME PARAMS FOR FS AS FOR MS??

        # ??? ChIN does not connect to FS and MS directly ???

        # Add targets for LTS neurons

        LTSgGABA = 1e-9  # !!! FIXME
        # LTSgNO = 1e-9

        # LTSDistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)"  # updated 2019-10-31
        LTSDistDepPruning = "1-exp(-(0.4*d/60e-6)**2)"  # using numexpr.evaluate now, so no np.

        # !!! Straub, Sabatini 2016
        # No LTS synapses within 70 micrometers of proximal MS dendrite
        # !!! ADD DISTANCE DEPENDENT PRUNING

        pfLTSdSPN = None
        pfLTSiSPN = None
        pfLTSChIN = None

        self.add_neuron_target(neuron_name="LTS",
                               target_name="dSPN",
                               connection_type="GABA",
                               dist_pruning=LTSDistDepPruning,
                               f1=1.0*0.3, soft_max=15, mu2=3, a3=0.3,
                               conductance=LTSgGABA,
                               cluster_synapses=False,
                               parameter_file=pfLTSdSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (3e-3, 1e3),
                                                         "tau2": (38e-3, 1e3)})
        # LTS -> SPN, rise time 3+/-0.1 ms, decay time 38+/-3.1 ms, Straub 2016

        self.add_neuron_target(neuron_name="LTS",
                               target_name="iSPN",
                               connection_type="GABA",
                               dist_pruning=LTSDistDepPruning,
                               f1=1.0*0.3, soft_max=15, mu2=3, a3=0.3,
                               conductance=LTSgGABA,
                               cluster_synapses=False,
                               parameter_file=pfLTSiSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (3e-3, 1e3),
                                                         "tau2": (38e-3, 1e3)})

        self.add_neuron_target(neuron_name="LTS",
                               target_name="ChIN",
                               connection_type="GABA",  # also NO, nitric oxide
                               dist_pruning=None,
                               f1=0.5, soft_max=10, mu2=3, a3=0.4,
                               conductance=LTSgGABA,
                               cluster_synapses=False,
                               parameter_file=pfLTSChIN,
                               mod_file="tmGabaA",
                               channel_param_dictionary=None)

    ############################################################################

    def define_GPe(self, num_neurons, d_min=None, neurons_dir=None):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        self.num_gpe_neurons = num_neurons
        self.num_neurons_total += num_neurons

        self.define_structure(struct_name="GPe",
                              struct_mesh="mesh/GPe-mesh.obj",
                              d_min=d_min)

        # !!! Need to add targets for neurons in GPe

    ############################################################################

    def define_GPi(self, num_neurons, d_min=None, neurons_dir=None):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        self.num_gpi_neurons = num_neurons
        self.num_neurons_total += num_neurons

        self.define_structure(struct_name="GPi",
                              struct_mesh="mesh/GPi-mesh.obj",
                              d_min=d_min)

        # !!! Need to add targets for neurons in GPi

    ############################################################################

    def define_STN(self, num_neurons, d_min=None, neurons_dir=None):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        self.num_stn_neurons = num_neurons
        self.num_neurons_total += num_neurons

        self.define_structure(struct_name="STN",
                              struct_mesh="mesh/STN-mesh.obj",
                              d_min=d_min)

        # !!! Need to add targets for neurons in STN

    ############################################################################

    def define_SNr(self, num_neurons, d_min=None, neurons_dir=None):

        if num_neurons <= 0:
            # No neurons, skipping
            return

        self.num_snr_neurons = num_neurons
        self.num_neurons_total += num_neurons

        self.define_structure(struct_name="SNr",
                              struct_mesh="mesh/SNr-mesh.obj",
                              mesh_bin_width=1e-4,
                              d_min=d_min)

        # !!! Need to add targets for neurons in SNr

    ############################################################################

    # TODO: neurons_dir not used here yet
    def define_cortex(self, num_neurons, d_min=None, neurons_dir=None):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        # Neurons with corticostriatal axons
        self.num_cortex_neurons = num_neurons

        self.num_neurons_total += num_neurons

        # Using start location of neuron  DOI: 10.25378/janelia.5521780 for centre
        # !!! If we use a larger mesh for cortex, we will need to reduce
        #     meshBinWidth to 1e-4 (or risk getting memory error)
        self.define_structure(struct_name="Cortex",
                              struct_mesh="cube",
                              struct_centre=np.array([7067e-6, 3007e-6, 2570e-6]),
                              side_len=200e-6,
                              mesh_bin_width=5e-5,
                              d_min=d_min)

        cortex_dir = os.path.join("$SNUDDA_DATA", "InputAxons", "Cortex", "Reg10")

        # Add cortex axon

        self.add_neurons("CortexAxon", cortex_dir, self.num_cortex_neurons,
                         model_type="virtual",
                         rotation_mode="",
                         volume_id="Cortex")

        # Define targets

        cortex_glut_cond = [1e-9, 0.1e-9]

        # We should have both ipsi and contra, M1 and S1 input, for now
        # picking one
        cortexSynParMS = os.path.join("$SNUDDA_DATA", "synapses", "striatum", "M1RH_Analysis_190925.h5-parameters-MS.json")
        cortexSynParFS = os.path.join("$SNUDDA_DATA", "synapses", "striatum", "M1RH_Analysis_190925.h5-parameters-FS.json")

        self.add_neuron_target(neuron_name="CortexAxon",
                               target_name="dSPN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               parameter_file=cortexSynParMS,
                               mod_file="tmGlut",
                               conductance=cortex_glut_cond,
                               cluster_synapses=False,
                               channel_param_dictionary=None)

        self.add_neuron_target(neuron_name="CortexAxon",
                               target_name="iSPN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               parameter_file=cortexSynParMS,
                               mod_file="tmGlut",
                               conductance=cortex_glut_cond,
                               cluster_synapses=False,
                               channel_param_dictionary=None)

        self.add_neuron_target(neuron_name="CortexAxon",
                               target_name="FS",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               parameter_file=cortexSynParFS,
                               mod_file="tmGlut",
                               conductance=cortex_glut_cond,
                               cluster_synapses=False,
                               channel_param_dictionary=None)

        # !!! No input for LTS and ChIN right now...

    ############################################################################

    # TODO: neurons_dir not used here yet
    def define_thalamus(self, num_neurons, d_min=None, neurons_dir=None):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        # Neurons with thalamustriatal axons
        self.num_thalamus_neurons = num_neurons

        self.num_neurons_total += num_neurons

        # Using start location of neuron DOI: 10.25378/janelia.5521765 for centre
        self.define_structure(struct_name="Thalamus",
                              struct_mesh="cube",
                              struct_centre=np.array([4997e-6, 4260e-6, 7019e-6]),
                              side_len=200e-6,
                              mesh_bin_width=5e-5,
                              d_min=d_min)

        # Define neurons

        thalamus_dir = os.path.join("$SNUDDA_DATA", "morphology", "InputAxons", "Thalamus", "Reg10")

        self.add_neurons("ThalamusAxon", thalamus_dir, self.num_thalamus_neurons,
                         model_type="virtual",
                         rotation_mode="",
                         volume_id="Thalamus")

        # Define targets

        thalamus_syn_par_ms = os.path.join("$SNUDDA_DATA", "synapses", "striatum",
                                        "TH_Analysis_191001.h5-parameters-MS.json")
        thalamus_syn_par_fs = os.path.join("$SNUDDA_DATA", "synapses", "striatum",
                                        "TH_Analysis_191001.h5-parameters-FS.json")

        thalamus_glut_cond = [1e-9, 0.1e-9]

        self.add_neuron_target(neuron_name="ThalamusAxon",
                               target_name="dSPN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               conductance=thalamus_glut_cond,
                               cluster_synapses=False,
                               parameter_file=thalamus_syn_par_ms,
                               mod_file="tmGlut",
                               channel_param_dictionary=None)

        self.add_neuron_target(neuron_name="ThalamusAxon",
                               target_name="iSPN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               conductance=thalamus_glut_cond,
                               cluster_synapses=False,
                               parameter_file=thalamus_syn_par_ms,
                               mod_file="tmGlut",
                               channel_param_dictionary=None)

        # Picked D1 parameters, lack
        self.add_neuron_target(neuron_name="ThalamusAxon",
                               target_name="FS",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               conductance=thalamus_glut_cond,
                               cluster_synapses=False,
                               parameter_file=thalamus_syn_par_fs,
                               mod_file="tmGlut",
                               channel_param_dictionary=None)

    ############################################################################

    @staticmethod
    def setup_random_seeds(random_seed=None):

        seed_types = ["init", "place", "detect", "project", "prune", "input", "simulate"]
        # print(f"Seeding with rand_seed={random_seed}")

        ss = np.random.SeedSequence(random_seed)
        all_seeds = ss.generate_state(len(seed_types))

        rand_seed_dict = collections.OrderedDict()

        if random_seed is not None:
            rand_seed_dict["masterseed"] = random_seed

        for st, s in zip(seed_types, all_seeds):
            rand_seed_dict[st] = s
            # print(f"Random seed {st} to {s}")

        init_rng = np.random.default_rng(rand_seed_dict["init"])

        return rand_seed_dict, init_rng


if __name__ == "__main__":

    full_striatum = True  # False #True

    # Striatum has about 1.73 million neurons in mouse

    # Rat data (Oorschot 1996 J Comp Neurol 366)
    # --> mouse estimated from scaling down to mouse from rat
    # Striatum: 2.79M --> 1.73M
    # GPe: 46,000 --> 28500
    # GPi: 3,200 --> 2000
    # STN: 13,600 --> 8400
    # SNRpc : 7,200 --> 4500
    # SNRpr : 26,300 --> 16300

    # --> SNr = 20800

    # Nd1=Nd2=9493
    # Nfsi=400
    # Nstn=97
    # Nta=82
    # Nti=247
    # Nsnr=189

    if full_striatum:
        struct_def = {"Striatum": 1730000,
                     "GPe": 28500,
                     "GPi": 2000,
                     "SNr": 20800,
                     "STN": 8400,
                     "Cortex": 1,
                     "Thalamus": 1}

        # !!! TEMP, only do stratium for now
        struct_def = {"Striatum": 1730000,
                     "GPe": 0,
                     "GPi": 0,
                     "SNr": 0,
                     "STN": 0,
                     "Cortex": 0,
                     "Thalamus": 0}


    else:
        struct_def = {"Striatum": 100000,
                     "GPe": 0,
                     "GPi": 0,
                     "SNr": 0,
                     "STN": 0,
                     "Cortex": 0,
                     "Thalamus": 0}

    nTotals = 0
    for x in struct_def:
        nTotals += struct_def[x]

    fName = os.path.join("config", f"basal-ganglia-config-{nTotals}.json")

    SnuddaInit(struct_def=struct_def, config_file=fName)
