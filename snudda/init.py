# Rewriting the create network config file to make it more general

# !!! Currently writing full path for the files in the snudda data directory
#     this makes moving between computers difficult. Fix...

#
# Add a function so that $SNUDDADATA refers to the base datapath for snudda
#

import numpy as np
import os.path
import glob
import collections
from .create_cube_mesh import create_cube_mesh
from .create_slice_mesh import create_slice_mesh

import json


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            # return super(NumpyEncoder, self).default(obj)
            return json.JSONEncoder.default(self, obj)


class SnuddaInit(object):

    def __init__(self, struct_def, config_file,
                 num_population_units=1, population_unit_centres="[[]]", population_unit_radius=None,
                 random_seed=None):

        print("CreateConfig")

        self.network_data = collections.OrderedDict([])

        self.network_data["RandomSeed"], self.init_rng = SnuddaInit.setup_random_seeds(random_seed)
        self.network_data["Volume"] = collections.OrderedDict([])
        self.num_neurons_total = 0
        self.config_file = config_file

        if config_file is not None:
            self.basePath = os.path.dirname(config_file)
        else:
            self.basePath = ""

        self.data_path = os.path.join(os.path.dirname(__file__), "data")

        # Population Units here refer to processing units, where the neurons within a Population Unit
        # might have different connectivity than neurons belonging to different population Units
        self.network_data["PopulationUnits"] = collections.OrderedDict([])
        self.network_data["PopulationUnits"]["nPopulationUnits"] = num_population_units

        use_random_population_units = False

        if use_random_population_units:
            self.network_data["PopulationUnits"]["method"] = "random"
        else:
            self.network_data["PopulationUnits"]["method"] = "populationUnitSpheres"

            # Centre of striatum mesh is [3540e-6,4645e-6,5081e-6]
            # - the population units will be shifted according to this coordinate

            if type(population_unit_centres) == str:
                population_unit_centres = eval(population_unit_centres)

            self.network_data["PopulationUnits"]["centres"] = population_unit_centres
            assert len(self.network_data["PopulationUnits"]["centres"]) == num_population_units, \
                (f"The number of Population Units ({num_population_units}) "
                 f"does not equal the number of centres ({len(self.network_data['PopulationUnits']['centres'])})")

            if population_unit_radius:
                self.network_data["PopulationUnits"]["radius"] = population_unit_radius * 1e-6
            else:
                self.network_data["PopulationUnits"]["radius"] = None

            print("Overriding the number of population units")

            self.network_data["PopulationUnits"]["nPopulationUnits"] \
                = len(self.network_data["PopulationUnits"]["centres"])

        self.network_data["Connectivity"] = dict([])
        self.network_data["Neurons"] = dict([])

        # self.neuronTargets = collections.OrderedDict([])

        print("Using " + str(num_population_units) + " Population Units")

        if num_population_units > 1:
            from scipy import spatial

            distance_between_population_units = spatial.distance.cdist(self.network_data["PopulationUnits"]["centres"],
                                                                       self.network_data["PopulationUnits"]["centres"],
                                                                       metric="euclidean")[
                np.triu_indices(len(self.network_data["PopulationUnits"]["centres"]), k=1)]

            print("Using radius of Population Unit Sphere " + str(
                np.ceil(self.network_data["PopulationUnits"]["radius"] * 1e6)) + " microns")

            print("Using distance between Population Unit Centres" + str(
                distance_between_population_units * 1e6) + " microns")

        self.n_population_units = num_population_units  # 5

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
                struct_func[sn](num_neurons=struct_def[sn])

            # Only write JSON file if the structDef was not empty
            self.write_json(self.config_file)
        else:
            print("No structDef defined, not writing JSON file in init")

    ############################################################################

    # meshBinWidth is used when voxelising the mesh to determine which part of
    # space is inside the mesh. For smaller structures we might need to use
    # a smaller meshBinWidth than the default 1e-4

    def define_structure(self,
                         struct_name,
                         struct_mesh,
                         d_min=15e-6,
                         struct_centre=None,
                         side_len=None,
                         slice_depth=None,
                         mesh_bin_width=None):

        if struct_mesh == "cube":
            assert slice_depth is None, "define_structure: sliceDepth is not used for cubes, please set to None"
            assert side_len is not None, "define_structure: cube needs sideLen specified"
            assert struct_centre is not None, "define_structure: cube needs a structCentre"

            struct_mesh = os.path.join(self.basePath, "mesh", f"{struct_name}-cube-mesh-{side_len}.obj")

            if mesh_bin_width is None:
                mesh_bin_width = side_len / 3.0
                print("Setting mesh_bin_width to " + str(mesh_bin_width))

            create_cube_mesh(file_name=struct_mesh,
                             centre_point=struct_centre,
                             side_len=side_len,
                             description=f"{struct_name} cube mesh, centre: {struct_centre}, side: {side_len}")

        elif struct_mesh == "slice":

            struct_mesh = os.path.join(self.basePath, "mesh", f"{struct_name}-slice-mesh-150mum-depth.obj")

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

            # Different pruning rules for within and between neuron units
            con_info["pruningOther"] = pruning_info_other

        # Json did not like tuples in keys, so we separate by comma
        nt_key = neuron_name + "," + target_name
        if nt_key not in self.network_data["Connectivity"]:
            self.network_data["Connectivity"][nt_key] = dict([])

        self.network_data["Connectivity"][nt_key][connection_type] = con_info

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
                    eval(axon_density[1])
                except:
                    print("!!! Axon density failed test: " + str(axon_density))
                    print("Inparameter: r = 1-D array of radius in meter")
            elif axon_density[0] == "xyz":
                x = np.linspace(axon_density[2][0], axon_density[2][1], 10)  # x,y,z used in eval below
                y = np.linspace(axon_density[2][2], axon_density[2][3], 10)
                z = np.linspace(axon_density[2][4], axon_density[2][5], 10)
                try:
                    eval(axon_density[1])
                except:
                    print("!!! Axon density failed test: " + str(axon_density))
                    print("Inparameters: x,y,z three 1-D arrays (units in meter)")
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    exit(-1)

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

                p_corner = eval(axon_density[1]) * (3e-6 ** 3)

                for P, xx, yy, zz in zip(p_corner, x, y, z):
                    print(name + " axon density P(" + str(xx) + "," + str(yy) + "," + str(zz) + ") = " + str(P))

                if (p_corner > 0.01).any():
                    print("Axon density too high at boundary!!")
                    print("Please increase bounding box")
                    exit(-1)

                # print(str(axonDensity[3]) + " " + str(name) \
                #      + " axon points to place")

        print("Adding neurons: " + str(name) + " from dir " + str(neuron_dir))
        # TODO: We should force users to use same name as the directory name
        # ie, fs/FSN_0 directory should be named FSN_0

        # Find which neurons are available in neuronDir
        dir_list = glob.glob(neuron_dir + "/*")
        neuron_file_list = []

        assert len(dir_list) > 0, "Neuron dir " + str(neuron_dir) + " is empty!"

        for d in dir_list:

            if os.path.isdir(d):
                par_file = os.path.join(d, "parameters.json")
                mech_file = os.path.join(d, "mechanisms.json")
                modulation_file = os.path.join(d, "modulation.json")
                if not os.path.exists(modulation_file):
                    modulation_file = None

                swc_file = glob.glob(os.path.join(d, "*swc"))
                hoc_file = glob.glob(os.path.join(d, "*hoc"))

                assert len(swc_file) == 1, "Morph dir " + d + " should contain one swc file"

                assert len(hoc_file) <= 1, "Morph dir " + d + " contains more than one hoc file"

                if len(hoc_file) == 0:
                    hoc_file = [None]

                neuron_file_list.append((d,
                                        swc_file[0],
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

            if not os.path.isfile(par_file) and model_type is not "virtual":
                print(f"Parameter file not found: {par_file}")

            if not os.path.isfile(mech_file) and model_type is not "virtual":
                print(f"Mechanism file not found: {mech_file}")

            if hoc_file is not None and not os.path.isfile(hoc_file):
                print(f"Hoc file not found: {hoc_file}")

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

    ############################################################################

    def write_json(self, filename):

        # Create directory if it does not already exist
        dir_name = os.path.dirname(filename)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        print(f"Writing {filename}")

        import json
        with open(filename, 'w') as f:
            json.dump(self.network_data, f, indent=4, cls=NumpyEncoder)

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

    def define_striatum(self, num_neurons=None,
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
                        # slice_depth=None,
                        cell_spec_dir=None,
                        neuron_density=80500):

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

        if volume_type == "mouseStriatum":
            self.define_structure(struct_name="Striatum",
                                  struct_mesh=os.path.join(self.data_path, "mesh", "Striatum-mesh.obj"),
                                  mesh_bin_width=1e-4)

        elif volume_type == "slice":
            self.define_structure(struct_name="Striatum",
                                  struct_mesh="slice",
                                  side_len=side_len)

        elif num_neurons <= 1e6:  # 1e6
            print("Using cube for striatum")
            # 1.73 million neurons, volume of allen striatal mesh is 21.5mm3
            striatum_volume = 1e-9 * num_neurons / neuron_density  # 80.5e3
            striatum_side_len = striatum_volume ** (1. / 3)
            striatum_centre = np.array([3540e-6, 4645e-6, 5081e-6])

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
                                  mesh_bin_width=mesh_bin_width)

        else:
            # Default, full size striatum
            self.define_structure(struct_name="Striatum",
                                  struct_mesh=os.path.join(self.data_path, "mesh", "Striatum-mesh.obj"),
                                  mesh_bin_width=1e-4)

        if cell_spec_dir is None:
            cs_dir = os.path.join(self.data_path, "cellspecs")
        else:
            cs_dir = cell_spec_dir

        FS_dir = os.path.join(cs_dir, "fs")
        dSPN_dir = os.path.join(cs_dir, "dspn")
        iSPN_dir = os.path.join(cs_dir, "ispn")
        ChIN_dir = os.path.join(cs_dir, "chin")
        LTS_dir = os.path.join(cs_dir, "lts")

        self.reg_size = 5

        if self.n_population_units == 1:
            self.population_unit_SPN_modifier = 1
        else:
            print("!!! OBS, modifying probaiblities within and between channe")
            self.population_unit_SPN_modifier = 0.2  # 0.2 = 20% within, 2 = 2x higher within
            print("populationUnitMSNmodifier: " + str(self.population_unit_SPN_modifier))

        # Add the neurons

        self.add_neurons(name="FSN", neuron_dir=FS_dir,
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
        ChIN_axon_density = ("r", "5000*1e12/3*np.exp(-r/120e-6)", 350e-6)

        self.add_neurons(name="ChIN", neuron_dir=ChIN_dir,
                         num_neurons=self.num_ChIN,
                         axon_density=ChIN_axon_density,
                         volume_id="Striatum")

        ############################################################################

        # Add LTS neuron

        # OBS, the SWC coordinates assume that the soma is centred at 0,0,0
        # Func type, Density function, [[xmin,xmax,ymin,ymax,zmin,zmax]], nAxonPoints

        # See plotLTSdensity.py
        LTS_axon_density = ("xyz",
                            "12*3000*1e12*( 0.25*np.exp(-(((x-200e-6)/100e-6)**2 + ((y-0)/50e-6)**2 + ((z-0)/30e-6)**2)) + 1*np.exp(-(((x-300e-6)/300e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/10e-6)**2)) + 1*np.exp(-(((x-700e-6)/100e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/15e-6)**2)) )",
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

        # FSDistDepPruning = "np.exp(-(0.3*d/60e-6)**2)"
        FS_dist_dep_pruning = "np.exp(-(0.5*d/60e-6)**2)"  # updated 2019-10-31
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
        pfFSdSPN = os.path.join(self.data_path, "synapses", "v2", "PlanertFitting-FD-tmgaba-fit.json")
        pfFSiSPN = os.path.join(self.data_path, "synapses", "v2", "PlanertFitting-FI-tmgaba-fit.json")

        # Increased from a3=0.1 to a3=0.7 to match FS-FS connectivity from Gittis
        self.add_neuron_target(neuron_name="FSN",
                               target_name="FSN",
                               connection_type="GABA",
                               dist_pruning=None,
                               f1=0.15, soft_max=5, mu2=2, a3=1,
                               conductance=FS_gGABA,
                               parameter_file=pfFSFS,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.33e-3, 1e3),
                                                         "tau2": (5.7e-3, 1e3)})
        # !!! Double check that channelParamDictionary works, and SI units gets
        # converted to natural units

        self.add_neuron_target(neuron_name="FSN",
                               target_name="dSPN",
                               connection_type="GABA",
                               dist_pruning=FS_dist_dep_pruning,
                               f1=0.5, soft_max=5, mu2=2, a3=1.0,
                               conductance=FS_gGABA,
                               parameter_file=pfFSdSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.2e-3, 1e3),
                                                         "tau2": (8e-3, 1e3)})

        self.add_neuron_target(neuron_name="FSN",
                               target_name="iSPN",
                               connection_type="GABA",
                               dist_pruning=FS_dist_dep_pruning,
                               f1=0.5, soft_max=5, mu2=2, a3=0.9,
                               conductance=FS_gGABA,
                               parameter_file=pfFSiSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (1.2e-3, 1e3),
                                                         "tau2": (8e-3, 1e3)})

        self.add_neuron_target(neuron_name="FSN",
                               target_name="LTS",
                               connection_type="GABA",
                               dist_pruning=None,
                               f1=0.15, soft_max=3, mu2=2, a3=1.0,
                               conductance=FS_to_LTS_gGABA,
                               parameter_file=pfFSLTS,
                               mod_file="tmGabaA",
                               channel_param_dictionary=None)

        # FS-FS gap junction, currently without pruning
        if True:
            self.add_neuron_target(neuron_name="FSN",
                                   target_name="FSN",
                                   connection_type="GapJunction",
                                   dist_pruning=None,
                                   f1=0.7, soft_max=8, mu2=2, a3=1.0,
                                   conductance=FS_gGapJunction,
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

        P11withinUnit = MSP11 * self.population_unit_SPN_modifier
        P11betweenUnit = MSP11 * (1 + (1 - self.population_unit_SPN_modifier) / self.n_population_units)
        P12withinUnit = MSP12 * self.population_unit_SPN_modifier
        P12betweenUnit = MSP12 * (1 + (1 - self.population_unit_SPN_modifier) / self.n_population_units)

        # pfdSPNdSPN = "synapses/v1/trace_table.txt-DD-model-parameters.json"
        # pfdSPNiSPN = "synapses/v1/trace_table.txt-DI-model-parameters.json"
        pfdSPNdSPN = os.path.join(self.data_path, "synapses", "v2", "PlanertFitting-DD-tmgaba-fit.json")
        pfdSPNiSPN = os.path.join(self.data_path, "synapses", "v2", "PlanertFitting-DI-tmgaba-fit.json")
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
        SPN2SPNdistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)"

        # Chuhma about 20pA response from 10% SPN, we need to reduce activity, try dist dep pruning
        # (already so few synapses and connectivity)
        SPN2ChINDistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)"

        # old f1 = 0.15
        self.add_neuron_target(neuron_name="dSPN",
                               target_name="dSPN",
                               connection_type="GABA",
                               dist_pruning=SPN2SPNdistDepPruning,
                               f1=0.38, soft_max=3, mu2=2.4,
                               a3=P11withinUnit,
                               a3_other=P11betweenUnit,
                               conductance=MSD1gGABA,
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
                               f1=0.20, soft_max=3, mu2=2.4,
                               a3=P12withinUnit,
                               a3_other=P12betweenUnit,
                               conductance=MSD1gGABA,
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
        P21withinUnit = MSP21 * self.population_unit_SPN_modifier
        P21betweenUnit = MSP21 * (1 + (1 - self.population_unit_SPN_modifier) / self.n_population_units)
        P22withinUnit = MSP22 * self.population_unit_SPN_modifier
        P22betweenUnit = MSP22 * (1 + (1 - self.population_unit_SPN_modifier) / self.n_population_units)

        pfiSPNdSPN = os.path.join(self.data_path, "synapses", "v2", "PlanertFitting-ID-tmgaba-fit.json")
        pfiSPNiSPN = os.path.join(self.data_path, "synapses", "v2", "PlanertFitting-II-tmgaba-fit.json")
        pfiSPNChIN = None

        # GABA decay från Taverna 2008

        # old f1 = 0.15
        self.add_neuron_target(neuron_name="iSPN",
                               target_name="dSPN",
                               connection_type="GABA",
                               dist_pruning=SPN2SPNdistDepPruning,
                               f1=0.3, soft_max=4, mu2=2.4,
                               a3=P21withinUnit,
                               a3_other=P21betweenUnit,
                               conductance=MSD2gGABA,
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
                                   parameter_file=pfChINLTS,
                                   mod_file="ACh",  # !!! DOES NOT YET EXIST --- FIXME
                                   channel_param_dictionary=None)

        # !!! USE SAME PARAMS FOR FS AS FOR MS??

        # ??? ChIN does not connect to FS and MS directly ???

        # Add targets for LTS neurons

        LTSgGABA = 1e-9  # !!! FIXME
        # LTSgNO = 1e-9

        LTSDistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)"  # updated 2019-10-31

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
                               f1=1.0, soft_max=15, mu2=3, a3=0.3,
                               conductance=LTSgGABA,
                               parameter_file=pfLTSdSPN,
                               mod_file="tmGabaA",
                               channel_param_dictionary={"tau1": (3e-3, 1e3),
                                                         "tau2": (38e-3, 1e3)})
        # LTS -> SPN, rise time 3+/-0.1 ms, decay time 38+/-3.1 ms, Straub 2016

        self.add_neuron_target(neuron_name="LTS",
                               target_name="iSPN",
                               connection_type="GABA",
                               dist_pruning=LTSDistDepPruning,
                               f1=1.0, soft_max=15, mu2=3, a3=0.3,
                               conductance=LTSgGABA,
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
                               parameter_file=pfLTSChIN,
                               mod_file="tmGabaA",
                               channel_param_dictionary=None)

    ############################################################################

    def define_GPe(self, num_neurons):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        self.num_GPe_neurons = num_neurons
        self.num_neurons_total += num_neurons

        self.define_structure(struct_name="GPe",
                              struct_mesh="mesh/GPe-mesh.obj")

        # !!! Need to add targets for neurons in GPe

    ############################################################################

    def define_GPi(self, num_neurons):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        self.num_GPi_neurons = num_neurons
        self.num_neurons_total += num_neurons

        self.define_structure(struct_name="GPi",
                              struct_mesh="mesh/GPi-mesh.obj")

        # !!! Need to add targets for neurons in GPi

    ############################################################################

    def define_STN(self, num_neurons):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        self.num_STN_neurons = num_neurons
        self.num_neurons_total += num_neurons

        self.define_structure(struct_name="STN",
                              struct_mesh="mesh/STN-mesh.obj")

        # !!! Need to add targets for neurons in STN

    ############################################################################

    def define_SNr(self, num_neurons):

        if num_neurons <= 0:
            # No neurons, skipping
            return

        self.num_SNr_neurons = num_neurons
        self.num_neurons_total += num_neurons

        self.define_structure(struct_name="SNr",
                              struct_mesh="mesh/SNr-mesh.obj",
                              mesh_bin_width=1e-4)

        # !!! Need to add targets for neurons in SNr

    ############################################################################

    def define_cortex(self, num_neurons):

        if num_neurons <= 0:
            # No neurons specified, skipping structure
            return

        # Neurons with corticostriatal axons
        self.num_Cortex_neurons = num_neurons

        self.num_neurons_total += num_neurons

        # Using start location of neuron  DOI: 10.25378/janelia.5521780 for centre
        # !!! If we use a larger mesh for cortex, we will need to reduce
        #     meshBinWidth to 1e-4 (or risk getting memory error)
        self.define_structure(struct_name="Cortex",
                              struct_mesh="cube",
                              struct_centre=np.array([7067e-6, 3007e-6, 2570e-6]),
                              side_len=200e-6,
                              mesh_bin_width=5e-5)

        cortex_dir = "morphology/InputAxons/Cortex/Reg10/"

        # Add cortex axon

        self.add_neurons("CortexAxon", cortex_dir, self.num_Cortex_neurons,
                         model_type="virtual",
                         rotation_mode="",
                         volume_id="Cortex")

        # Define targets

        cortex_glut_cond = [1e-9, 0.1e-9]

        # We should have both ipsi and contra, M1 and S1 input, for now
        # picking one
        cortexSynParMS = os.path.join(self.data_path, "synapses", "v2", "M1RH_Analysis_190925.h5-parameters-MS.json")
        cortexSynParFS = os.path.join(self.data_path, "synapses", "v2", "M1RH_Analysis_190925.h5-parameters-FS.json")

        self.add_neuron_target(neuron_name="CortexAxon",
                               target_name="dSPN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               parameter_file=cortexSynParMS,
                               mod_file="tmGlut",
                               conductance=cortex_glut_cond,
                               channel_param_dictionary=None)

        self.add_neuron_target(neuron_name="CortexAxon",
                               target_name="iSPN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               parameter_file=cortexSynParMS,
                               mod_file="tmGlut",
                               conductance=cortex_glut_cond,
                               channel_param_dictionary=None)

        self.add_neuron_target(neuron_name="CortexAxon",
                               target_name="FSN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               parameter_file=cortexSynParFS,
                               mod_file="tmGlut",
                               conductance=cortex_glut_cond,
                               channel_param_dictionary=None)

        # !!! No input for LTS and ChIN right now...

    ############################################################################

    def define_thalamus(self, num_neurons):

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
                              mesh_bin_width=5e-5)

        # Define neurons

        thalamus_dir = "morphology/InputAxons/Thalamus/Reg10/"

        self.add_neurons("ThalamusAxon", thalamus_dir, self.num_thalamus_neurons,
                         model_type="virtual",
                         rotation_mode="",
                         volume_id="Thalamus")

        # Define targets

        thalamusSynParMS = os.path.join(self.data_path, "synapses", "v2", "TH_Analysis_191001.h5-parameters-MS.json")
        thalamusSynParFS = os.path.join(self.data_path, "synapses", "v2", "TH_Analysis_191001.h5-parameters-FS.json")

        ThalamusGlutCond = [1e-9, 0.1e-9]

        self.add_neuron_target(neuron_name="ThalamusAxon",
                               target_name="dSPN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               conductance=ThalamusGlutCond,
                               parameter_file=thalamusSynParMS,
                               mod_file="tmGlut",
                               channel_param_dictionary=None)

        self.add_neuron_target(neuron_name="ThalamusAxon",
                               target_name="iSPN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               conductance=ThalamusGlutCond,
                               parameter_file=thalamusSynParMS,
                               mod_file="tmGlut",
                               channel_param_dictionary=None)

        # Picked D1 parameters, lack
        self.add_neuron_target(neuron_name="ThalamusAxon",
                               target_name="FSN",
                               connection_type="AMPA_NMDA",
                               dist_pruning=None,
                               f1=1.0, soft_max=3, mu2=2.4, a3=None,
                               conductance=ThalamusGlutCond,
                               parameter_file=thalamusSynParFS,
                               mod_file="tmGlut",
                               channel_param_dictionary=None)

    ############################################################################

    @staticmethod
    def setup_random_seeds(random_seed=None):

        seed_types = ["init", "place", "detect", "prune", "input", "simulate"]
        print(f"Seeding with rand_seed={random_seed}")

        ss = np.random.SeedSequence(random_seed)
        all_seeds = ss.generate_state(len(seed_types))

        rand_seed_dict = collections.OrderedDict()

        if random_seed is not None:
            rand_seed_dict["masterseed"] = random_seed

        for st, s in zip(seed_types, all_seeds):
            rand_seed_dict[st] = s
            print(f"Random seed {st} to {s}")

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
        structDef = {"Striatum": 1730000,
                     "GPe": 28500,
                     "GPi": 2000,
                     "SNr": 20800,
                     "STN": 8400,
                     "Cortex": 1,
                     "Thalamus": 1}

        # !!! TEMP, only do stratium for now
        structDef = {"Striatum": 1730000,
                     "GPe": 0,
                     "GPi": 0,
                     "SNr": 0,
                     "STN": 0,
                     "Cortex": 0,
                     "Thalamus": 0}


    else:
        structDef = {"Striatum": 100000,
                     "GPe": 0,
                     "GPi": 0,
                     "SNr": 0,
                     "STN": 0,
                     "Cortex": 0,
                     "Thalamus": 0}

    nTotals = 0
    for x in structDef:
        nTotals += structDef[x]

    fName = os.path.join("config", f"basal-ganglia-config-{nTotals}.json")

    SnuddaInit(struct_def=structDef, config_file=fName, num_population_units=1)
