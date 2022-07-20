from snudda.init import SnuddaInit
import numpy as np
import os


class neuromodulationInit(SnuddaInit):

    def __init__(self, network_path=None,
                 struct_def=None, config_file=None,
                 random_seed=None):

        super(neuromodulationInit, self).__init__(network_path=network_path,
                                                  struct_def=struct_def, config_file=config_file,
                                                  random_seed=random_seed)

        self.num_LTS = None
        self.num_ChIN = None
        self.num_iSPN = None
        self.num_dSPN = None
        self.num_FS = None
        struct_func = {"Striatum": self.define_striatum,
                       "GPe": self.define_GPe,
                       "GPi": self.define_GPi,
                       "STN": self.define_STN,
                       "SNr": self.define_SNr,
                       "Cortex": self.define_cortex,
                       "Thalamus": self.define_thalamus,
                       "SNc": self.define_snc}

        if struct_def:
            for sn in struct_def:
                print("Adding " + sn + " with " + str(struct_def[sn]) + " neurons")
                struct_func[sn](num_neurons=struct_def[sn])

            # Only write JSON file if the structDef was not empty
            self.write_json(self.config_file)
        else:
            pass
            # print("No structDef defined, not writing JSON file in init")

    def define_snc(self, nNeurons, neuron_dir):

        if (nNeurons <= 0):
            # No neurons specified, skipping structure
            return

        # Neurons with corticostriatal axons
        self.nSNc = nNeurons


        striatum_volume = 1e-9 * 10 / 8e4  # 80.5e3
        striatum_side_len = striatum_volume ** (1. / 3)

        mesh_bin_width = striatum_side_len

        self.define_structure(struct_name="SNc",
                              struct_mesh="cube",
                              struct_centre=np.array([3540e-6, 4645e-6, 5081e-6]),
                              mesh_bin_width=mesh_bin_width,
                              side_len=striatum_side_len)

        SNcDir = os.path.abspath(os.path.join(neuron_dir, 'SNc'))

        SNcaxonDensity = ("r", "50000*1e12/3*exp(-r/120e-6)", 500e-6)

        # Add Dopaminergic Axon

        self.add_neurons(name="DopaminergicAxon", neuron_dir=SNcDir, num_neurons=self.nSNc, \
                         axon_density=SNcaxonDensity,
                         volume_id="SNc")

        # Define targets

        DopaminergicCond = [1e-9, 0.1e-9]

        self.add_neuron_target(neuron_name="DopaminergicAxon",
                               target_name="dSPN",
                               connection_type="Dopamine",
                               dist_pruning=None,
                               f1=None, soft_max=None, mu2=None, a3=None,
                               conductance=DopaminergicCond,
                               parameter_file=None,
                               mod_file="concDA",
                               channel_param_dictionary=None)

        self.add_neuron_target(neuron_name="DopaminergicAxon",
                               target_name="iSPN",
                               connection_type="Dopamine",
                               dist_pruning=None,
                               f1=None, soft_max=None, mu2=None, a3=None,
                               conductance=DopaminergicCond,
                               parameter_file=None,
                               mod_file="concDA",
                               channel_param_dictionary=None)

    def define_striatum_neuromodulation(self,
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
                                        neurons_dir=None,
                                        neuron_density=80500,
                                        population_unit_SPN_modifier=1):

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
                print("Striatum should have " + str(num_neurons) + " but " + str(self.num_neurons_total)
                      + " are being requested, check fractions set for defineStriatum.")

        if volume_type == "mouseStriatum":
            self.define_structure(struct_name="Striatum",
                                  struct_mesh=os.path.join("$DATA", "mesh", "Striatum-d.obj"),
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
                                  struct_mesh=os.path.join("$DATA", "mesh", "Striatum-d.obj"),
                                  mesh_bin_width=1e-4)

        if neurons_dir is None:
            neurons_dir = os.path.join("$DATA", "neurons")

        FS_dir = os.path.join(neurons_dir, "striatum", "fs")
        dSPN_dir = os.path.join(neurons_dir, "striatum", "dspn")
        iSPN_dir = os.path.join(neurons_dir, "striatum", "ispn")
        ChIN_dir = os.path.join(neurons_dir, "striatum", "chin")
        LTS_dir = os.path.join(neurons_dir, "striatum", "lts")

        self.reg_size = 5

        if "PopulationUnits" in self.network_data and "Striatum" in self.network_data["PopulationUnits"]:
            num_population_units = len(self.network_data["PopulationUnits"]["Striatum"]["unitID"])
        else:
            num_population_units = 1

        if num_population_units <= 1:
            population_unit_SPN_modifier = 1

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

        ChIN_axon_density = ("r", "50000*1e12/3*exp(-r/120e-6)", 500e-6)

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
        pfFSdSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-FD-tmgaba-fit.json")
        pfFSiSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-FI-tmgaba-fit.json")

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

        # !!! TODO: When this runs we do not know how many population units will be added...

        P11withinUnit = MSP11 * population_unit_SPN_modifier
        P11betweenUnit = MSP11
        P12withinUnit = MSP12 * population_unit_SPN_modifier
        P12betweenUnit = MSP12

        # pfdSPNdSPN = "synapses/v1/trace_table.txt-DD-model-parameters.json"
        # pfdSPNiSPN = "synapses/v1/trace_table.txt-DI-model-parameters.json"
        pfdSPNdSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-DD-tmgaba-fit.json")
        pfdSPNiSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-DI-tmgaba-fit.json")
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
        P21withinUnit = MSP21 * population_unit_SPN_modifier
        P21betweenUnit = MSP21
        P22withinUnit = MSP22 * population_unit_SPN_modifier
        P22betweenUnit = MSP22

        pfiSPNdSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-ID-tmgaba-fit.json")
        pfiSPNiSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-II-tmgaba-fit.json")
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

        self.add_neuron_target(neuron_name="ChIN",
                               target_name="dSPN",
                               connection_type="Acetylcholine",
                               dist_pruning=None,
                               f1=None, soft_max=None, mu2=None, a3=None,  # SM 15
                               conductance=ChINgGABA,
                               parameter_file=pfChINdSPN,
                               mod_file="concACh",
                               channel_param_dictionary=None)

        # TEST SETTING THIS TO ACh (SHOULD BE GABA), will this change?
        # !!!

        self.add_neuron_target(neuron_name="ChIN",
                               target_name="iSPN",
                               connection_type="Acetylcholine",
                               dist_pruning=None,
                               f1=0.5, soft_max=10, mu2=10, a3=0.1,  # SM 12
                               conductance=ChINgGABA,
                               parameter_file=pfChINiSPN,
                               mod_file="concACh",
                               channel_param_dictionary=None)
        # ================================================================

        # We got an increasing connection distribution with distance, looks fishy
        # !!! Should be ACh, lets try set it to GABA and see if that changes things
        # --- trying same pruning as for ChIN to MSD2

        self.add_neuron_target(neuron_name="ChIN",
                               target_name="LTS",
                               connection_type="Acetylcholine",
                               dist_pruning=None,
                               f1=0.5, soft_max=None, mu2=10, a3=None,  # SM 12
                               conductance=ChINgACh,
                               parameter_file=pfChINLTS,
                               mod_file="concACh",  # !!! DOES NOT YET EXIST --- FIXME
                               channel_param_dictionary=None)

        self.add_neuron_target(neuron_name="ChIN",
                               target_name="FSN",
                               connection_type="Acetylcholine",
                               dist_pruning=None,
                               f1=0.5, soft_max=None, mu2=10, a3=None,  # SM 12
                               conductance=ChINgACh,
                               parameter_file=pfChINLTS,
                               mod_file="concACh",  # !!! DOES NOT YET EXIST --- FIXME
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
