import os
import glob
import collections
import json

import h5py

from snudda.CreateCubeMesh import CreateCubeMesh
from snudda.Neuron_morphology import NeuronMorphology
from snudda.init import SnuddaInit
from snudda.input import SnuddaInput
from snudda.load import SnuddaLoad
from snudda.simulate import SnuddaSimulate
from snudda.core import Snudda
import numpy as np
import timeit

import matplotlib.pyplot as plt


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        return json.JSONEncoder.default(self, obj)


class InputScaling(object):

    def __init__(self, network_dir, cellspec_dir):

        self.network_dir = network_dir
        self.cellspec_dir = cellspec_dir

        self.neuron_types = None
        self.neuron_id = None
        self.input_info = None
        self.neuron_info = None
        self.init_rng = None

        # TODO: Check baseline, close to threshold...
        # TODO: Check when at tonic activity, how sharp short burst can we get without depolarisation block
        self.frequency_range = np.arange(0, 5, 1) # np.arange(0, 50, 1)  # to 50 Hz, or maybe 100Hz...
        self.input_duration = 3.0  # 10.0
        self.max_time = self.input_duration * len(self.frequency_range)

        if not os.path.isdir(self.network_dir):
            os.mkdir(self.network_dir)

        assert os.path.isdir(self.cellspec_dir), f"Cellspec directory {self.cellspec_dir} does not exist."

        self.network_config_file_name = os.path.join(self.network_dir, "network-config.json")
        self.network_file = os.path.join(self.network_dir, "network-pruned-synapses.hdf5")
        self.input_config_file = os.path.join(self.network_dir, "input-config.json")
        self.input_spikes_file = os.path.join(self.network_dir, 'input.hdf5')

        self.output_spike_file = os.path.join(self.network_dir, 'output_spikes.txt')
        self.output_volt_file = os.path.join(self.network_dir, 'output_volt.txt')

        self.core = Snudda(self.network_dir)

    # Writes config files


    def setup_network(self):

        # TODO: num_replicas should be set by a parameter, it affects how many duplicates of each neuron
        # and thus how many steps we have between n_min and n_max number of inputs specified.
        config_def = self.create_network_config(num_replicas=10)

        print(f"Writing network config file to {self.network_config_file_name}")
        with open(self.network_config_file_name, "w") as f:
            json.dump(config_def, f, indent=2, cls=NumpyEncoder)

        CreateCubeMesh("data/mesh/InputTestMesh.obj", [0, 0, 0], 1e-3,
                       description="Mesh file used for Input Scaling")

        from snudda.place import SnuddaPlace
        from snudda.detect import SnuddaDetect
        from snudda.prune import SnuddaPrune

        position_file = os.path.join(self.network_dir, "network-neuron-positions.hdf5")
        SnuddaPlace(config_file=self.network_config_file_name, verbose=True).write_data_HDF5(position_file)

        SnuddaDetect(config_file=self.network_config_file_name,
                     position_file=position_file,
                     save_file=os.path.join(self.network_dir, "voxels", "network-putative-synapses.hdf5"))
        SnuddaPrune(work_history_file=os.path.join(self.network_dir, "log", "network-detect-worklog.hdf5"))

        # TODO: Skip placing neurons that will not receive any inputs or distribute any inputs

    def setup_input(self, input_type=None):

        synapse_density_cortical_input = "1.15*0.05/(1+np.exp(-(d-30e-6)/5e-6))"
        synapse_density_thalamic_input = "0.05*np.exp(-d/200e-6)"
        #  synapse_density_thalamic_input = "(d > 100e-6)*1"  # TEST!!

        cortical_SPN_synapse_parameter_file = "data/synapses/v2/M1RH_Analysis_190925.h5-parameters-MS.json"
        thalamic_SPN_synapse_parameter_file = "data/synapses/v2/TH_Analysis_191001.h5-parameters-MS.json"
        cortical_FS_synapse_parameter_file = "data/synapses/v2/M1RH_Analysis_190925.h5-parameters-FS.json"
        thalamic_FS_synapse_parameter_file = "data/synapses/v2/TH_Analysis_191001.h5-parameters-FS.json"
        cortical_ChIN_synapse_parameter_file = "data/synapses/v2/M1RH_Analysis_190925.h5-parameters-CHAT.json"
        thalamic_ChIN_synapse_parameter_file = "data/synapses/v2/TH_Analysis_191001.h5-parameters-CHAT.json"
        cortical_LTS_synapse_parameter_file = "data/synapses/v2/M1RH_Analysis_190925.h5-parameters-LTS.json"

        if input_type == 'cortical':
            synapse_density = synapse_density_cortical_input
            synapse_parameter_file = { "dspn": cortical_SPN_synapse_parameter_file,
                                       "ispn": cortical_SPN_synapse_parameter_file,
                                       "fs": cortical_FS_synapse_parameter_file,
                                       "lts": cortical_LTS_synapse_parameter_file,
                                       "chin": cortical_ChIN_synapse_parameter_file }
            print("Using cortical synapse density for input.")
        elif input_type == 'thalamic':
            synapse_density = synapse_density_thalamic_input
            synapse_parameter_file = { "dspn": thalamic_SPN_synapse_parameter_file,
                                       "ispn": thalamic_SPN_synapse_parameter_file,
                                       "fs": thalamic_FS_synapse_parameter_file,
                                       "chin": thalamic_ChIN_synapse_parameter_file }
            print("Using thalamic synapse density for input")
        else:
            synapse_density = "1"
            synapse_parameter_file = {}
            print("No density profile used for input.")

        self.create_input_config(input_config_file=self.input_config_file,
                                 input_type=input_type,
                                 input_frequency=list(self.frequency_range),  #[1.0],
                                 n_input_min=0,
                                 n_input_max=1000,
                                 synapse_conductance=0.5e-9,
                                 synapse_density=synapse_density,
                                 input_duration=self.input_duration,
                                 synapse_parameter_file=synapse_parameter_file)

        SnuddaInput(input_config_file=self.input_config_file,
                    hdf5_network_file=os.path.join(self.network_dir, 'network-pruned-synapses.hdf5'),
                    spike_data_filename=self.input_spikes_file,
                    time=self.max_time)

    def analyse_results(self):

        data = self.load_data()

        # We need to extract which cell models correspond to what traces, and also what input the got
        # Best to do it from the data in the files that were used.

        pass

    def load_data(self, skip_time=0.0):

        network_file = os.path.join(self.network_dir, "network-pruned-synapses.hdf5")
        network_info = SnuddaLoad(network_file)

        spike_data_file = os.path.join(self.network_dir, "output_spikes.txt")
        n_neurons = network_info.data["nNeurons"]
        spike_data = self.load_spike_data(spike_data_file, n_neurons)

        # We need to figure out what neuronID correspond to that morphologies
        # Then figure out what input frequencies the different runs had

        neuron_id_list = [x["neuronID"] for x in network_info.data["neurons"]]
        neuron_id_name_pairs = [(x["neuronID"], x["name"]) for x in network_info.data["neurons"]]

        # For each morphology-model we have a list of the run with that model
        neuron_lookup = dict()

        for neuron_id, neuron_name in neuron_id_name_pairs:
            if neuron_name in neuron_lookup:
                neuron_lookup[neuron_name].append(neuron_id)
            else:
                neuron_lookup[neuron_name] = [neuron_id]

        # Next identify number of inputs each run had
        input_config = self.load_input_config()

        n_inputs_lookup = dict()
        for neuron_label in input_config.keys():
            neuron_id = int(neuron_label)
            n_inputs = 0
            for input_types in input_config[neuron_label].keys():
                n_inputs += input_config[neuron_label]["nInputs"]

            n_inputs_lookup[neuron_id] = n_inputs

        frequency_data = dict()
        for neuron_name in neuron_lookup.keys():
            frequency_data[neuron_name] = dict()

        for neuron_id in neuron_id_list:
            n_inputs = n_inputs_lookup[neuron_id]
            frequency_data[neuron_id][n_inputs] = self.extract_spikes(spike_data=spike_data,
                                                                      config_data=input_config,
                                                                      neuron_id=neuron_id,
                                                                      skip_time=skip_time)

        # TO BE CONTINUED... next we need to plot this data.

        # !!! NEED TO RETURN DATA

    # TODO: We should set skip_time to 1 second, 0 for now while testing
    def extract_spikes(self, spike_data, config_data, neuron_id, skip_time=0.0):

        # This function checks the config_data to find out what time ranges to extract

        assert skip_time >= 0, "Time skipped at beginning of stimulus must be positive or zero"
        assert len(config_data[str(neuron_id)].values()) == 1, "Analysis can only handle one input type at a time"
        input_type = list(config_data[str(neuron_id)].keys())[0]

        cfg_data = config_data[str(neuron_id)][input_type]

        input_frequency = []
        output_frequency = []

        for start_time, end_time, input_freq in zip(cfg_data["start"], cfg_data["end"], cfg_data["frequency"]):

            assert start_time + skip_time < end_time, "Too large skip time, no data to analyse"
            spike_idx = np.where(start_time + skip_time <= spike_data[neuron_id] & spike_data[neuron_id] <= end_time)
            output_freq = len(spike_idx)

            input_frequency.append(input_freq)
            output_frequency.append(output_freq)

        input_frequency = np.array(input_frequency)
        output_frequency = np.array(output_frequency)

        return input_frequency, output_frequency

    def load_input_config(self):

        input_config_file = os.path.join(self.network_dir, "input-config.json")
        with open(input_config_file) as f:
            input_config = json.load(f)

        return input_config

    def load_spike_data(self, file_name, n_cells):

        data = np.genfromtxt(self.fileName, delimiter='\t')
        spike_time = data[:, 0] * 1e-3
        spike_id = data[:, 1].astype(int)

        spike_times = dict()
        for nid in range(0, n_cells):
            spike_times[nid] = []

        for sid, t in zip(spike_id, spike_time):
            spike_times[sid].append(t)

        for nid in range(0, n_cells):
            spike_times[nid] = np.array(sorted(spike_times[nid]))

        return spike_times

    # TODO: Extract spiking frequency (skip first second for each interval to let network settle)
    # TODO: Create summary graphs

    # This loops through all neuron directories in cellspec in preparation of writing a network config file
    def gather_all_neurons(self):
        all_neurons = collections.OrderedDict()

        neuron_type_dir = [d for d in glob.glob(os.path.join(self.cellspec_dir, '*')) if os.path.isdir(d)]

        self.neuron_types = []

        for ntd in neuron_type_dir:

            neuron_type = os.path.basename(os.path.normpath(ntd))
            self.neuron_types.append(neuron_type)

            neuron_dir = [d for d in glob.glob(os.path.join(ntd, '*')) if os.path.isdir(d)]
            neuron_ctr = 0

            for nd in neuron_dir:
                neuron_info = collections.OrderedDict()

                # Find neuron morphology swc file, obs currently assume lowercase(!)
                neuron_morph_list = glob.glob(os.path.join(nd, '*swc'))

                parameter_file = os.path.join(nd, "parameters.json")
                mechanism_file = os.path.join(nd, "mechanisms.json")
                modulation_file = os.path.join(nd, "modulation.json")  # Optional

                if len(neuron_morph_list) == 0:
                    assert (not os.path.isfile(parameter_file) and not os.path.isfile(mechanism_file)), \
                        f"Directory {nd} has parameter.json or mechanism.json but no swc file."

                    # No swc file, skipping directory
                    continue

                # Check if empty neuron_morph_list, or if more than one morphology
                assert len(neuron_morph_list) == 1, f"Should only be one swc file in {nd}"
                assert os.path.isfile(parameter_file), f"Missing parameter file {parameter_file}"
                assert os.path.isfile(mechanism_file), f"Missing mechanism file {mechanism_file}"

                neuron_info["morphology"] = neuron_morph_list[0]
                neuron_info["parameters"] = parameter_file
                neuron_info["mechanisms"] = mechanism_file

                # Modulation file is optional
                if os.path.isfile(modulation_file):
                    neuron_info["modulation"] = modulation_file

                neuron_ctr += 1
                neuron_name = f"{neuron_type}_{neuron_ctr}"

                all_neurons[neuron_name] = neuron_info

            if neuron_ctr > 0:
                print(f"Found {neuron_ctr} neurons in {ntd}")

        return all_neurons

    @staticmethod
    def has_axon(swc_file):
        nm = NeuronMorphology(swc_filename=swc_file)

        return len(nm.axon) > 0

    def create_network_config(self, num_replicas=10, random_seed=None):

        config_def = collections.OrderedDict()
        config_def["RandomSeed"], self.init_rng = SnuddaInit.setup_random_seeds(random_seed)

        volume_def = collections.OrderedDict()
        vol_name = "InputTest"
        volume_def[vol_name] = collections.OrderedDict()

        volume_def[vol_name]["type"] = "mesh"
        volume_def[vol_name]["dMin"] = 15e-6
        volume_def[vol_name]["meshFile"] = "data/mesh/InputTestMesh.obj"
        volume_def[vol_name]["meshBinWidth"] = 100e-6

        config_def["Volume"] = volume_def
        config_def["Connectivity"] = dict()  # Unconnected

        neuron_def = self.gather_all_neurons()

        fake_axon_density = ["r", "1", 10e-6]

        for n in neuron_def.keys():
            neuron_def[n]["num"] = num_replicas
            neuron_def[n]["volumeID"] = vol_name
            neuron_def[n]["rotationMode"] = "random"
            neuron_def[n]["hoc"] = None

            if not self.has_axon(neuron_def[n]["morphology"]):
                print(f"Morphology {neuron_def[n]['morphology']} has no axon, faking it.")
                # We will have no connections in this test network, so add empty density
                neuron_def[n]["axonDensity"] = fake_axon_density

        config_def["Neurons"] = neuron_def

        return config_def

    def load_network_info(self, network_file):

        self.neuron_info = SnuddaLoad(network_file).data["neurons"]

    def collect_neurons(self):

        if not self.neuron_info:
            self.load_network_info(network_file=self.network_file)

        neuron_id = np.array([x["neuronID"] for x in self.neuron_info])
        neuron_type = np.array([x["type"] for x in self.neuron_info])

        neuron_sets = dict()

        for nt in set(neuron_type):
            idx = np.where(neuron_type == nt)
            neuron_id_list = np.sort(neuron_id[idx])

            neuron_sets[nt] = neuron_id_list

        return neuron_sets

    # synapse_density and synapse_conductance can be either values (then same for all neurons)
    # or dictionaries, with neuron_type as key.

    def create_input_config(self,
                            input_config_file,
                            input_type,
                            n_input_min, n_input_max, input_frequency,
                            synapse_density,
                            synapse_conductance,
                            synapse_parameter_file,
                            input_duration=10.0):

        neuron_sets = self.collect_neurons()
        n_inputs = dict()

        for neuron_type in neuron_sets:

            # For each neuron model we will have num_replicas copies (see other part of code), and this
            # will determine how many steps we have between n_input_min and n_input_max

            neuron_id_list = neuron_sets[neuron_type]
            num_range = np.linspace(n_input_min, n_input_max, num=len(neuron_id_list)).astype(int)

            for neuron_id, num_input in zip(neuron_id_list, num_range):

                if type(synapse_density) == dict:
                    sd = synapse_density[neuron_type]
                else:
                    sd = synapse_density

                if type(synapse_parameter_file) == dict:
                    if neuron_type not in synapse_parameter_file:
                        print(f"No parameter file for {neuron_type} {input_type} input, excluding from run.")
                        continue

                    spf = synapse_parameter_file[neuron_type]
                else:
                    spf = synapse_parameter_file

                if type(synapse_conductance) == dict:
                    sc = synapse_conductance[neuron_type]
                else:
                    sc = synapse_conductance

                self.add_input(input_target=neuron_id,
                               input_type=input_type,
                               input_frequency=input_frequency,
                               input_duration=input_duration,
                               input_density=sd,
                               num_input=num_input,
                               input_conductance=sc,
                               synapse_parameter_file=spf)

        with open(input_config_file, "w") as f:
            json.dump(self.input_info, f, indent=4, cls=NumpyEncoder)

    def add_input(self, input_target, input_type, input_frequency, input_duration,
                  input_density, num_input, input_conductance,
                  synapse_parameter_file):

        if type(input_target) != str:
            input_target = str(input_target)

        if not self.input_info:
            self.input_info = collections.OrderedDict()

        assert type(input_frequency) == list, f"add_input: input_freq should be a list: {input_frequency}"
        n_steps = len(input_frequency)
        input_start = input_duration * np.arange(0, n_steps)
        input_end = input_duration * np.arange(1, n_steps+1)

        self.input_info[input_target] = collections.OrderedDict()
        self.input_info[input_target][input_type] = collections.OrderedDict()

        self.input_info[input_target][input_type]["generator"] = "poisson"
        self.input_info[input_target][input_type]["type"] = "AMPA_NMDA"
        self.input_info[input_target][input_type]["synapseDensity"] = input_density
        self.input_info[input_target][input_type]["nInputs"] = num_input
        self.input_info[input_target][input_type]["frequency"] = input_frequency
        self.input_info[input_target][input_type]["start"] = input_start
        self.input_info[input_target][input_type]["end"] = input_end
        self.input_info[input_target][input_type]["populationUnitCorrelation"] = 0.0
        self.input_info[input_target][input_type]["jitter"] = 0.0
        self.input_info[input_target][input_type]["conductance"] = input_conductance
        self.input_info[input_target][input_type]["modFile"] = "tmGlut"
        self.input_info[input_target][input_type]["parameterFile"] = synapse_parameter_file

    def plot_generated_input(self, num_bins=50):
        # This function just checks that we have reasonable spikes generated

        input_spike_data = h5py.File(self.input_spikes_file, 'r')
        network_data = h5py.File(self.network_file, 'r')

        neuron_type = np.array([x.decode().split("_")[0].lower() for x in network_data["network/neurons/name"]])
        neuron_id = np.array([x for x in network_data["network/neurons/neuronID"]])

        # Plot the input spikes
        for nt in set(neuron_type):
            neuron_idx = np.where(neuron_type == nt)
            distance_to_soma = dict()

            if len(neuron_idx) == 0:
                continue

            fig, ax = plt.subplots()

            for nid in neuron_id[neuron_idx]:
                for input_type in input_spike_data["input"][str(nid)]:
                    spikes = input_spike_data["input"][str(nid)][input_type]["spikes"][:].ravel()
                    spikes = spikes[spikes >= 0]  # Negative -1 is filler values, remove them.
                    ax.hist(spikes, num_bins, histtype="step")

                    if input_type not in distance_to_soma:
                        distance_to_soma[input_type] = []

                    distance_to_soma[input_type].append(input_spike_data["input"][str(nid)][input_type]["distanceToSoma"][:])

            plt.title(f"Input to {nt}")
            plt.xlabel("Time (s)")
            plt.ylabel("Count")
            plt.ion()
            plt.show()
            plt.pause(0.001)

            fig_dir = os.path.join(self.network_dir, "figures")
            if not os.path.isdir(fig_dir):
                os.mkdir(fig_dir)

            fig_name1 = os.path.join(self.network_dir, "figures", f"{nt}-binned-input.pdf")
            plt.savefig(fig_name1)

            fig2, ax2 = plt.subplots()
            leg = []
            for input_type in distance_to_soma:
                leg.append(input_type)
                ax2.hist(np.concatenate(distance_to_soma[input_type])*1e6, 50, histtype="step")

            plt.legend(leg)
            plt.title(f"Synapses onto {nt}")
            plt.xlabel(f"Distance to soma (micrometers)")
            plt.ylabel(f"Count")
            plt.show()
            plt.pause(0.001)

            fig_name2 = os.path.join(self.network_dir, "figures", f"{nt}-input-distance-to-soma.pdf")
            plt.savefig(fig_name2)

        input_spike_data.close()
        network_data.close()


    def simulate(self):

        from neuron import h  # , gui
        start = timeit.default_timer()

        pc = h.ParallelContext()

        sim = SnuddaSimulate(network_file=self.network_file,
                             input_file=self.input_spikes_file,
                             log_file=None,  # Set log file?
                             verbose=True)

        sim.add_external_input()
        sim.check_memory_status()

        sim.add_recording()

        t_sim = self.max_time * 1000  # Convert from s to ms for Neuron simulator

        sim.check_memory_status()
        print("Running simulation for " + str(t_sim) + " ms.")
        sim.run(t_sim)  # In milliseconds

        print("Simulation done, saving output")
        sim.write_spikes(self.output_spike_file )
        sim.write_voltage(self.output_volt_file)

        stop = timeit.default_timer()
        if sim.pc.id() == 0:
            print("Program run time: " + str(stop - start))


if __name__ == "__main__":
    from argparse import ArgumentParser, RawTextHelpFormatter

    parser = ArgumentParser("Input Scaling", formatter_class=RawTextHelpFormatter)
    parser.add_argument("action", choices=["setup", "simulate", "analyse"], help="Action to run.")
    parser.add_argument("networkPath", help="Network path")
    parser.add_argument("cellspecspath", help="Cellspecs path")

    args = parser.parse_args()

    # TODO: Let the user choose input type, duration for each "run", frequency range, number of input range

    input_scaling = InputScaling(args.networkPath, args.cellspecspath)

    if args.action == "setup":
        input_scaling.setup_network()
        input_scaling.setup_input(input_type="thalamic")
        #input_scaling.setup_input(input_type="cortical")

    elif args.action == "simulate":
        print("Run simulation...")
        print("Tip, to run in parallel on your local machine use: "
              "mpiexec -n 4 python3 tuning/input_scaling.py simulate <yournetworkhere>")
        input_scaling.simulate()

    elif args.action == "analyse":
        input_scaling.plot_generated_input()

    else:
        print(f"Unknown action {args.action}")


    # python3 tuning/input_scaling.py setup networks/input_scaling_v1/ data/cellspecs-v2/
    # mpiexec -n 4 python3 tuning/input_scaling.py simulate networks/input_scaling_v1/ data/cellspecs-v2/ < input.txt &> output-tuning.txt &