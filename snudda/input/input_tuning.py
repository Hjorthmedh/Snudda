import ast
import os
import sys
import glob
import collections
import json

import h5py

from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.place.create_cube_mesh import create_cube_mesh
from snudda.neurons.neuron_morphology import NeuronMorphology
from snudda.utils import SnuddaLoadNetworkSimulation
from snudda.utils.reposition_neurons import RepositionNeurons
from snudda.init.init import SnuddaInit
from snudda.input.input import SnuddaInput
from snudda.utils.load import SnuddaLoad
from snudda.simulate.simulate import SnuddaSimulate
from snudda.core import Snudda
import numpy as np
import timeit

import matplotlib.pyplot as plt

from snudda.utils.snudda_path import snudda_isdir, snudda_parse_path, snudda_simplify_path


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        return json.JSONEncoder.default(self, obj)


class InputTuning(object):

    def __init__(self, network_path):

        self.network_path = network_path
        self.neurons_path = None

        self.neuron_types = None
        self.neuron_id = None
        self.input_info = None
        self.neuron_info = None
        self.init_rng = None

        # TODO: Check baseline, close to threshold...
        # TODO: Check when at tonic activity, how sharp short burst can we get without depolarisation block
        self.frequency_range = None
        self.input_duration = None
        self.max_time = None  # self.input_duration * len(self.frequency_range)

        if not os.path.isdir(self.network_path):
            os.mkdir(self.network_path)

        self.network_config_file_name = os.path.join(self.network_path, "network-config.json")
        self.network_file = os.path.join(self.network_path, "network-synapses.hdf5")
        self.input_config_file = os.path.join(self.network_path, "input_config.json")
        self.input_spikes_file = os.path.join(self.network_path, 'input.hdf5')

        self.core = Snudda(self.network_path)

    # Writes config files

    def setup_network(self, neurons_path=None, num_replicas=10, neuron_types=None,
                      parameter_key=None, morphology_key=None, modulation_key=None,
                      single_neuron_path=None):

        if not morphology_key and not parameter_key and not modulation_key:
            all_combinations = True
        else:
            all_combinations = False

        # TODO: num_replicas should be set by a parameter, it affects how many duplicates of each neuron
        # and thus how many steps we have between n_min and n_max number of inputs specified.
        config_def = self.create_network_config(neurons_path=neurons_path,
                                                num_replicas=num_replicas,
                                                neuron_types=neuron_types,
                                                single_neuron_path=single_neuron_path,
                                                parameter_key=parameter_key,
                                                morphology_key=morphology_key,
                                                modulation_key=modulation_key,
                                                all_combinations=all_combinations)

        print(f"Writing network config file to {self.network_config_file_name}")
        with open(self.network_config_file_name, "w") as f:
            json.dump(config_def, f, indent=2, cls=NumpyEncoder)

        create_cube_mesh(os.path.join("data", "mesh", "InputTestMesh.obj"), [0, 0, 0], 1e-3,
                         description="Mesh file used for Input Scaling")

        # Write the neurons path to file
        self.write_tuning_info()

        from snudda.place.place import SnuddaPlace
        from snudda.detect.detect import SnuddaDetect
        from snudda.detect.prune import SnuddaPrune

        sp = SnuddaPlace(network_path=self.network_path)
        sp.parse_config(resort_neurons=False)  # By not resorting neurons, we have original order
        sp.write_data()

        sd = SnuddaDetect(network_path=self.network_path)
        sd.detect()

        sp = SnuddaPrune(network_path=self.network_path)
        sp.prune()

        # TODO: Skip placing neurons that will not receive any inputs or distribute any inputs

    def setup_input(self, input_type=None, num_input_min=100, num_input_max=1000,
                    input_duration=10,
                    input_frequency_range=None):

        if not input_frequency_range:
            input_frequency_range = [1.0]

        self.frequency_range = np.array(input_frequency_range)
        self.input_duration = input_duration
        self.max_time = self.input_duration * len(self.frequency_range)

        synapse_density_cortical_input = "1.15*0.05/(1+exp(-(d-30e-6)/5e-6))"
        synapse_density_thalamic_input = "0.05*exp(-d/200e-6)"
        #  synapse_density_thalamic_input = "(d > 100e-6)*1"  # TEST!!

        cortical_SPN_synapse_parameter_file = "$DATA/synapses/striatum/M1RH_Analysis_190925.h5-parameters-MS.json"
        thalamic_SPN_synapse_parameter_file = "$DATA/synapses/striatum/TH_Analysis_191001.h5-parameters-MS.json"
        cortical_FS_synapse_parameter_file = "$DATA/synapses/striatum/M1RH_Analysis_190925.h5-parameters-FS.json"
        thalamic_FS_synapse_parameter_file = "$DATA/synapses/striatum/TH_Analysis_191001.h5-parameters-FS.json"
        cortical_ChIN_synapse_parameter_file = "$DATA/synapses/striatum/M1RH_Analysis_190925.h5-parameters-CHAT.json"
        thalamic_ChIN_synapse_parameter_file = "$DATA/synapses/striatum/TH_Analysis_191001.h5-parameters-CHAT.json"
        cortical_LTS_synapse_parameter_file = "$DATA/synapses/striatum/M1RH_Analysis_190925.h5-parameters-LTS.json"

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
            print("No density profile used for input    .")

        self.create_input_config(input_config_file=self.input_config_file,
                                 input_type=input_type,
                                 input_frequency=list(self.frequency_range),  #[1.0],
                                 n_input_min=num_input_min,
                                 n_input_max=num_input_max,
                                 synapse_conductance=0.5e-9,
                                 synapse_density=synapse_density,
                                 input_duration=self.input_duration,
                                 synapse_parameter_file=synapse_parameter_file)

        si = SnuddaInput(input_config_file=self.input_config_file,
                         hdf5_network_file=os.path.join(self.network_path, 'network-synapses.hdf5'),
                         spike_data_filename=self.input_spikes_file,
                         time=self.max_time)
        si.generate()

        # Info we need to run right duration of simulation
        self.write_tuning_info()

    def analyse_results(self, show_plots=False):

        frequency_data, voltage_data = self.load_data()
        self.plot_frequency_data(frequency_data, show_plots=show_plots)
        self.plot_frequency_data_alt(frequency_data, show_plots=show_plots)
        self.plot_volt_data(voltage_data, show_plots=show_plots)

        print(f"To plot traces:\n" 
              f"python3 plotting/Network_plot_traces.py {self.network_path}output_volt.txt " 
              f"{self.network_path}network-synapses.hdf5 ")

    def load_data(self, skip_time=0.0):

        network_file = os.path.join(self.network_path, "network-synapses.hdf5")
        network_info = SnuddaLoad(network_file)

        input_data = h5py.File(self.input_spikes_file, "r")

        output_data_loader = SnuddaLoadNetworkSimulation(network_path=self.network_path)
        spike_data = output_data_loader.get_spikes()

        # cell_id = output_data_loader.get_id_of_neuron_type()
        volt, time = output_data_loader.get_voltage()

        # We need to figure out what neuronID correspond to that morphologies
        # Then figure out what input frequencies the different runs had

        neuron_id_list = output_data_loader.get_id_of_neuron_type()
        neuron_name_list = output_data_loader.get_neuron_name()

        # For each morphology-model we have a list of the run with that model
        neuron_id_lookup = dict()

        for neuron_id, neuron_name in zip(neuron_id_list, neuron_name_list):
            if neuron_name in neuron_id_lookup:
                neuron_id_lookup[neuron_name].append(neuron_id)
            else:
                neuron_id_lookup[neuron_name] = [neuron_id]

        # Next identify number of inputs each run had
        input_config = self.load_input_config()

        n_inputs_lookup = dict()

        for neuron_label in input_data["input"]:
            neuron_id = int(neuron_label)
            n_inputs = 0

            for input_type in input_data["input"][neuron_label]:
                n_inputs += input_data["input"][neuron_label][input_type]["spikes"].shape[0]

            n_inputs_lookup[neuron_id] = n_inputs

        frequency_data = dict()
        voltage_data = dict()
        
        for neuron_name in neuron_id_lookup.keys():
            frequency_data[neuron_name] = dict()
            voltage_data[neuron_name] = dict()

        for neuron_name in neuron_name_list:
            for neuron_id in neuron_id_lookup[neuron_name]:

                if neuron_id not in n_inputs_lookup:
                    print(f"No inputs for neuron_id={neuron_id}, ignoring. Please update your setup.")
                    continue

                n_inputs = n_inputs_lookup[neuron_id]
                frequency_data[neuron_name][n_inputs] = self.extract_spikes(spike_data=spike_data,
                                                                            config_data=input_config,
                                                                            neuron_id=neuron_id,
                                                                            skip_time=skip_time)
                voltage_data[neuron_name][n_inputs] = self.extract_voltage(volt=volt, 
                                                                           time=time,
                                                                           config_data=input_config,
                                                                           neuron_id=neuron_id,
                                                                           skip_time=skip_time)

        # TODO: Load voltage trace and warn for depolarisation blocking

        return frequency_data, voltage_data

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


            assert start_time + skip_time < end_time, \
                f"Too large skip time, no data to analyse. start {start_time} + skip {skip_time} < end {end_time}"

            try:
                spike_idx = np.where((start_time + skip_time <= spike_data[neuron_id])
                                     & (spike_data[neuron_id] <= end_time))[0]
                output_freq = len(spike_idx) / (end_time - start_time)
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

            input_frequency.append(input_freq)
            output_frequency.append(output_freq)

        input_frequency = np.array(input_frequency)
        output_frequency = np.array(output_frequency)

        return input_frequency, output_frequency, input_type

    def load_input_config(self):

        input_config_file = os.path.join(self.network_path, "input_config.json")
        with open(input_config_file) as f:
            input_config = json.load(f, object_pairs_hook=collections.OrderedDict)

        return input_config

    def extract_voltage(self, volt, time, config_data, neuron_id, skip_time=0.0):

        assert skip_time >= 0
        assert len(config_data[str(neuron_id)].values()) == 1

        input_type = list(config_data[str(neuron_id)].keys())[0]
        cfg_data = config_data[str(neuron_id)][input_type]

        input_frequency = []
        mean_voltage = []
        max_voltage = []

        for start_time, end_time, input_freq in zip(cfg_data["start"], cfg_data["end"], cfg_data["frequency"]):

            assert start_time + skip_time < end_time, "Too large skip time, no data to analyse"

            idx = np.where((start_time + skip_time <= time) & (time <= end_time))[0]
            v = volt[neuron_id][idx]

            input_frequency.append(input_freq)
            mean_voltage.append(np.mean(v))
            max_voltage.append(np.max(v))
  
        input_frequency = np.array(input_frequency)
        mean_voltage = np.array(mean_voltage)
        max_voltage = np.array(max_voltage)

        return input_frequency, mean_voltage, max_voltage, input_type

    def load_spike_data(self, file_name, n_cells):

        data = np.genfromtxt(file_name, delimiter='\t')
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

    def plot_volt_data(self, volt_data, show_plots=True):

        for neuron_name in volt_data:
            fig, ax = plt.subplots()
            legend_text = []
            input_type_all = None

            cmap = plt.get_cmap('tab20', len(volt_data[neuron_name]))
            ax.set_prop_cycle('color', [cmap(i) for i in range(0, len(volt_data[neuron_name]))])

            for num_input in volt_data[neuron_name]:
                input_freq, mean_voltage, max_voltage, input_type = volt_data[neuron_name][num_input]

                if input_type_all:
                    assert input_type == input_type_all, "All input types must be the same for neuron"
                else:
                    input_type_all = input_type

                legend_text.append(f"n={num_input}")
                ax.plot(input_freq, mean_voltage)

            plt.title(f"{neuron_name} receiving {input_type_all} input")
            plt.xlabel("Input frequency (per synapse)")
            plt.ylabel("Mean voltage")
            ax.legend(legend_text)

            if show_plots:
                plt.ion()
                plt.show()
                plt.pause(0.001)

            fig_name = os.path.join(self.network_path, "figures", f"input-scaling-freq-{neuron_name}-mean-voltage.pdf")
            if not os.path.exists(os.path.dirname(fig_name)):
                os.mkdir(os.path.dirname(fig_name))

            plt.savefig(fig_name, dpi=300)

            if not show_plots:
                plt.close()

    # TODO: Extract spiking frequency (skip first second for each interval to let network settle)
    # TODO: Create summary graphs

    def plot_frequency_data(self, frequency_data, show_plots=True):

        for neuron_name in frequency_data:
            fig, ax = plt.subplots()
            legend_text = []
            input_type_all = None

            cmap = plt.get_cmap('tab20', len(frequency_data[neuron_name]))
            ax.set_prop_cycle('color', [cmap(i) for i in range(0, len(frequency_data[neuron_name]))])

            for num_input in frequency_data[neuron_name]:
                input_freq, output_freq, input_type = frequency_data[neuron_name][num_input]

                if input_type_all:
                    assert input_type == input_type_all, "All input types must be the same for neuron"
                else:
                    input_type_all = input_type

                legend_text.append(f"n={num_input}")
                ax.plot(input_freq, output_freq)

            plt.title(f"{neuron_name} receiving {input_type_all} input")
            plt.xlabel("Input frequency (per synapse)")
            plt.ylabel("Firing frequency")
            ax.legend(legend_text)

            if show_plots:
                plt.ion()
                plt.show()
                plt.pause(0.001)

            fig_name = os.path.join(self.network_path, "figures", f"input-scaling-freq-{neuron_name}.pdf")
            if not os.path.exists(os.path.dirname(fig_name)):
                os.mkdir(os.path.dirname(fig_name))

            plt.savefig(fig_name, dpi=300)

            if not show_plots:
                plt.close()

    def plot_frequency_data_alt(self, frequency_data, show_plots=True):

        _freq_data = dict()
        _all_num_inputs = []
        _all_input_freq = []

        input_type_all = None

        # Put everything in a dictionary structure first
        for neuron_name in frequency_data:
            _freq_data[neuron_name] = dict()

            for num_inputs in frequency_data[neuron_name]:
                input_freq, output_freq, input_type = frequency_data[neuron_name][num_inputs]

                if input_type_all is not None:
                    assert input_type == input_type_all, "All neurons must have same input type for this plot"
                else:
                    input_type_all = input_type

                for in_freq, out_freq in zip(input_freq, output_freq):
                    if in_freq not in _freq_data[neuron_name]:
                        _freq_data[neuron_name][in_freq] = dict()

                    _freq_data[neuron_name][in_freq][num_inputs] = out_freq
                    _all_input_freq.append(in_freq)

                _all_num_inputs.append(num_inputs)

        # Extract the relevant data, make sure it is sorted in right order also
        input_freq_list = np.array(sorted(list(set(_all_input_freq))))
        num_input_list = np.array(sorted(list(set(_all_num_inputs))))

        freq_data = dict()
        for neuron_name in frequency_data:
            freq_data[neuron_name] = dict()
            for input_freq in input_freq_list:
                num_in_list = []  # Number of inputs
                out_f_list = []
                for num_inputs in num_input_list:
                    if input_freq in _freq_data[neuron_name] and \
                            num_inputs in _freq_data[neuron_name][input_freq]:
                        num_in_list.append(num_inputs)
                        out_f_list.append(_freq_data[neuron_name][input_freq][num_inputs])

                freq_data[neuron_name][input_freq] = np.array(num_in_list), np.array(out_f_list)

        for neuron_name in frequency_data:
            fig, ax = plt.subplots()
            legend_text = []

            cmap = plt.get_cmap("tab20", len(freq_data[neuron_name]))
            ax.set_prop_cycle('color', [cmap(i) for i in range(0, len(freq_data[neuron_name]))])

            for input_freq in freq_data[neuron_name]:
                num_input, output_freq = freq_data[neuron_name][input_freq]
                legend_text.append(f"{input_freq} Hz input")
                ax.plot(num_input, output_freq)

            plt.title(f"{neuron_name} receiving {input_type_all} input")
            plt.xlabel("Number of synapses")
            plt.ylabel("Firing frequency")
            ax.legend(legend_text)

            if show_plots:
                plt.ion()
                plt.show()
                plt.pause(0.001)

            fig_name = os.path.join(self.network_path, "figures", f"input-scaling-ninputs-{neuron_name}.pdf")
            if not os.path.exists(os.path.dirname(fig_name)):
                os.mkdir(os.path.dirname(fig_name))

            plt.savefig(fig_name, dpi=300)

            if not show_plots:
                plt.close()

    def get_neuron_info(self, neuron_path):

        neuron_path = snudda_parse_path(neuron_path)

        neuron_info = collections.OrderedDict()

        # Find neuron morphology swc file, obs currently assume lowercase(!)
        neuron_morph = SnuddaInit.get_morphologies(neuron_path)

        parameter_file = os.path.join(neuron_path, "parameters.json")
        mechanism_file = os.path.join(neuron_path, "mechanisms.json")
        modulation_file = os.path.join(neuron_path, "modulation.json")  # Optional
        meta_file = os.path.join(neuron_path, "meta.json")

        # Check if empty neuron_morph_list, or if more than one morphology
        assert os.path.isfile(parameter_file), f"Missing parameter file {parameter_file}"
        assert os.path.isfile(mechanism_file), f"Missing mechanism file {mechanism_file}"

        neuron_info["morphology"] = snudda_simplify_path(neuron_morph)
        neuron_info["parameters"] = snudda_simplify_path(parameter_file)
        neuron_info["mechanisms"] = snudda_simplify_path(mechanism_file)
        neuron_info["meta"] = snudda_simplify_path(meta_file)

        # Modulation file is optional
        if os.path.isfile(modulation_file):
            neuron_info["modulation"] = snudda_simplify_path(modulation_file)

        return neuron_info

    # This loops through all single neuron directories in neurons_path
    # in preparation of writing a network config file

    def gather_all_neurons(self, neuron_types=None, all_combinations=True):
        all_neurons = collections.OrderedDict()

        assert snudda_isdir(self.neurons_path), f"Neurons directory {self.neurons_path} does not exist."

        neuron_type_dir = [d for d in glob.glob(os.path.join(snudda_parse_path(self.neurons_path), '*'))
                           if snudda_isdir(d)]

        self.neuron_types = []

        if neuron_types is not None:
            if type(neuron_types) == str:
                neuron_types = [neuron_types]
            neuron_types = [x.lower() for x in neuron_types]

        for ntd in neuron_type_dir:

            neuron_type = os.path.basename(os.path.normpath(ntd))

            if neuron_types is not None:
                if neuron_type.lower() not in neuron_types:
                    print(f"Skipping neuron type {neuron_type}")
                    continue

            self.neuron_types.append(neuron_type)

            neuron_dir = [d for d in glob.glob(os.path.join(ntd, '*')) if os.path.isdir(d)]
            neuron_ctr = 0

            for nd in neuron_dir:
                neuron_info = self.get_neuron_info(nd)

                if all_combinations:
                    assert "meta" in neuron_info and neuron_info["meta"], \
                        f"meta.json required for all_combinations=True. {os.path.dirname(neuron_info['parameters'])}"

                    neuron_info_combination_list = self.get_all_combinations(neuron_info)
                    for ni in neuron_info_combination_list:
                        n_name = os.path.basename(os.path.dirname(ni["parameters"]))
                        param_key = ni["parameterKey"]
                        morph_key = ni["morphologyKey"]
                        short_name = n_name[:min(10, len(n_name))]
                        neuron_name = f"{neuron_type}_{short_name}_{param_key}_{morph_key}".replace("-","_")
                        all_neurons[neuron_name] = ni
                        neuron_ctr += 1
                else:
                    neuron_name = os.path.basename(os.path.dirname(neuron_info["parameters"]))
                    neuron_ctr += 1

                    all_neurons[neuron_name] = neuron_info

            if neuron_ctr > 0:
                print(f"Found {neuron_ctr} neuron models in {ntd}")

        assert len(all_neurons) > 0, (f"No neurons selected. Did you specify an incorrect neuronType? {neuron_types}"
                                      f"\nSee skipped neurons above error message for available ones.")

        return all_neurons

    @staticmethod
    def get_all_combinations(neuron_info):

        assert "meta" in neuron_info and neuron_info["meta"]

        pm_list = []

        with open(snudda_parse_path(neuron_info["parameters"]), "r") as pf:
            param_data = json.load(pf)

        with open(snudda_parse_path(neuron_info["meta"]), "r") as mf:
            meta_data = json.load(mf)

        for p_idx, p_key in enumerate(param_data):
            assert p_key in meta_data, f"parameter key {p_key} missing in {neuron_info['meta_file']}"
            for m_idx, m_key in enumerate(meta_data[p_key]):
                ni = neuron_info.copy()
                ni["morphologyKey"] = m_key
                ni["parameterKey"] = p_key
                ni["morphology_file"] = os.path.join(ni["morphology"], meta_data[p_key][m_key]["morphology"])
                ni["parameter_id"] = p_idx
                ni["morphology_id"] = m_idx

                pm_list.append(ni)

        return pm_list

    @staticmethod
    def has_axon(neuron_info):

        nm = NeuronPrototype(neuron_name="JJJ",
                             neuron_path=None,
                             morphology_path=neuron_info["morphology"],
                             parameter_path=neuron_info["parameters"],
                             mechanism_path=neuron_info["mechanisms"],
                             virtual_neuron=False,
                             axon_stump_id_flag=False)
        nm.instantiate()
        return nm.all_have_axon()

    def create_network_config(self,
                              neurons_path=None,
                              num_replicas=10,
                              random_seed=None,
                              neuron_types=None,
                              single_neuron_path=None,
                              parameter_key=None,
                              morphology_key=None,
                              modulation_key=None,
                              all_combinations=True):

        """ Create a network with num_replicas number of replicas of each neuron found in
            neuron_path (e.g. Snudda/data/neurons). Alternatively if only one neuron model
            is required, use single_neuron_path pointing to neuron
            (e.g. Snudda/data/neurons/striatum/fs/fs_model_1)
            """

        self.neurons_path = neurons_path

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

        if single_neuron_path:
            # Override and only get one neuron
            assert type(neuron_types) == str, "neuron_types must be string if single_neuron_path set"
            neuron_def = collections.OrderedDict()
            neuron_def[neuron_types] = self.get_neuron_info(single_neuron_path)

            if parameter_key:
                neuron_def[neuron_types]["parameterKey"] = parameter_key

            if morphology_key:
                neuron_def[neuron_types]["morphologyKey"] = morphology_key

            if modulation_key:
                neuron_def[neuron_types]["modulationKey"] = modulation_key

        else:
            neuron_def = self.gather_all_neurons(neuron_types=neuron_types, all_combinations=all_combinations)

        fake_axon_density = ["r", "1", 10e-6]

        for n in neuron_def.keys():
            neuron_def[n]["num"] = num_replicas
            neuron_def[n]["volumeID"] = vol_name
            neuron_def[n]["rotationMode"] = "random"
            neuron_def[n]["hoc"] = None

            if not self.has_axon(neuron_def[n]):
                print(f"One or more of morphologies {neuron_def[n]['morphology']} has no axon, faking it.")
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

        assert n_input_min > 0, "No point using n_input_min=0, please instead use input_frequency 0."

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

            fig_dir = os.path.join(self.network_path, "figures")
            if not os.path.isdir(fig_dir):
                os.mkdir(fig_dir)

            fig_name1 = os.path.join(self.network_path, "figures", f"{nt}-binned-input.pdf")
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

            fig_name2 = os.path.join(self.network_path, "figures", f"{nt}-input-distance-to-soma.pdf")
            plt.savefig(fig_name2)

        input_spike_data.close()
        network_data.close()

    def simulate(self, mech_dir=None):

        from snudda.core import Snudda
        Snudda.compile_mechanisms(mech_dir=mech_dir)

        # Get info so we can set max_time correctly
        self.read_tuning_info()

        from mpi4py import MPI
        from neuron import h  # , gui
        start = timeit.default_timer()

        pc = h.ParallelContext()

        sim = SnuddaSimulate(network_file=self.network_file,
                             input_file=self.input_spikes_file)
        sim.setup()
        sim.add_external_input()
        sim.check_memory_status()

        sim.add_recording()

        t_sim = self.max_time * 1000  # Convert from s to ms for Neuron simulator

        sim.check_memory_status()
        print("Running simulation for " + str(t_sim) + " ms.")
        sim.run(t_sim)  # In milliseconds

        print("Simulation done, saving output")
        sim.write_output()

        stop = timeit.default_timer()
        if sim.pc.id() == 0:
            print(f"Program run time: {stop - start:.1f}s")

    def read_tuning_info(self):
        tuning_info_file = os.path.join(self.network_path, "tuning-info.json")

        if not os.path.exists(tuning_info_file):
            print("No tuning info file exists.")
            return

        try:
            with open(tuning_info_file, 'rt') as f:
                tuning_meta_data = json.load(f, object_pairs_hook=collections.OrderedDict)

            # max_time is the important one, we want to make sure we simulate correct duration without having the user
            # provide the parameter twice
            self.max_time = tuning_meta_data["MaxTime"]
            self.input_duration = tuning_meta_data["InputDuration"]
            self.frequency_range = tuning_meta_data["FrequencyRange"]
            self.neurons_path = tuning_meta_data["NeuronsDirectory"]
        except:
            print(f"Failed to read {tuning_info_file}")

    def write_tuning_info(self):
        tuning_meta_data = collections.OrderedDict()
        tuning_meta_data["InputDuration"] = self.input_duration
        tuning_meta_data["MaxTime"] = self.max_time
        tuning_meta_data["FrequencyRange"] = self.frequency_range
        tuning_meta_data["NeuronsDirectory"] = self.neurons_path

        tuning_info_file = os.path.join(self.network_path, "tuning-info.json")
        with open(tuning_info_file, "wt") as f:
            json.dump(tuning_meta_data, f, indent=4, cls=NumpyEncoder)


if __name__ == "__main__":

    if '-python' in sys.argv:
        print("Network_simulate.py called through nrniv, fixing arguments")
        pythonidx = sys.argv.index('-python')
        if len(sys.argv) > pythonidx:
            sys.argv = sys.argv[pythonidx + 1:]
    
    from argparse import ArgumentParser, RawTextHelpFormatter

    parser = ArgumentParser("Input Scaling", formatter_class=RawTextHelpFormatter)
    parser.add_argument("action", choices=["setup", "simulate", "analyse"], help="Action to run.")
    parser.add_argument("networkPath", help="Network path")
    parser.add_argument("--neurons", help="Neurons path")
    parser.add_argument("--inputType", help="Type of external input",
                        choices=["thalamic", "cortical"], default="thalamic")
    parser.add_argument("--numInputSteps", type=int, help="Number of steps for number of inputs to neurons",
                        default=10)
    parser.add_argument("--numInputMin", type=int, help="Minimum number of synaptic inputs of inputType", default=100)
    parser.add_argument("--numInputMax", type=int, help="Maximum number of synaptic inputs of inputType", default=1000)
    parser.add_argument("--inputDuration", type=float, default=10.0,
                        help="Duration of each frequency test, longer need for irregularly firing neurons")
    parser.add_argument("--inputFrequency", type=str, default="[0,1,2,5]",
                        help="Input frequency, float or list of floats")
    parser.add_argument("--neuronType", default=None, type=str,
                        help="Optional, if only we want to simulate one neuron type, eg. FS")

    args = parser.parse_args()

    # TODO: Let the user choose input type, duration for each "run", frequency range, number of input range

    input_scaling = InputTuning(args.networkPath)

    if args.action == "setup":
        input_frequency = ast.literal_eval(args.inputFrequency)
        if type(input_frequency) != list:
            input_frequency = np.array(list(input_frequency))

        input_scaling.setup_network(neurons_path=args.neurons,
                                    num_replicas=args.numInputSteps,
                                    neuron_types=args.neuronType)
        input_scaling.setup_input(input_type=args.inputType,
                                  num_input_min=args.numInputMin,
                                  num_input_max=args.numInputMax,
                                  input_duration=args.inputDuration,
                                  input_frequency_range=input_frequency)

        print("Tip, to run in parallel on your local machine use: "
              "mpiexec -n 4 python3 tuning/input_tuning.py simulate <yournetworkhere>")

    elif args.action == "simulate":
        print("Run simulation...")
        print("Tip, to run in parallel on your local machine use: "
              "mpiexec -n 4 python3 tuning/input_tuning.py simulate <yournetworkhere>")
        input_scaling.simulate()

    elif args.action == "analyse":
        # input_scaling.plot_generated_input()
        input_scaling.analyse_results()


    else:
        print(f"Unknown action {args.action}")


    # python3 input_tuning/input_tuning.py setup networks/input_scaling_v1/ data/neurons/striatum/
    # mpiexec -n 4 python3 input_tuning/input_tuning.py simulate networks/input_scaling_v1/ < input.txt &> output-tuning.txt &

    #
