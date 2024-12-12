import ast
import collections
import glob
import json
import os
import sys
import timeit

# Now locally importing matplotlib.pyplot in functions, since Dardel (parallel computer) could not handle it

# Must be run before NEURON import to run in parallel
from mpi4py import MPI

import h5py
import numpy as np
import copy

from snudda.core import Snudda
from snudda.init.init import SnuddaInit
from snudda.input.input import SnuddaInput
from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.place.create_cube_mesh import create_cube_mesh
from snudda.simulate.simulate import SnuddaSimulate
from snudda.utils import SnuddaLoadSimulation
from snudda.utils.load import SnuddaLoad
from snudda.utils.snudda_path import snudda_isdir, snudda_parse_path, snudda_simplify_path, get_snudda_data


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

    def __init__(self, network_path, snudda_data=None, rc=None, input_seed_list=None):

        self.network_path = network_path
        self.neurons_path = None
        self.snudda_data = get_snudda_data(snudda_data=snudda_data)

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
        self.num_replicas = None

        self.rc = rc

        if not os.path.isdir(self.network_path):
            os.makedirs(self.network_path)

        self.network_config_file_name = os.path.join(self.network_path, "network-config.json")
        self.network_file = os.path.join(self.network_path, "network-synapses.hdf5")
        self.input_config_file = os.path.join(self.network_path, "input_config.json")

        self.input_seed_list = input_seed_list

        if self.input_seed_list is None:
            self.input_spikes_file = os.path.join(self.network_path, 'input.hdf5')
            self.output_file = os.path.join(self.network_path, "simulations", "output.hdf5")
        else:
            self.input_spikes_file = [os.path.join(self.network_path, f"input-{s}.hdf5") for s in self.input_seed_list]
            self.output_file = [os.path.join(self.network_path, "simulations", f"output-{s}.hdf5") for s in self.input_seed_list]

        self.core = Snudda(self.network_path)
        self.init_helper = SnuddaInit(network_path=self.network_path, snudda_data=self.snudda_data)

    # Writes config files

    def setup_network(self, neurons_path=None, num_replicas=10, neuron_types=None,
                      parameter_key=None, morphology_key=None, modulation_key=None,
                      reaction_diffusion_file=None,
                      single_neuron_path=None, network_random_seed=None):

        if not morphology_key and not parameter_key and not modulation_key:
            all_combinations = True
        else:
            all_combinations = False

        self.num_replicas = num_replicas

        # TODO: num_replicas should be set by a parameter, it affects how many duplicates of each neuron
        # and thus how many steps we have between n_min and n_max number of inputs specified.
        config_def = self.create_network_config(neurons_path=neurons_path,
                                                snudda_data=self.snudda_data,
                                                num_replicas=num_replicas,
                                                neuron_types=neuron_types,
                                                reaction_diffusion_file=reaction_diffusion_file,
                                                single_neuron_path=single_neuron_path,
                                                parameter_key=parameter_key,
                                                morphology_key=morphology_key,
                                                modulation_key=modulation_key,
                                                all_combinations=all_combinations,
                                                random_seed=network_random_seed)

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

        sd = SnuddaDetect(network_path=self.network_path, rc=self.rc)
        sd.detect()

        sp = SnuddaPrune(network_path=self.network_path, rc=self.rc)
        sp.prune()

        # TODO: Skip placing neurons that will not receive any inputs or distribute any inputs

    def setup_input(self, input_type=None,
                    num_input_min=100, num_input_max=1000, num_input_steps=None,
                    input_duration=10,
                    input_frequency_range=None,
                    input_correlation=None,
                    use_meta_input=True, generate=True, clear_old_input=True):

        if clear_old_input:
            self.input_info = None

        if not input_frequency_range:
            input_frequency_range = [1.0]

        if not num_input_steps:
            num_input_steps = self.num_replicas

        self.frequency_range = np.array(input_frequency_range)
        self.input_duration = input_duration
        self.max_time = self.input_duration * len(self.frequency_range)

        synapse_density_cortical_input = "1.15*0.05/(1+exp(-(d-30e-6)/5e-6))"
        synapse_density_thalamic_input = "0.05*exp(-d/200e-6)"
        #  synapse_density_thalamic_input = "(d > 100e-6)*1"  # TEST!!

        # TODO: These should be read from JSON file, so user can add additional neuron and input types
        cortical_SPN_synapse_parameter_file = "$DATA/synapses/striatum/M1RH_Analysis_190925.h5-parameters-MS.json"
        cortical_contralateral_SPN_synapse_parameter_file = "$DATA/synapses/striatum/M1LH_Analysis_191001.h5-parameters-MS.json"
        thalamic_SPN_synapse_parameter_file = "$DATA/synapses/striatum/TH_Analysis_191001.h5-parameters-MS.json"

        cortical_FS_synapse_parameter_file = "$DATA/synapses/striatum/M1RH_Analysis_190925.h5-parameters-FS.json"
        cortical_contralateral_FS_synapse_parameter_file = "$DATA/synapses/striatum/M1LH_Analysis_191001.h5-parameters-FS.json"
        thalamic_FS_synapse_parameter_file = "$DATA/synapses/striatum/TH_Analysis_191001.h5-parameters-FS.json"

        cortical_ChIN_synapse_parameter_file = "$DATA/synapses/striatum/M1RH_Analysis_190925.h5-parameters-CHAT.json"
        thalamic_ChIN_synapse_parameter_file = "$DATA/synapses/striatum/TH_Analysis_191001.h5-parameters-CHAT.json"

        cortical_LTS_synapse_parameter_file = "$DATA/synapses/striatum/M1RH_Analysis_190925.h5-parameters-LTS.json"

        if "cortical_contralateral" in input_type.lower():
            synapse_density = synapse_density_cortical_input
            synapse_parameter_file = {"dspn": cortical_contralateral_SPN_synapse_parameter_file,
                                      "ispn": cortical_contralateral_SPN_synapse_parameter_file,
                                      "fs": cortical_contralateral_FS_synapse_parameter_file}

        elif 'cortical' in input_type.lower():
            synapse_density = synapse_density_cortical_input
            synapse_parameter_file = {"dspn": cortical_SPN_synapse_parameter_file,
                                      "ispn": cortical_SPN_synapse_parameter_file,
                                      "fs": cortical_FS_synapse_parameter_file,
                                      "lts": cortical_LTS_synapse_parameter_file,
                                      "chin": cortical_ChIN_synapse_parameter_file}
            print("Using cortical synapse density for input.")
        elif 'thalamic' in input_type.lower():
            synapse_density = synapse_density_thalamic_input
            synapse_parameter_file = {"dspn": thalamic_SPN_synapse_parameter_file,
                                      "ispn": thalamic_SPN_synapse_parameter_file,
                                      "fs": thalamic_FS_synapse_parameter_file,
                                      "chin": thalamic_ChIN_synapse_parameter_file}
            print("Using thalamic synapse density for input")
        else:
            synapse_density = "1"
            synapse_parameter_file = {}
            print("No density profile used for input    .")

        if "cluster" in input_type.lower():
            cluster_size = 4
            cluster_spread = 20e-6
        else:
            cluster_size = None
            cluster_spread = None

        self.create_input_config(input_config_file=self.input_config_file,
                                 input_type=input_type,
                                 input_frequency=list(self.frequency_range),  # [1.0],
                                 input_correlation=input_correlation,
                                 n_input_min=num_input_min,
                                 n_input_max=num_input_max,
                                 num_input_steps=num_input_steps,
                                 synapse_conductance=0.5e-9,
                                 synapse_density=synapse_density,
                                 input_duration=self.input_duration,
                                 synapse_parameter_file=synapse_parameter_file,
                                 cluster_size=cluster_size,
                                 cluster_spread=cluster_spread)

        if generate:
            self.generate_input_helper(use_meta_input=use_meta_input)

        # Info we need to run right duration of simulation
        self.write_tuning_info()

    def generate_input_helper(self, use_meta_input=True):
        if self.input_seed_list is None:
            si = SnuddaInput(input_config_file=self.input_config_file,
                             hdf5_network_file=os.path.join(self.network_path, 'network-synapses.hdf5'),
                             spike_data_filename=self.input_spikes_file,
                             time=self.max_time,
                             logfile=os.path.join(self.network_path, "log", "input.txt"),
                             use_meta_input=use_meta_input, rc=self.rc)
            si.generate()
        else:
            for ctr, (seed, input_file) in enumerate(zip(self.input_seed_list, self.input_spikes_file)):
                print(f"Iteration {ctr}/{len(self.input_seed_list)} (seed: {seed})")

                si = SnuddaInput(input_config_file=self.input_config_file,
                                 hdf5_network_file=os.path.join(self.network_path, 'network-synapses.hdf5'),
                                 spike_data_filename=input_file,
                                 time=self.max_time,
                                 logfile=os.path.join(self.network_path, "log", f"input-{seed}.txt"),
                                 use_meta_input=use_meta_input, rc=self.rc, random_seed=seed)
                si.generate()

    def setup_background_input(self, input_types=["cortical_background", "thalamic_background"],
                               input_density=["1.15*0.05/(1+exp(-(d-30e-6)/5e-6))", "0.05*exp(-d/200e-6)"],
                               input_fraction=[0.5, 0.5],
                               num_input_min=10, num_input_max=500,
                               input_frequency=[1, 1], input_duration=10,
                               generate_input=True):

        """ Tries to find the maximum number of synapses that will not make the neuron spike. """

        if np.sum(input_fraction) != 1:
            # If [0.5, 0.5] and 100 inputs, means 50 and 50 synapses for the two inputs
            input_fraction = np.array(input_fraction) / np.sum(input_fraction)
            print(f"Adjusting input fraction to make total 1: {input_fraction}")

        num_input_steps = self.num_replicas
        self.input_duration = input_duration
        self.max_time = self.input_duration   # Only one frequency for background, we are varying number of inputs

        # TODO: Setup network, run neurons at one frequency.

        # Make sure input is cleared on first iteration, and generated on last iteration
        generate_flag = np.zeros((len(input_frequency),), dtype=bool)
        generate_flag[-1] = generate_input
        clear_flag = np.zeros((len(input_frequency),), dtype=bool)
        clear_flag[0] = True

        for input_type, density, fraction, frequency, gf, cf \
            in zip(input_types, input_density, input_fraction, input_frequency,
                   generate_flag, clear_flag):

            self.setup_input(input_type=input_type,
                             num_input_min=np.round(fraction*num_input_min),
                             num_input_max=np.round(fraction*num_input_max),
                             input_duration=input_duration,
                             input_frequency_range=[frequency],
                             use_meta_input=False,
                             generate=gf,
                             clear_old_input=cf)

    def setup_input_verification(self, input_type="cortical", neuron_type="dSPN",
                                 input_frequency_range=None,
                                 input_duration=10,
                                 generate=True, seed_list=None):

        # Here we use the META input, and replace the frequency, start and end variables
        # input_type can be "cortical", "thalamic", "cortical_background", "thalamic_background"

        if not input_frequency_range:
            input_frequency_range = [1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0]

        self.frequency_range = np.array(input_frequency_range)
        self.input_duration = input_duration
        self.max_time = self.input_duration * len(self.frequency_range)

        n_steps = len(self.frequency_range)
        input_start = input_duration * np.arange(0, n_steps)
        input_end = input_duration * np.arange(1, n_steps + 1)

        input_target = neuron_type
        self.input_info = collections.OrderedDict()
        self.input_info[input_target] = collections.OrderedDict()

        self.input_info[input_target][input_type] = collections.OrderedDict()

        self.input_info[input_target][input_type]["generator"] = "poisson"
        self.input_info[input_target][input_type]["start"] = input_start
        self.input_info[input_target][input_type]["end"] = input_end
        self.input_info[input_target][input_type]["frequency"] = self.frequency_range

        with open(self.input_config_file, "w") as f:
            json.dump(self.input_info, f, indent=4, cls=NumpyEncoder)

        if generate:
            self.generate_input_helper(use_meta_input=True)

        # Info we need to run right duration of simulation
        self.write_tuning_info()

    def analyse_results(self, input_type='', show_plots=False):

        frequency_data, voltage_data = self.load_data()

        self.plot_frequency_data(frequency_data, show_plots=show_plots, input_type_name=input_type)
        self.plot_frequency_data_alt(frequency_data, show_plots=show_plots, input_type_name=input_type)
        self.plot_volt_data(voltage_data, show_plots=show_plots, input_type_name=input_type)
        self.plot_volt_vs_ninputs(voltage_data, show_plots=show_plots, input_type_name=input_type)

        print(f"To plot traces:\n"
              f"python3 plotting/Network_plot_traces.py {self.network_path}output_volt.txt "
              f"{self.network_path}network-synapses.hdf5 ")

    def find_signal_strength(self, requested_frequency=10.0, skip_time=0.0, show_plot=True, quiet_load=False):

        spike_data = dict()
        depolarisation_blocks = dict()

        for idx in range(len(self.input_seed_list)):
            network_info, input_config, _, neuron_id_lookup, neuron_name_list, \
                spike_data[idx], _, time, depolarisation_blocks[idx] = \
                self.load_data_helper(idx=idx, load_input=False, quiet_load=quiet_load)

        input_config_info = dict()

        duration = np.max(time) - skip_time
        requested_spikes = requested_frequency * duration

        bad_neuron_list = []

        for neuron_name in neuron_id_lookup.keys():
            neuron_id = neuron_id_lookup[neuron_name]

            spike_count = np.zeros((len(neuron_id), len(self.input_seed_list)), dtype=int)
            depol_block_flag = np.zeros((len(neuron_id), len(self.input_seed_list)), dtype=bool)

            for idx in range(len(self.input_seed_list)):
                spike_count[:, idx], depol_block_flag[:, idx] = \
                    self.extract_background_spikes(spike_data=spike_data[idx],
                                                   neuron_id=neuron_id,
                                                   skip_time=skip_time,
                                                   depolarisation_blocks=depolarisation_blocks[idx])

            spike_count_mean = np.mean(spike_count, axis=1)

            best_config, neuron_info, depol_block = self.get_best_config(data=spike_count_mean, requested_value=requested_spikes,
                                                                         neuron_id=neuron_id,
                                                                         input_config=input_config,
                                                                         network_info=network_info,
                                                                         depol_block_flag=depol_block_flag)

            if depol_block:
                bad_neuron_list.append(neuron_name)

            input_config_info[neuron_name] = (best_config,
                                              snudda_parse_path(os.path.join(neuron_info["neuron_path"], "meta.json"),
                                                                snudda_data=self.snudda_data),
                                              neuron_info["parameter_key"],
                                              neuron_info["morphology_key"])

            if neuron_info["parameter_key"] not in neuron_name \
                    or neuron_info["morphology_key"] not in neuron_name:
                print(f"{neuron_name = }, {neuron_info['parameter_key'] = }, {neuron_info['morphology_key'] = }")
                print(f"Did you accidentally use different networks for the runs?")
                import pdb
                pdb.set_trace()

            # Just an idiot check to make sure all neurons we are comparing are the same
            for nid in neuron_id:
                assert network_info.data["neurons"][neuron_id[0]]["name"] == network_info.data["neurons"][nid]["name"]

            # Check if spike frequency is too low, if so give image a different name

            move_bad = False

            if depol_block:
                label = f"signal-{requested_frequency}-Hz-BLOCKED"
                move_bad = True
            else:
                label = f"signal-{requested_frequency}-Hz"

            if np.max(spike_count_mean) < requested_spikes:
                label = f"{label}-TOO-LOW-FREQ"
                move_bad = True

            self.plot_signal_info(neuron_id=neuron_id, neuron_info=neuron_info, best_config=best_config,
                                  spike_count=spike_count, input_config=input_config, max_time=np.max(time),
                                  requested_frequency=requested_frequency, depol_block_flag=depol_block_flag,
                                  label=label, show_plot=show_plot, move_bad=move_bad)

        for name in bad_neuron_list:
            print(f"Found early depolarisation block: {name}")

        return input_config_info

    def get_best_config(self, data, requested_value, neuron_id, input_config, network_info,
                        depol_block_flag):

        all_morph_keys = np.array([network_info.data["neurons"][nid]["morphology_key"] for nid in neuron_id])
        all_param_keys = np.array([network_info.data["neurons"][nid]["parameter_key"] for nid in neuron_id])

        assert (all_param_keys == all_param_keys[0]).all(), \
            f"Internal error: All parameter_key should be the same {all_param_keys}"
        assert (all_morph_keys == all_morph_keys[0]).all(), \
            f"Internal error: All morphology_keys should be the same {all_morph_keys}"

        idx_above = np.argmax(data > requested_value)

        if data[idx_above] < requested_value:
            # If we did not reach the requested value, pick the closest point.
            idx_above = idx_below = np.argmin(abs(data - requested_value))

        elif idx_above == 0:
            idx_below = 0
        else:
            idx_below = idx_above - 1

        config_above = input_config[str(neuron_id[idx_above])]
        config_below = input_config[str(neuron_id[idx_below])]

        input_type = list(config_above.keys())
        assert len(input_type) == 1, f"Interpolation can only handle one input type. Found {input_type}"
        input_type = input_type[0]

        n_syn_above = config_above[input_type]["num_inputs"]
        n_syn_below = config_below[input_type]["num_inputs"]

        value_above = data[idx_above]
        value_below = data[idx_below]

        depol_block = depol_block_flag[:idx_above, :].any()

        assert value_below <= value_above

        if n_syn_above == n_syn_below:
            n_syn = n_syn_above
        else:
            n_syn = int(np.round(n_syn_below
                                 + (n_syn_above - n_syn_below) * (requested_value - value_below)
                                 / (value_above - value_below)))

        try:
            assert n_syn_below <= n_syn <= n_syn_above, f"NOT TRUE: n_syn_below {n_syn_below} <= n_syn {n_syn} <= n_syn_above {n_syn_above}"
        except:
            import traceback
            import pdb
            print(traceback.format_exc())
            pdb.set_trace()

        best_config = copy.deepcopy(config_above)
        best_config[input_type]["num_inputs"] = n_syn

        neuron_info = network_info.data["neurons"][neuron_id[idx_above]]

        return best_config, neuron_info, depol_block

    def find_highest_non_spiking_background_input(self, skip_time=0.0, show_plot=True, quiet_load=False):

        spike_data = dict()
        depolarisation_blocks = dict()

        for idx in range(len(self.input_seed_list)):
            network_info, input_config, _, neuron_id_lookup, neuron_name_list, \
                spike_data[idx], _, time, depolarisation_blocks[idx] = \
                self.load_data_helper(idx=idx, load_input=False, quiet_load=quiet_load)

        input_config_info = dict()

        for neuron_name in neuron_id_lookup.keys():
            neuron_id = neuron_id_lookup[neuron_name]

            spike_count = np.zeros((len(neuron_id), len(self.input_seed_list)), dtype=int)
            depol_block_flag = np.zeros((len(neuron_id), len(self.input_seed_list)), dtype=bool)

            for idx in range(len(self.input_seed_list)):
                spike_count[:, idx], depol_block_flag[:, idx] = \
                    self.extract_background_spikes(spike_data=spike_data[idx],
                                                   neuron_id=neuron_id,
                                                   skip_time=skip_time,
                                                   depolarisation_blocks=depolarisation_blocks[idx])

            spike_count_sum = np.sum(spike_count, axis=1)

            # TODO: Need to handle if there are spikes in all traces
            best_ctr = 0
            for ctr, sc in enumerate(spike_count_sum):
                if sc > 0:
                    break
                else:
                    best_ctr = ctr

            best_neuron_id = neuron_id[best_ctr]

            neuron_info = network_info.data["neurons"][best_neuron_id]
            input_config_info[neuron_name] = (input_config[str(best_neuron_id)],
                                              snudda_parse_path(os.path.join(neuron_info["neuron_path"], "meta.json"),
                                                                snudda_data=self.snudda_data),
                                              neuron_info["parameter_key"],
                                              neuron_info["morphology_key"])

            # Just an idiot check to make sure all neurons we are comparing are the same
            for nid in neuron_id:
                assert network_info.data["neurons"][neuron_id[0]]["name"] == network_info.data["neurons"][nid]["name"]

            self.plot_background_info(neuron_id=neuron_id, neuron_info=neuron_info, best_neuron_id=best_neuron_id,
                                      spike_count=spike_count, input_config=input_config, depol_block_flag=depol_block_flag,
                                      max_time=np.max(time), label="background-inputs", show_plot=show_plot)

        return input_config_info

    def plot_background_info(self, neuron_id, neuron_info, best_neuron_id, spike_count, input_config,
                             max_time, skip_time=0, label="background-inputs", show_plot=True,
                             depol_block_flag=None):

        import matplotlib.pyplot as plt

        n_inputs_total = np.zeros((len(neuron_id),), dtype=int)
        fig_dir = os.path.join(self.network_path, "figures")

        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)

        fig_name = os.path.join(fig_dir, f"{neuron_info['morphology_key']}-{neuron_info['parameter_key']}-{neuron_info['name']}-{label}.png")

        for ctr, nid in enumerate(neuron_id):
            # Get total input.
            for input_conf in input_config[str(nid)].values():
                n_inputs_total[ctr] += input_conf["num_inputs"]

        best_idx = np.where(neuron_id == best_neuron_id)[0]

        plt.figure()
        plt.plot(n_inputs_total, spike_count/(max_time-skip_time), 'k.')

        if depol_block_flag is not None:
            # Mark the depolarisation blocks with a red circle.
            bad_idx = np.where(depol_block_flag)
            plt.plot(n_inputs_total[bad_idx[0]], spike_count[bad_idx]/(max_time-skip_time), 'ro')

        plt.plot(n_inputs_total[best_idx], spike_count[best_idx]/(max_time-skip_time), 'b*')
        plt.xlabel("Total number of synapses")
        plt.ylabel("Spike frequency")
        plt.title(f"Neuron {neuron_info['name']}")

        print(f"Saving figure to {fig_name}")
        plt.savefig(fig_name)

        if show_plot:
            plt.ion()
            plt.show()

    def plot_signal_info(self, neuron_id, neuron_info, best_config, spike_count, input_config,
                         max_time, requested_frequency, skip_time=0,
                         label="background-inputs", show_plot=True,
                         depol_block_flag=None, move_bad=False):

        import matplotlib.pyplot as plt

        n_inputs_total = np.zeros((len(neuron_id),), dtype=int)
        fig_dir = os.path.join(self.network_path, "figures")

        if move_bad:
            fig_dir = os.path.join(fig_dir, "_BAD")

        if not os.path.isdir(fig_dir):
            os.makedirs(fig_dir)

        fig_name = os.path.join(fig_dir, f"{neuron_info['morphology_key']}-{neuron_info['parameter_key']}-{neuron_info['name']}-{label}.png")

        for ctr, nid in enumerate(neuron_id):
            # Get total input.
            for input_conf in input_config[str(nid)].values():
                n_inputs_total[ctr] += input_conf["num_inputs"]

        plt.figure()
        plt.plot(n_inputs_total, spike_count/(max_time-skip_time), 'k.')

        if depol_block_flag is not None:
            # Mark the depolarisation blocks with a red circle.
            bad_idx = np.where(depol_block_flag)
            plt.plot(n_inputs_total[bad_idx[0]], spike_count[bad_idx]/(max_time-skip_time), 'ro')

        input_type = list(best_config.keys())
        assert len(input_type) == 1
        input_type = input_type[0]
        best_n_syn = best_config[input_type]["num_inputs"]

        y_lim = plt.gca().get_ylim()
        plt.plot([best_n_syn, best_n_syn], y_lim, 'r-')  # best n_inputs

        x_lim = plt.gca().get_xlim()
        plt.plot(x_lim, [requested_frequency, requested_frequency], 'k--')

        # plt.plot(n_inputs_total[best_idx], spike_count[best_idx]/(max_time-skip_time), 'r*')
        plt.xlabel("Total number of synapses")
        plt.ylabel("Spike frequency")
        plt.title(f"Neuron {neuron_info['name']} ({best_n_syn} {input_type} synapses)")

        print(f"Saving figure to {fig_name}")
        plt.savefig(fig_name)

        if show_plot:
            plt.ion()
            plt.show()

    def update_meta(self, input_config_info, overwrite=True, set_frequency=None):

        last_meta_file = None
        meta_data = None

        for input_config, meta_file, parameter_key, morphology_key in input_config_info.values():

            if meta_file != last_meta_file:
                # Write old data to file
                if last_meta_file is not None:
                    print(f"Writing {last_meta_file}")
                    with open(last_meta_file, "w") as f:
                        json.dump(meta_data, f, indent=4)

                # Load new data
                with open(meta_file, "r") as f:
                    meta_data = json.load(f)
                last_meta_file = meta_file

            new_config = copy.deepcopy(input_config)
            try:
                for input_name in new_config.keys():
                    if set_frequency is not None:
                        # Override the frequency (this is useful to have a signal, but not have it active)
                        new_config[input_name]["frequency"] = set_frequency
                    else:
                        # Use the frequency from the input tuning
                        new_config[input_name]["frequency"] = new_config[input_name]["frequency"][0]

                    del new_config[input_name]["start"]
                    del new_config[input_name]["end"]

                    # If parameter_file and parameter-ist both are given, only keep the latter
                    if "parameter_file" in new_config[input_name] and "parameter_list" in new_config[input_name]:
                        del new_config[input_name]["parameter_list"]
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

            if parameter_key not in meta_data or morphology_key not in meta_data[parameter_key]:
                print(f"Parameter key {parameter_key}, morphology key {morphology_key} not found in {meta_file} -- was it manually removed?")
                continue

            # print(f"Writing {parameter_key = }, {morphology_key = }")
            if overwrite:
                meta_data[parameter_key][morphology_key]["input"] = new_config
            else:
                meta_data[parameter_key][morphology_key]["input"] |= new_config

        # Write last iteration to file
        if last_meta_file is not None:
            print(f"Writing {last_meta_file}")
            with open(last_meta_file, "w") as f:
                json.dump(meta_data, f, indent=4)

    def load_data(self, skip_time=0.0, quiet_load=False):

        input_data = dict()
        spike_data = dict()
        volt = dict()
        time = dict()
        depolarisation_blocks = dict()

        if self.input_seed_list is None:
            network_info, input_config, input_data, neuron_id_lookup, neuron_name_list, \
            spike_data, volt, time, depolarisation_blocks \
                = self.load_data_helper(idx=None, quiet_load=quiet_load)

            first_input_data = input_data["input"]
        else:
            for idx in range(len(self.input_seed_list)):
                network_info, input_config, input_data[idx], neuron_id_lookup, neuron_name_list, \
                    spike_data[idx], volt[idx], time[idx], depolarisation_blocks[idx] \
                    = self.load_data_helper(idx=idx, quiet_load=quiet_load)

            first_input_data = input_data[0]["input"]

        n_inputs_lookup = dict()

        for neuron_label in first_input_data:
            neuron_id = int(neuron_label)
            n_inputs = 0

            for input_type in first_input_data[neuron_label]:
                n_inputs += first_input_data[neuron_label][input_type]["spikes"].shape[0]

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

                if self.input_seed_list is None:

                    n_inputs = n_inputs_lookup[neuron_id]

                    if n_inputs not in frequency_data[neuron_name]:
                        frequency_data[neuron_name][n_inputs] = dict()

                    frequency_data[neuron_name][n_inputs] = self.extract_frequencies(spike_data=spike_data,
                                                                                     config_data=input_config,
                                                                                     neuron_id=neuron_id,
                                                                                     skip_time=skip_time)
                    voltage_data[neuron_name][n_inputs] = self.extract_voltage(volt=volt,
                                                                               time=time,
                                                                               config_data=input_config,
                                                                               neuron_id=neuron_id,
                                                                               skip_time=skip_time)

                else:

                    for idx in range(len(self.input_seed_list)):

                        n_inputs = n_inputs_lookup[neuron_id]

                        if n_inputs not in frequency_data[neuron_name]:
                            frequency_data[neuron_name][n_inputs] = dict()

                        frequency_data[neuron_name][n_inputs][idx] = self.extract_frequencies(spike_data=spike_data,
                                                                                              config_data=input_config,
                                                                                              neuron_id=neuron_id,
                                                                                              skip_time=skip_time)
                        voltage_data[neuron_name][n_inputs][idx] = self.extract_voltage(volt=volt,
                                                                                        time=time,
                                                                                        config_data=input_config,
                                                                                        neuron_id=neuron_id,
                                                                                        skip_time=skip_time)

        # TODO: Load voltage trace and warn for depolarisation blocking

        return frequency_data, voltage_data

    def load_data_helper(self, idx=None, load_input=True, quiet_load=False):

        network_file = os.path.join(self.network_path, "network-synapses.hdf5")
        network_info = SnuddaLoad(network_file)
        input_data = None

        if idx is None:
            output_file = self.output_file
            if load_input and os.path.isfile(self.input_spikes_file):
                input_data = h5py.File(self.input_spikes_file, "r")
        else:
            output_file = self.output_file[idx]
            if load_input and os.path.isfile(self.input_spikes_file[idx]):
                input_data = h5py.File(self.input_spikes_file[idx], "r")

        output_data_loader = SnuddaLoadSimulation(network_path=self.network_path,
                                                  network_simulation_output_file=output_file,
                                                  do_test=True, quiet_load=quiet_load)
        spike_data = output_data_loader.get_spikes()

        # cell_id = output_data_loader.get_id_of_neuron_type()
        volt = output_data_loader.get_voltage()
        time = output_data_loader.get_time()
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

        depolarisation_blocks = output_data_loader.get_depolarisation_dictionary()

        return network_info, input_config, input_data, neuron_id_lookup, neuron_name_list, spike_data, volt, time, depolarisation_blocks

    def is_depolarisation_blocked(self, neuron_id, time_range, depolarisation_blocks):

        # Checks if there are any depolarisation blocks for neuron_id within the time_range

        is_blocked = False

        if neuron_id in depolarisation_blocks:

            for start_block, end_block in depolarisation_blocks[neuron_id]:

                if time_range[0] < end_block and time_range[1] > start_block:
                    is_blocked = True

        return is_blocked

    def plot_depolarisation_blocked_neurons(self, freq_bin=10):

        import matplotlib.pyplot as plt

        spike_data = dict()
        volt = dict()
        depolarisation_blocks = dict()

        if type(self.input_spikes_file) == list:
            idx_list = np.arange(0, len(self.input_spikes_file))

            for idx in idx_list:
                network_info, input_config, _, neuron_id_lookup, neuron_name_list, \
                    spike_data[idx], volt[idx], time, depolarisation_blocks[idx] = \
                    self.load_data_helper(idx=idx, load_input=False, quiet_load=True)

        else:
            network_info, input_config, _, neuron_id_lookup, neuron_name_list, \
                spike_data[0], volt[0], time, depolarisation_blocks[0] = \
                self.load_data_helper(load_input=False, quiet_load=True)

        bad_neurons = []

        for idx in depolarisation_blocks:
            bad_neurons += list(depolarisation_blocks[idx].keys())

        for neuron_id in bad_neurons:
            plt.figure()

            for ctr, idx in enumerate(volt):
                plt.plot(time*1e3, volt[idx][neuron_id]*1e3)
                plt.plot(spike_data[idx][neuron_id]*1e3,
                         np.full(spike_data[idx][neuron_id].shape, 40+ctr), 'k.')
                if neuron_id in depolarisation_blocks[idx]:
                    for bad_start, bad_end in depolarisation_blocks[idx][neuron_id]:
                        plt.plot([bad_start*1e3, bad_end*1e3], [50, 50], 'r')

            full_morph_key = network_info.data["neurons"][neuron_id]["morphology_key"]
            full_param_key = network_info.data["neurons"][neuron_id]["parameter_key"]
            neuron_type = network_info.data["neurons"][neuron_id]["type"]

            if freq_bin is not None:
                n_bins = int(np.floor(np.max(time) / freq_bin)) + 1
                freq_data = np.zeros((n_bins, len(spike_data.keys())))

                for idx in spike_data:
                    try:
                        unique, counts = np.unique(np.floor(spike_data[idx][neuron_id].flatten() / freq_bin).astype(int), return_counts=True)
                    except:
                        import traceback
                        print(traceback.format_exc())
                        import pdb
                        pdb.set_trace()

                    for val, ctr in zip(unique, counts):
                        freq_data[val, idx] = ctr / freq_bin

                freq_vals = np.mean(freq_data, axis=1)
                for idx, fv in enumerate(freq_vals):
                    plt.text((idx+0.5)*freq_bin*1e3, -10, f"{fv:.1f}Hz",
                             horizontalalignment="center")

            plt.title(f"{neuron_type}, param: {full_param_key}, morph: {full_morph_key} (id {neuron_id})")
            plt.xlabel("Time (ms)")
            plt.ylabel("Voltage (mV)")

            fig_path = os.path.join(self.network_path, "figures", "_bad",
                                    f"{full_morph_key}-{full_param_key}-{neuron_type}-BAD-trace.png")

            if not os.path.exists(os.path.dirname(fig_path)):
                os.makedirs(os.path.dirname(fig_path))

            plt.savefig(fig_path, dpi=300)

            plt.ion()
            plt.show()

    def plot_voltage_trace(self, morphology_key, parameter_key, mp_idx=None, time_range=None):

        import matplotlib.pyplot as plt

        # Find all neurons with the morphology key, and parameter key

        # If idx is given pick the n:th trace specified to plot, if not plot all traces.
        # idx can be a list of multiple traces.

        spike_data = dict()
        volt = dict()
        depolarisation_blocks = dict()

        if type(self.input_spikes_file) == list:
            idx_list = np.arange(0, len(self.input_spikes_file))

            for idx in idx_list:
                network_info, input_config, _, neuron_id_lookup, neuron_name_list, \
                    spike_data[idx], volt[idx], time, depolarisation_blocks[idx] = \
                    self.load_data_helper(idx=idx, load_input=False, quiet_load=True)

        else:
            network_info, input_config, _, neuron_id_lookup, neuron_name_list, \
                spike_data[0], volt[0], time, depolarisation_blocks[0] = \
                self.load_data_helper(load_input=False, quiet_load=True)

        morph_flag = np.array([morphology_key in n["morphology_key"] for n in network_info.data["neurons"]], dtype=bool)
        param_flag = np.array([parameter_key in n["parameter_key"] for n in network_info.data["neurons"]], dtype=bool)

        match_idx = np.where(np.logical_and(morph_flag, param_flag))[0]

        if mp_idx is None:
            plot_idx = match_idx
        else:
            plot_idx = match_idx[mp_idx]

        full_morph_key = network_info.data["neurons"][plot_idx[0]]["morphology_key"]
        full_param_key = network_info.data["neurons"][plot_idx[0]]["parameter_key"]
        neuron_type = network_info.data["neurons"][plot_idx[0]]["type"]

        fig = plt.figure()
        for r_idx in volt.keys():
            for ctr, idx in enumerate(plot_idx):
                spike_times = spike_data[r_idx][idx]

                if time_range is None:
                    plt.plot(time*1e3, volt[r_idx][idx]*1e3)
                    plt.plot(spike_times*1e3, np.full(spike_times.shape, 40+ctr), 'k.')
                else:
                    t_idx = np.logical_and(time_range[0] <= time, time <= time_range[1])
                    plt.plot(time[t_idx]*1e3, volt[r_idx][idx][t_idx]*1e3)
                    s_idx = np.where(np.logical_and(time_range[0] <= spike_times,
                                                    spike_times <= time_range[1]))
                    plt.plot(spike_times[s_idx]*1e3, np.full(spike_times[s_idx].shape, 40+ctr), 'k.')

                if idx in depolarisation_blocks[r_idx]:
                    bad_ranges = depolarisation_blocks[r_idx][idx]
                    for bad_start, bad_end in bad_ranges:
                        if time_range is None or (bad_start <= time_range[1] and bad_end >= time_range[0]):
                            plt.plot([bad_start*1e3, bad_end*1e3], [50, 50], 'r')

        plt.title(f"{neuron_type}, param: {full_param_key}, morpg: {full_morph_key}")
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")

        fig_path = os.path.join(self.network_path, "figures", f"Trace-{neuron_type}-{full_param_key}-{full_morph_key}.png")
        plt.savefig(fig_path)

        plt.ion()
        plt.show()

    # TODO: We should set skip_time to 1 second, 0 for now while testing
    def extract_frequencies(self, spike_data, config_data, neuron_id, skip_time=0.0):

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

    def extract_background_spikes(self, spike_data, neuron_id, skip_time=0.0, depolarisation_blocks=None):

        spike_count = np.zeros((len(neuron_id), ), dtype=int)

        for ctr, n_id in enumerate(neuron_id):
            spike_count[ctr] = len(np.where(spike_data[n_id] >= skip_time)[0])

        if depolarisation_blocks is not None:
            depol_block_flag = np.zeros((len(neuron_id),), dtype=bool)
            for ctr, n_id in enumerate(neuron_id):
                if n_id in depolarisation_blocks:
                    depol_block_flag[ctr] = True

            return spike_count, depol_block_flag

        return spike_count

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

    def plot_volt_data(self, volt_data, show_plots=True, input_type_name=''):

        import matplotlib.pyplot as plt

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
                ax.plot(input_freq, mean_voltage, linestyle='-', marker='o')

            plt.title(f"{neuron_name} receiving {input_type_all} input")
            plt.xlabel("Input frequency (per synapse)")
            plt.ylabel("Mean membrane potential (V)")
            ax.legend(legend_text)

            if show_plots:
                plt.ion()
                plt.show()
                plt.pause(0.001)

            # fig_name = os.path.join(self.network_path, "figures", f"input-scaling-freq-{neuron_name}-mean-voltage-{input_type_name}.pdf")
            fig_name = os.path.join(self.network_path, "figures",
                                    f"input-scaling-freq-{neuron_name}-mean-voltage-{input_type_name}.png")
            if not os.path.exists(os.path.dirname(fig_name)):
                os.mkdir(os.path.dirname(fig_name))

            plt.savefig(fig_name, dpi=300)

            if not show_plots:
                plt.close()

    def plot_volt_vs_ninputs(self, volt_data, show_plots=True, input_type_name=''):

        import matplotlib.pyplot as plt

        for neuron_name in volt_data:
            fig, ax = plt.subplots()
            legend_text = []
            input_type_all = None

            cmap = plt.get_cmap('tab20', len(volt_data[neuron_name]))
            ax.set_prop_cycle('color', [cmap(i) for i in range(0, len(volt_data[neuron_name]))])
            num_inputs = []
            mean_voltages = []
            for num_input in volt_data[neuron_name]:
                input_freq, mean_voltage, max_voltage, input_type = volt_data[neuron_name][num_input]
                num_inputs.append(num_input)
                mean_voltages.append(mean_voltage)
                if input_type_all:
                    assert input_type == input_type_all, "All input types must be the same for neuron"
                else:
                    input_type_all = input_type

            # legend_text.append(f"input freq.={input_freq[0]}Hz")#Currently assumes one frequency used. Fix later if necessary
            # ax.legend(legend_text)
            ax.plot(num_inputs, mean_voltages, linestyle='None', marker='o')
            # plt.title(f"{neuron_name} receiving {input_type_all} input")
            plt.xlabel("Number of synapses")
            plt.ylabel("Mean membrane potential (V)")

            if show_plots:
                plt.ion()
                plt.show()
                plt.pause(0.001)

            # fig_name = os.path.join(self.network_path, "figures", f"input-scaling-ninputs-{neuron_name}-mean-voltage-{input_type_name}.pdf")
            fig_name = os.path.join(self.network_path, "figures",
                                    f"input-scaling-ninputs-{neuron_name}-mean-voltage-{input_type_name}.png")
            if not os.path.exists(os.path.dirname(fig_name)):
                os.mkdir(os.path.dirname(fig_name))

            plt.savefig(fig_name, bbox_inches='tight', dpi=300)

            if not show_plots:
                plt.close()

    # TODO: Extract spiking frequency (skip first second for each interval to let network settle)
    # TODO: Create summary graphs

    def plot_frequency_data(self, frequency_data, show_plots=True, input_type_name=''):

        import matplotlib.pyplot as plt

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

            fig_name = os.path.join(self.network_path, "figures",
                                    f"input-scaling-freq-{neuron_name}-{input_type_name}.png")
            if not os.path.exists(os.path.dirname(fig_name)):
                os.mkdir(os.path.dirname(fig_name))

            plt.savefig(fig_name, dpi=300)

            if not show_plots:
                plt.close()

    def plot_frequency_data_alt(self, frequency_data, show_plots=True, input_type_name=''):

        import matplotlib.pyplot as plt

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
            plt.ylabel("Spikes/s")
            ax.legend(legend_text)

            if show_plots:
                plt.ion()
                plt.show()
                plt.pause(0.001)

            # fig_name = os.path.join(self.network_path, "figures", f"input-scaling-ninputs-{neuron_name}-{input_type_name}.pdf")
            fig_name = os.path.join(self.network_path, "figures",
                                    f"input-scaling-ninputs-{neuron_name}-{input_type_name}.png")
            if not os.path.exists(os.path.dirname(fig_name)):
                os.mkdir(os.path.dirname(fig_name))

            plt.savefig(fig_name, dpi=300)

            if not show_plots:
                plt.close()

    def plot_verify_frequency_distribution(self, input_type="cortical"):

        import matplotlib.pyplot as plt

        network_info, input_config, input_data, neuron_id_lookup, neuron_name_list, \
            spike_data, volt, time, depolarisation_blocks = self.load_data_helper()

        # First find out what time ranges we need to look at, and what the input frequency is for those
        neuron_type = list(input_config.keys())
        assert len(neuron_type) == 1, f"Plot only supports one neuron type at a time, neuron_type = {neuron_type}"
        neuron_type = neuron_type[0]

        assert input_type in input_config[neuron_type], f"{input_type} not in input_config: {input_config}"
        start_times = np.array(input_config[neuron_type][input_type]["start"])
        end_times = np.array(input_config[neuron_type][input_type]["end"])
        input_freq = np.array(input_config[neuron_type][input_type]["frequency"])
        output_freq_list = []

        depol_blocked_freqs = []
        depol_blocked_lookup = dict()

        for neuron_id in sorted(spike_data.keys()):
            out_freq = []
            for start_t, end_t, in_freq in zip(start_times, end_times, input_freq):
                n_spikes = np.sum(np.logical_and(start_t <= spike_data[neuron_id], spike_data[neuron_id] < end_t))
                f = n_spikes / (end_t - start_t)
                out_freq.append(f)

                if self.is_depolarisation_blocked(neuron_id=neuron_id, time_range=(start_t, end_t),
                                                  depolarisation_blocks=depolarisation_blocks):
                    depol_blocked_freqs.append((in_freq, f))

                    if in_freq not in depol_blocked_lookup:
                        depol_blocked_lookup[in_freq] = [f]
                    else:
                        depol_blocked_lookup[in_freq].append(f)

            output_freq_list.append(out_freq)

        output_freq = np.array(output_freq_list).T

        if len(depol_blocked_freqs) > 0:
            plt.figure()
            plt.plot(input_freq, output_freq, 'k')

            depol_in_freq, depol_out_freq = zip(*depol_blocked_freqs)
            plt.plot(depol_in_freq, depol_out_freq, 'r*')

            plt.xlabel("Input frequency")
            plt.ylabel("Output frequency")
            plt.ion()
            plt.show()

        one_mat = np.ones((2,))

        bad_idx = set()

        fig, ax = plt.subplots(len(input_freq), 1, figsize=(6, 12))

        for idx, (in_freq, out_freq) in enumerate(zip(input_freq, output_freq)):

            # ax[idx].hist(out_freq, bins=50)
            counts, bins = np.histogram(out_freq, bins=50)
            ax[idx].stairs(counts, bins, color="black")

            if in_freq in depol_blocked_lookup:
                depol_counts, depol_bins = np.histogram(depol_blocked_lookup[in_freq], bins=bins)
                ax[idx].stairs(depol_counts, depol_bins, color="red")

            yl = ax[idx].get_ylim()
            mean_freq = np.mean(out_freq)
            std_freq = np.std(out_freq)

            ax[idx].plot(one_mat*mean_freq, yl, 'k-')
            ax[idx].plot(one_mat*(mean_freq + 2*std_freq), yl, 'k--')
            ax[idx].plot(one_mat*(mean_freq - 2*std_freq), yl, 'k--')
            ax[idx].plot([in_freq, in_freq], yl, 'b-')

            bad_idx = bad_idx.union(set(np.where(np.logical_or(out_freq < mean_freq - 2*std_freq,
                                                               out_freq > mean_freq + 2*std_freq))[0]))

            ax[idx].set_xlabel("Spiking frequency (Hz)")

        plt.ion()
        plt.show()

        for idx in bad_idx:
            bad_neuron = network_info.data["neurons"][idx]
            print(f"Frequency outliers: {bad_neuron['name']} ({idx}) -- {output_freq[:, idx]} Hz\n{bad_neuron['neuron_path']}")

    def get_neuron_info(self, neuron_path):

        neuron_path = snudda_parse_path(neuron_path, self.snudda_data)

        neuron_info = collections.OrderedDict()

        # Find neuron morphology swc file, obs currently assume lowercase(!)
        neuron_morph = self.init_helper.get_morphologies(neuron_path)

        parameter_file = os.path.join(neuron_path, "parameters.json")
        mechanism_file = os.path.join(neuron_path, "mechanisms.json")
        modulation_file = os.path.join(neuron_path, "modulation.json")  # Optional
        meta_file = os.path.join(neuron_path, "meta.json")

        # Check if empty neuron_morph_list, or if more than one morphology
        assert os.path.isfile(parameter_file), f"Missing parameter file {parameter_file}"
        assert os.path.isfile(mechanism_file), f"Missing mechanism file {mechanism_file}"

        # TODO: We need to have neuron_info["neuron_path"][name_of_neuron] = neuron_path
        # The below line is wrong, but it is FRIDAY evening so will fix it monday...
        neuron_info["neuron_path"] = snudda_simplify_path(neuron_path, self.snudda_data)

        # OBS, these are not used by snudda when placing neurons, it is just for internal bookkeeping of input tuning
        neuron_info["morphology"] = snudda_simplify_path(neuron_morph, self.snudda_data)
        neuron_info["parameters"] = snudda_simplify_path(parameter_file, self.snudda_data)
        neuron_info["mechanisms"] = snudda_simplify_path(mechanism_file, self.snudda_data)
        neuron_info["meta"] = snudda_simplify_path(meta_file, self.snudda_data)

        # Modulation file is optional
        if os.path.isfile(modulation_file):
            neuron_info["modulation"] = snudda_simplify_path(modulation_file, self.snudda_data)

        return neuron_info

    # This loops through all single neuron directories in neurons_path
    # in preparation of writing a network config file

    def gather_all_neurons(self, neuron_types=None, all_combinations=True):
        all_neurons = collections.OrderedDict()

        assert snudda_isdir(self.neurons_path, self.snudda_data), \
            f"Neurons directory {self.neurons_path} does not exist."

        neuron_type_dir = [d for d in glob.glob(os.path.join(snudda_parse_path(self.neurons_path, self.snudda_data),
                                                             '*'))
                           if snudda_isdir(d, self.snudda_data)]

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
                        param_key = ni["parameter_key"]
                        morph_key = ni["morphology_key"]
                        short_name = n_name[:min(10, len(n_name))]
                        neuron_name = f"{neuron_type.replace('-','_')}_{short_name}_{param_key}_{morph_key}".replace("-", "_")
                        all_neurons[neuron_name] = ni
                        neuron_ctr += 1
                else:
                    neuron_name = os.path.basename(os.path.dirname(neuron_info["parameters"])).replace("-", "_")
                    neuron_ctr += 1

                    all_neurons[neuron_name] = neuron_info

            if neuron_ctr > 0:
                print(f"Found {neuron_ctr} neuron models in {ntd}")

        assert len(all_neurons) > 0, (f"No neurons selected. Did you specify an incorrect neuronType? {neuron_types}"
                                      f"\nSee skipped neurons above error message for available ones.")

        return all_neurons

    def get_all_combinations(self, neuron_info):

        assert "meta" in neuron_info and neuron_info["meta"]

        pm_list = []

        with open(snudda_parse_path(neuron_info["parameters"], self.snudda_data), "r") as pf:
            param_data = json.load(pf)

        with open(snudda_parse_path(neuron_info["meta"], self.snudda_data), "r") as mf:
            meta_data = json.load(mf)

        for p_idx, p_key in enumerate(param_data):
            if p_key not in meta_data:
                print(f"parameter key {p_key} missing in {neuron_info['meta']}")
                # Skip this key
                continue

            for m_idx, m_key in enumerate(meta_data[p_key]):
                ni = copy.deepcopy(neuron_info)
                ni["morphology_key"] = m_key
                ni["parameter_key"] = p_key
                ni["morphology_file"] = os.path.join(ni["morphology"], meta_data[p_key][m_key]["morphology"])
                ni["parameter_id"] = p_idx
                ni["morphology_id"] = m_idx

                pm_list.append(ni)

        return pm_list

    def has_axon(self, neuron_info):

        nm = NeuronPrototype(neuron_name="JJJ",
                             neuron_path=None,
                             snudda_data=self.snudda_data,
                             morphology_path=neuron_info["morphology"],
                             parameter_path=neuron_info["parameters"],
                             mechanism_path=neuron_info["mechanisms"],
                             virtual_neuron=False)
        nm.instantiate()
        return nm.all_have_axon()

    def create_network_config(self,
                              neurons_path=None,
                              snudda_data=None,
                              num_replicas=10,
                              random_seed=None,
                              neuron_types=None,
                              reaction_diffusion_file=None,
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

        if snudda_data is None:
            snudda_data = self.snudda_data

        config_def = dict()
        config_def["snudda_data"] = snudda_data
        config_def["random_seed"], self.init_rng = SnuddaInit.setup_random_seeds(random_seed)

        region_def = dict()
        vol_name = "InputTest"
        region_def[vol_name] = dict()
        region_def[vol_name]["volume"] = dict()

        region_def[vol_name]["volume"]["type"] = "mesh"
        region_def[vol_name]["volume"]["d_min"] = 15e-6
        region_def[vol_name]["volume"]["mesh_file"] = "data/mesh/InputTestMesh.obj"
        region_def[vol_name]["volume"]["num_putative_points"] = 100000
        region_def[vol_name]["connectivity"] = dict()  # Unconnected

        if single_neuron_path:
            # Override and only get one neuron
            assert type(neuron_types) == str, "neuron_types must be string if single_neuron_path set"
            neuron_def = collections.OrderedDict()
            neuron_def[neuron_types] = self.get_neuron_info(single_neuron_path)

            if parameter_key:
                neuron_def[neuron_types]["parameter_key"] = parameter_key

            if morphology_key:
                neuron_def[neuron_types]["morphology_key"] = morphology_key

            if modulation_key:
                neuron_def[neuron_types]["modulation_key"] = modulation_key

        else:
            neuron_def = self.gather_all_neurons(neuron_types=neuron_types, all_combinations=all_combinations)

        if reaction_diffusion_file is not None:
            for neuron_key in neuron_def.keys():
                neuron_def[neuron_key]["reaction_diffusion"] = reaction_diffusion_file

        # Just generate a set of points
        region_def[vol_name]["volume"]["n_putative_points"] = max(len(neuron_def.keys())*5, 10000)

        fake_axon_density = ["r", "1", 10e-6]

        for n in neuron_def.keys():
            neuron_def[n]["num_neurons"] = num_replicas
            neuron_def[n]["volume_id"] = vol_name
            neuron_def[n]["rotation_mode"] = "random"
            neuron_def[n]["hoc"] = None

            if isinstance(neuron_def[n]["neuron_path"], str):
                neuron_def[n]["neuron_path"] = {n: neuron_def[n]["neuron_path"]}

            if not self.has_axon(neuron_def[n]):
                print(f"One or more of morphologies {neuron_def[n]['morphology']} has no axon, faking it.")
                # We will have no connections in this test network, so add empty density
                neuron_def[n]["axon_density"] = fake_axon_density

        region_def[vol_name]["neurons"] = neuron_def

        config_def["regions"] = region_def

        return config_def

    def load_network_info(self, network_file):

        self.neuron_info = SnuddaLoad(network_file).data["neurons"]

    def collect_neurons(self):

        if not self.neuron_info:
            self.load_network_info(network_file=self.network_file)

        neuron_id = np.array([x["neuron_id"] for x in self.neuron_info])
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
                            n_input_min, n_input_max, num_input_steps,
                            input_frequency,
                            input_correlation,
                            synapse_density,
                            synapse_conductance,
                            synapse_parameter_file,
                            cluster_size,
                            cluster_spread,
                            input_duration=10.0):

        # assert n_input_min > 0, "No point using n_input_min=0, please instead use input_frequency 0."

        neuron_sets = self.collect_neurons()
        n_inputs = dict()

        for neuron_type in neuron_sets:

            # For each neuron model we will have num_replicas copies (see other part of code), and this
            # will determine how many steps we have between n_input_min and n_input_max

            neuron_id_list = neuron_sets[neuron_type]
            n_unique = len(neuron_id_list) / num_input_steps

            assert n_unique % 1 == 0, f"Internal error, every model version should exist {num_input_steps} times"

            num_range = np.linspace(n_input_min, n_input_max, num=num_input_steps).astype(int).reshape((1, num_input_steps)).repeat(int(n_unique), axis=0).flatten()

            for neuron_id, num_input in zip(neuron_id_list, num_range):

                if type(synapse_density) == dict:
                    sd = synapse_density[neuron_type]
                else:
                    sd = synapse_density

                if type(synapse_parameter_file) == dict:
                    if neuron_type in synapse_parameter_file:
                        spf = synapse_parameter_file[neuron_type]
                    else:
                        spf = None
                        print(f"No parameter file for {neuron_type} {input_type} input")
                else:
                    spf = synapse_parameter_file

                if type(synapse_conductance) == dict:
                    sc = synapse_conductance[neuron_type]
                else:
                    sc = synapse_conductance

                self.add_input(input_target=neuron_id,
                               input_type=input_type,
                               input_frequency=input_frequency,
                               input_correlation=input_correlation,
                               input_duration=input_duration,
                               input_density=sd,
                               num_input=num_input,
                               input_conductance=sc,
                               synapse_parameter_file=spf,
                               cluster_size=cluster_size,
                               cluster_spread=cluster_spread)

        with open(input_config_file, "w") as f:
            json.dump(self.input_info, f, indent=4, cls=NumpyEncoder)

    def add_input(self, input_target, input_type, input_frequency, input_correlation,
                  input_duration,
                  input_density, num_input, input_conductance,
                  cluster_size, cluster_spread,
                  synapse_parameter_file):

        if type(input_target) != str:
            input_target = str(input_target)

        if not self.input_info:
            self.input_info = collections.OrderedDict()

        assert type(input_frequency) == list, f"add_input: input_freq should be a list: {input_frequency}"
        n_steps = len(input_frequency)
        input_start = input_duration * np.arange(0, n_steps)
        input_end = input_duration * np.arange(1, n_steps + 1)

        if input_target not in self.input_info:
            # By only creating dictionary if it does not exist, we can have multiple inputs
            self.input_info[input_target] = collections.OrderedDict()

        self.input_info[input_target][input_type] = collections.OrderedDict()

        self.input_info[input_target][input_type]["generator"] = "poisson"
        self.input_info[input_target][input_type]["type"] = "AMPA_NMDA"
        self.input_info[input_target][input_type]["synapse_density"] = input_density
        self.input_info[input_target][input_type]["num_inputs"] = num_input
        self.input_info[input_target][input_type]["frequency"] = input_frequency
        if input_correlation is not None:
            self.input_info[input_target][input_type]["correlation"] = input_correlation
        self.input_info[input_target][input_type]["start"] = input_start
        self.input_info[input_target][input_type]["end"] = input_end
        self.input_info[input_target][input_type]["population_unit_correlation"] = 0.0
        self.input_info[input_target][input_type]["jitter"] = 0.0
        self.input_info[input_target][input_type]["conductance"] = input_conductance
        self.input_info[input_target][input_type]["mod_file"] = "tmGlut"
        if synapse_parameter_file is not None:
            self.input_info[input_target][input_type]["parameter_file"] = synapse_parameter_file

        if cluster_size is not None:
            self.input_info[input_target][input_type]["cluster_size"] = cluster_size

        if cluster_spread is not None:
            self.input_info[input_target][input_type]["cluster_spread"] = cluster_spread

    def plot_generated_input(self, num_bins=50):
        # This function just checks that we have reasonable spikes generated

        import matplotlib.pyplot as plt

        input_spike_data = h5py.File(self.input_spikes_file, 'r')
        network_data = h5py.File(self.network_file, 'r')

        neuron_type = np.array([x.decode().split("_")[0].lower() for x in network_data["network/neurons/name"]])
        neuron_id = np.array([x for x in network_data["network/neurons/neuron_id"]])

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

                    distance_to_soma[input_type].append(
                        input_spike_data["input"][str(nid)][input_type]["distance_to_soma"][:])

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
                ax2.hist(np.concatenate(distance_to_soma[input_type]) * 1e6, 50, histtype="step")

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

    def simulate(self, mech_dir=None, sample_dt=0.005):

        if self.input_seed_list is None:
            self.simulate_helper(mech_dir=mech_dir, sample_dt=sample_dt)
        else:
            for ctr, (input_file, output_file) in enumerate(zip(self.input_spikes_file, self.output_file)):
                print(f"Iteration: {ctr+1}/{len(self.input_spikes_file)}")
                print(f"Input file: {input_file}\nOutput file: {output_file}")
                self.simulate_helper(mech_dir=mech_dir, sample_dt=sample_dt,
                                     input_spikes_file=input_file, output_file=output_file)

    def simulate_helper(self, mech_dir=None, sample_dt=0.01, input_spikes_file=None, output_file=None):

        if input_spikes_file is None:
            input_spikes_file = self.input_spikes_file

        if output_file is None:
            output_file = self.output_file

        from snudda.core import Snudda
        Snudda.compile_mechanisms(mech_dir=mech_dir)

        # Get info so we can set max_time correctly
        self.read_tuning_info()

        from neuron import h  # , gui
        start = timeit.default_timer()

        pc = h.ParallelContext()

        sim = SnuddaSimulate(network_file=self.network_file,
                             input_file=input_spikes_file,
                             output_file=output_file,
                             sample_dt=sample_dt)
        sim.setup()
        sim.add_external_input()
        sim.check_memory_status()

        # sim.add_volt_recording()
        sim.add_volt_recording_soma()

        t_sim = self.max_time * 1000  # Convert from s to ms for Neuron simulator

        sim.check_memory_status()
        print(f"Running simulation for {t_sim} ms.")
        sim.run(t_sim)  # In milliseconds

        print("Simulation done, saving output")
        sim.write_output()

        stop = timeit.default_timer()
        if sim.pc.id() == 0:
            print(f"Program run time: {stop - start:.1f}s")

        print("About to clear memory")
        sim.check_memory_status()
        sim.clear_neuron()
        print("Memory cleared")
        sim.check_memory_status()
        sim = None

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
            self.max_time = tuning_meta_data["max_time"]
            self.input_duration = tuning_meta_data["input_duration"]
            self.frequency_range = tuning_meta_data["frequency_range"]
            self.neurons_path = tuning_meta_data["neurons_directory"]
        except:
            print(f"Failed to read {tuning_info_file}")

    def write_tuning_info(self):
        tuning_meta_data = collections.OrderedDict()
        tuning_meta_data["input_duration"] = self.input_duration
        tuning_meta_data["max_time"] = self.max_time
        tuning_meta_data["frequency_range"] = self.frequency_range
        tuning_meta_data["neurons_directory"] = self.neurons_path

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
    parser.add_argument("--mechDir", type=str, help="Path to mechanisms", default=None)
    parser.add_argument("--input_type", help="Type of external input",
                        choices=["thalamic", "cortical", "corticothalamic"],
                        default="thalamic")  # only use corticalthalamic in analyse
    parser.add_argument("--numInputSteps", type=int, help="Number of steps for number of inputs to neurons",
                        default=10)
    parser.add_argument("--numInputMin", type=int, help="Minimum number of synaptic inputs of input_type", default=100)
    parser.add_argument("--numInputMax", type=int, help="Maximum number of synaptic inputs of input_type", default=1000)
    parser.add_argument("--inputDuration", type=float, default=10.0,
                        help="Duration of each frequency test, longer need for irregularly firing neurons")
    parser.add_argument("--inputFrequency", type=str, default="[0,1,2,5]",
                        help="Input frequency, float or list of floats")
    parser.add_argument("--neuronType", default=None, type=str,
                        help="Optional, if only we want to simulate one neuron type, eg. FS")
    parser.add_argument("--singleNeuronType", default=None, type=str,
                        help="Optional, if only we want to simulate one neuron subtype, eg. FS_1")
    parser.add_argument("--meta_input", action="store_true", default=False)
    parser.add_argument("--seed_list", type=str, default=None)
    parser.add_argument("--no_downsampling", action="store_true")

    args = parser.parse_args()

    # TODO: Let the user choose input type, duration for each "run", frequency range, number of input range

    if args.seed_list is not None:
        seed_list = ast.literal_eval(args.seed_list)
    else:
        seed_list = None

    input_scaling = InputTuning(args.networkPath, input_seed_list=seed_list)

    if args.action == "setup":
        input_frequency = ast.literal_eval(args.inputFrequency)
        if type(input_frequency) != list:
            input_frequency = np.array(list(input_frequency))

        input_scaling.setup_network(neurons_path=args.neurons,
                                    num_replicas=args.numInputSteps,
                                    neuron_types=args.neuronType,
                                    single_neuron_path=args.singleNeuronType)
        input_scaling.setup_input(input_type=args.input_type,
                                  num_input_min=args.numInputMin,
                                  num_input_max=args.numInputMax,
                                  num_replicas=args.numInputSteps,
                                  input_duration=args.inputDuration,
                                  input_frequency_range=input_frequency,
                                  use_meta_input=args.no_meta_input)

        print("Tip, to run in parallel on your local machine use: "
              "mpiexec -n 4 python3 tuning/input_tuning.py simulate <yournetworkhere>")

    if args.action == "setup_background":
        input_frequency = ast.literal_eval(args.inputFrequency)

        if type(input_frequency) != list:
            input_frequency = np.array(list(input_frequency))

        if len(input_frequency) != 1:
            raise ValueError("input_frequency must only be one value when doing setup_background")

        input_scaling.setup_network(neurons_path=args.neurons,
                                    num_replicas=args.numInputSteps,
                                    neuron_types=args.neuronType,
                                    single_neuron_path=args.singleNeuronType)

        print(f"Setting up background input, will do cortical and thalamic background 50-50 at {input_frequency}")

        input_scaling.setup_background_input(input_types=["cortical_background", "thalamic_background"],
                                             input_density=["1.15*0.05/(1+exp(-(d-30e-6)/5e-6))", "0.05*exp(-d/200e-6)"],
                                             input_fraction=[0.5, 0.5],
                                             num_input_min=args.numInputMin,
                                             num_input_max=args.numInputMax,
                                             input_duration=args.inputDuration,
                                             input_frequency=[input_frequency, input_frequency])

    elif args.action == "simulate":
        print("Run simulation...")
        print("Tip, to run in parallel on your local machine use: "
              "mpiexec -n 4 python3 tuning/input_tuning.py simulate <yournetworkhere>")

        if args.no_downsampling:
            sample_dt = None
        else:
            sample_dt = 0.01

        input_scaling.simulate(mech_dir=args.mechDir, sample_dt=sample_dt)

    elif args.action == "analyse":
        # input_scaling.plot_generated_input()
        input_scaling.analyse_results(input_type=args.input_type)

    elif args.action == "analyse_background":
        input_scaling.analyse_background()

    else:
        print(f"Unknown action {args.action}")

    # python3 input_tuning/input_tuning.py setup networks/input_scaling_v1/ data/neurons/striatum/
    # mpiexec -n 4 python3 input_tuning/input_tuning.py simulate networks/input_scaling_v1/ < input.txt &> output-tuning.txt &

    #
