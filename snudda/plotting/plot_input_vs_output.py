# python3 plot_traces.py save/simulation/network-output.hdf5 save/network-synapses.hdf5


import sys
import os,json
from collections import OrderedDict
import h5py
import numpy as np
from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation
import re
import ntpath
import time
import matplotlib.pyplot as plt

class PlotInputVsOutput:
    def __init__(self, network_path, output_file, network_file=None, input_file=None):
        self.network_path = network_path
        self.output_file = output_file
        self.network_file = network_file
        self.time = []
        self.voltage = dict([])
        self.neuron_name_remap = {"FSN": "FS"}
        try:
            self.ID = int(re.findall('\d+', ntpath.basename(output_file))[0])
        except:
            print("Unable to guess ID, using 666.")
            self.ID = 666

        if network_file is None and "simulation" in output_file:
            network_path = os.path.dirname(os.path.dirname(output_file))
            network_file = os.path.join(network_path, "network-synapses.hdf5")
            if os.path.exists(network_file):
                self.network_file = network_file

        if network_file is not None:
            network_path = os.path.dirname(os.path.dirname(network_file))
        else:
            network_path = None

        if self.network_file is not None:
            print(f"Loading network info from {self.network_file}")
            self.network_info = SnuddaLoad(self.network_file)
        else:
            self.network_info = None

        self.output_load = SnuddaLoadNetworkSimulation(network_simulation_output_file=output_file,
                                                       network_path=network_path)
        self.spike_time = self.output_load.get_spikes()

        self.voltage = self.output_load.get_voltage()
        self.time = self.output_load.get_time()
    def make_figures_directory(self):

        fig_dir = os.path.join(self.network_path, "figures")

        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
    def get_input_vs_output(self, experiment_config_file=None):
        self.make_figures_directory()
        with open(experiment_config_file, "r") as f:
            exp_proto = json.load(f, object_pairs_hook=OrderedDict)
        start_times = exp_proto["iSPN"]["CorticalSignal"]["start"]
        end_times = exp_proto["iSPN"]["CorticalSignal"]["end"]
        in_freqs = exp_proto["iSPN"]["CorticalSignal"]["frequency"]
        num_neurons = len(self.spike_time.keys())


        out_freqs_all = []
        for start_time, end_time in zip(start_times, end_times):
            out_freqs = []
            for nrni in range(0, num_neurons):
                spikes = self.spike_time[nrni][0]
                spikes_within = spikes[(spikes >= start_time) & (spikes < end_time)]
                spike_rate = len(spikes_within)/(end_time-start_time)
                out_freqs.append(spike_rate)
            out_freqs_all.append(out_freqs)
        return in_freqs, out_freqs_all

    def plot_cur_input_vs_output(self,experiment_config_file=None):
        self.make_figures_directory()
        #experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "experiment_config","experiment-config-6-bobek.json")
        with open(experiment_config_file, "r") as f:
            exp_proto = json.load(f, object_pairs_hook=OrderedDict)
        end_time=exp_proto["meta"]["simulationDuration"]
        neuron_subtypes=self.output_load.get_neuron_name()#['dSPN_1', 'dSPN_0', 'dSPN_3', 'dSPN_2']
        para_keys=self.output_load.get_parameter_keys()
        morph_keys=self.output_load.get_morphology_keys()
        #self.network_info.get_neuron_id_of_type('dSPN')#array([0, 1, 2, 3])
        num_subtypes=len(exp_proto["currentInjection"])
        spike_rate_lists=[]
        for i in np.arange(0,num_subtypes):
            start_times=exp_proto["currentInjection"][i]["start"]
            end_times=exp_proto["currentInjection"][i]["end"]
            cur_steps=exp_proto["currentInjection"][i]["amplitude"]
            n_steps=len(cur_steps)
            nrn_id=exp_proto["currentInjection"][i]["neuronID"]
            spikes=self.spike_time[nrn_id][0]
            spike_rate_list=[]
            for start_time,end_time in zip(start_times, end_times):
                #step_period=(self.time>start_time) and (self.time<end_time)
                #membr_pot=self.voltage[step_period]
                spike_rate_list.append(len(spikes[(spikes >= start_time) & (spikes < end_time)]))
            spike_rate_lists.append(spike_rate_list)
        neuron_subtypes = np.array(neuron_subtypes)
        unique_neuron_subtypes=np.unique(neuron_subtypes)
        spike_rate_lists=np.array(spike_rate_lists)

        for unique_neuron_subtype in unique_neuron_subtypes:
            fig, ax = plt.subplots()
            bool_vec=unique_neuron_subtype == neuron_subtypes
            sub_spike_rate_lists=spike_rate_lists[bool_vec]
            for sub_spike_rate_list in sub_spike_rate_lists:
                ax.plot(cur_steps,sub_spike_rate_list,marker='.',color='gray')
            avg_spike_rate_list=np.mean(sub_spike_rate_lists,0)
            ax.plot(cur_steps, avg_spike_rate_list, marker='.',color='black')
            ax.set_ylabel('Spikes/s')
            ax.set_xlabel('Current injection (A)')
            figure_path = os.path.join(self.network_path, "figures", unique_neuron_subtype + '_' +
                                       "spike_rate_vs_current_input.png")
            fig.savefig(figure_path, dpi=300)
            plt.close(fig)

if __name__ == "__main__":
    print("Calling this way is not ready")
