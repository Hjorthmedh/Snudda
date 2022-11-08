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
        print("Entered PlotInputVsOutput")
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
            network_file = os.path.join(self.network_path, "network-synapses.hdf5")
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
        print("End of PlotInputVsOutput")
    def make_figures_directory(self):

        fig_dir = os.path.join(self.network_path, "figures")

        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
    def plot_cur_input_vs_output(self,experiment_config_file=None):
        print("Entered plot_cur_input_vs_output")
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
        print("Gathering lists of spike rates for each neuron subtype.")
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
        print("Gathering completed.", flush=True)

        for unique_neuron_subtype in unique_neuron_subtypes:
            print("Plotting for: " + unique_neuron_subtype, flush=True)
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
    def plot_input_output_corr(self,experiment_config_file=None):
        print("Entered plot_cur_input_vs_output2")
        self.make_figures_directory()
        #experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "experiment_config","experiment-config-6-bobek.json")
        with open(experiment_config_file, "r") as f:
            exp_proto = json.load(f, object_pairs_hook=OrderedDict)
        start_times = exp_proto["iSPN"]["CorticalSignal"]["start"]
        end_times = exp_proto["iSPN"]["CorticalSignal"]["end"]
        correlation = exp_proto["iSPN"]["CorticalSignal"]["populationUnitCorrelation"]
        spike_rate_lists = []
        spike_rate_list = []
        num_neurons = len(self.spike_time.keys())
        binSize = 5
        dt = 1#This is sample size actually?
        T = (np.array(end_times) - np.array(start_times))[0]*1000
        binKernel = np.ones(int(binSize / dt))

        # spikes1 = self.spike_time[0][0]
        # spikes1=spikes1*1000
        # binspk1 = self.spk2bin(spikes1, np.zeros(int(T / dt) + 1), binKernel, dt)
        # spikes2 = self.spike_time[1][0]
        # spikes2 = spikes2 * 1000
        # binspk2 = self.spk2bin(spikes2, np.zeros(int(T / dt) + 1), binKernel, dt)
        # spike_mat = np.array([binspk1,binspk2])
        # cor = self.spikeCoin(spike_mat)
        # np.corrcoef(binspk1, binspk2)
        cors = []
        reset_time = T#First stimulus starts after initial pause
        for start_time, end_time in zip(start_times, end_times):
            reset_time = start_time*1000#+= (start_time - end_time_prev + end_time - start_time) * 1000  # Next correlation stim starts after stimulus and pause period
            #binspks=np.array([])
            binspks=[]
            #spike_rate_list.append(len(spikes1[(spikes1 >= start_time) & (spikes1 < end_time)]))
            fig2, ax2 = plt.subplots()
            binspk_sum=np.array([0])
            for nrni in range(0, num_neurons):
                spikes = self.spike_time[nrni][0]
                spikes_within = spikes[(spikes >= start_time) & (spikes < end_time)]
                spikes_within = spikes_within * 1000 #Seconds to millisecons
                spikes_within = spikes_within - reset_time
                binspk = self.spk2bin(spikes_within, np.zeros(int(T / dt) + 1), binKernel, dt)
                binspks.append(binspk)
                #binspk_sum=binspk_sum+binspk
                #ax2.plot(binspk)# color='black'
            binspks = np.array(binspks)
            ax2.plot(sum(binspks))  # color='black'
            #cor = self.spikeCoin(binspks)
            #cors.append(cor[1])
            cor = np.corrcoef(binspks)
            cor = np.array(cor)
            cormean = cor[np.triu_indices_from(cor, 1)].mean()
            cors.append(cormean)
            print("midloop")

        fig, ax = plt.subplots()

        ax.plot(correlation,cors, marker='.',color='black')
        #ax.set_ylabel('Spikes/s')
        #ax.set_xlabel('Current injection (A)')
        figure_path = os.path.join(self.network_path, "figures", "input_output_correlation.png")
        fig.savefig(figure_path, dpi=300)
        plt.close(fig)
        print("end")
        # neuron_subtypes=self.output_load.get_neuron_name()#['dSPN_1', 'dSPN_0', 'dSPN_3', 'dSPN_2']
        # para_keys=self.output_load.get_parameter_keys()
        # morph_keys=self.output_load.get_morphology_keys()
        # #self.network_info.get_neuron_id_of_type('dSPN')#array([0, 1, 2, 3])
        # num_subtypes=len(exp_proto["currentInjection"])
        # spike_rate_lists=[]
        # print("Gathering lists of spike rates for each neuron subtype.")
        # for i in np.arange(0,num_subtypes):
        #     start_times=exp_proto["currentInjection"][i]["start"]
        #     end_times=exp_proto["currentInjection"][i]["end"]
        #     cur_steps=exp_proto["currentInjection"][i]["amplitude"]
        #     n_steps=len(cur_steps)
        #     nrn_id=exp_proto["currentInjection"][i]["neuronID"]
        #     spikes=self.spike_time[nrn_id][0]
        #     spike_rate_list=[]
        #     for start_time,end_time in zip(start_times, end_times):
        #         #step_period=(self.time>start_time) and (self.time<end_time)
        #         #membr_pot=self.voltage[step_period]
        #         spike_rate_list.append(len(spikes[(spikes >= start_time) & (spikes < end_time)]))
        #     spike_rate_lists.append(spike_rate_list)
        # neuron_subtypes = np.array(neuron_subtypes)
        # unique_neuron_subtypes=np.unique(neuron_subtypes)
        # spike_rate_lists=np.array(spike_rate_lists)
        # print("Gathering completed.", flush=True)
        #
        # for unique_neuron_subtype in unique_neuron_subtypes:
        #     print("Plotting for: " + unique_neuron_subtype, flush=True)
        #     fig, ax = plt.subplots()
        #     bool_vec=unique_neuron_subtype == neuron_subtypes
        #     sub_spike_rate_lists=spike_rate_lists[bool_vec]
        #     for sub_spike_rate_list in sub_spike_rate_lists:
        #         ax.plot(cur_steps,sub_spike_rate_list,marker='.',color='gray')
        #     avg_spike_rate_list=np.mean(sub_spike_rate_lists,0)
        #     ax.plot(cur_steps, avg_spike_rate_list, marker='.',color='black')
        #     ax.set_ylabel('Spikes/s')
        #     ax.set_xlabel('Current injection (A)')
        #     figure_path = os.path.join(self.network_path, "figures", unique_neuron_subtype + '_' +
        #                                "spike_rate_vs_current_input.png")
        #     fig.savefig(figure_path, dpi=300)
        #     plt.close(fig)
        return cors, correlation

    def spk2bin(self, signal, template, kernel, dt):
        """Binning spikes from spiking times.

        Parameters
        ----------
        signal : ndarray
            Spiking times.
        template : ndarray
            Binning template.
        kernel : ndarray
            Spike convolution kernel.

        Returns
        -------
        ndarray
            Binned signal from spiking times.

        """
        if signal is []:
            return template
        else:
            for s in signal:
                template[int(s/dt)-1] = 1
            return np.convolve(template, kernel, 'same')
    def spikeCoin(self, spikeMat):
        """Cout coincidence of several spike trains.

        Parameters
        ----------
        spikeMat : numpy matrix
            A matrix of binned spiking actitivities.

        Returns
        -------
        type
            Description of returned object.

        """
        # Get number of trials
        N, _ = spikeMat.shape

        # Calculate coincidence
        coV_ori = np.dot(spikeMat, np.transpose(spikeMat))
        coV = coV_ori.copy()

        # Normalize to geometric mean energy
        for i in range(N):
            for j in range(N):
                energy = coV_ori[i,i]*coV_ori[j,j]
                if energy:
                    coV[i,j] /= np.sqrt(energy)
                else:
                    coV[i,j] = 0

        # Correlation
        cor = (np.sum(coV) - np.trace(coV))/N/(N-1)

        return coV, cor

if __name__ == "__main__":
    print("entered mainscript of plot_cur_input_vs_output.py")
    from argparse import ArgumentParser, RawTextHelpFormatter
    parser = ArgumentParser("Plot network", formatter_class=RawTextHelpFormatter)
    parser.add_argument("network_path",type=str, help="Network path")
    parser.add_argument("--network_file",type=str, default="network-synapses.hdf5")
    parser.add_argument("--network_output",type=str, default="output.hdf5")
    parser.add_argument("--experiment_config_file", type=str , help="")
    args = parser.parse_args()
    print("making object of PlotInputVsOutput")
    pio = PlotInputVsOutput(args.network_path, args.network_output, network_file=args.network_file)
    print("entering plot_cur_input_vs_output()")
    pio.plot_cur_input_vs_output(experiment_config_file=args.experiment_config_file)
    
    
