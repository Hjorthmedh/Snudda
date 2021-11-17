# In Chuhma et al 2011 they stimulate 10% of MS population and record the
# response in MS, ChIN and FS.
#
# Let us do the same thing but in our model network.
#
# Before you run this network you need to create a network with reduced MS
# density, no point in simulating the silent 90%.
#
# Suggestion: create a slice that is 0.6 x 0.6 mm in length, and 0.15 mm deep.
#
#
# Example usage:
#
# python3 snudda_model_current_injections.py setup Chuhma2011 networks/Chuhma2011-v15 
#
# mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_model_current_injections.py run Chuhma2011 networks/Chuhma2011-v15/ 
#
# python3 snudda_model_current_injections.py analyse Chuhma2011 networks/Chuhma2011-v15/ 
#
# OR
#
#
# python3 snudda_model_current_injections.py setup Straub2016FS networks/Straub2016FS-v9
# mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_model_current_injections.py run Straub2016FS networks/Straub2016FS-v9
# python3 snudda_model_current_injections.py analyse Straub2016FS networks/Straub2016FS-v9
#
#

import os
import sys

from snudda.simulate.simulate import SnuddaSimulate
from snudda.utils.load import SnuddaLoad
from snudda.init.init import SnuddaInit
from snudda import Snudda

import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import neuron


class SnuddaModelCurrentInjections(object):

    def __init__(self,
                 sim_name,
                 sim_type,
                 cur_inj=10e-9):

        self.sim_name = sim_name
        self.sim_type = sim_type
        self.snudda_sim = None
        self.exp_data_dict = dict()
        self.exp_trace_dict = dict()

        if sim_type == "Chuhma2011":
            self.t_inj = 0.3
            self.inj_duration = 1e-3
            self.cur_inj = 10e-9
            self.t_window = 0.03
            self.sim_end = self.t_inj + self.t_window * 2
            self.hold_v = -60e-3
            self.gaba_rev = 2e-3  # 144mM inside, 133.5mM outside, 32-33C --NERNST--> 2mV
            self.n_nrns = 100  # How many we measure from of each type

        elif sim_type == "Straub2016FS" or sim_type == "Straub2016LTS":
            self.t_inj = 0.5
            self.inj_duration = 1e-3
            self.cur_inj = 10e-9
            self.t_window = 0.03
            self.sim_end = self.t_inj + self.t_window * 2
            self.hold_v = -70e-3
            self.gaba_rev = -0.3e-3  # Out: 133.5 mM chloride, In 131.5 mM, Temperature 33-34 C
            self.n_nrns = 30
        elif sim_type == "Szydlowski2013":

            assert False, "Szydlowski did not stimulate all FS or did they when recording LTS?"
            self.t_inj = 0.5
            self.inj_duration = 1e-3
            self.cur_inj = 10e-9
            self.t_window = 0.03
            self.sim_end = self.t_inj + self.t_window * 2
            self.hold_v = -70e-3
            self.GABA_rev = -30e-3
            self.n_nrns = 30

        else:
            print(f"Unknown simType: {sim_type}")

        self.plot_exp_trace = False

    ############################################################################

    def define_network(self, sim_name, sim_type=None):

        if sim_type is None:
            sim_type = self.sim_type

        config_name = os.path.join(sim_name, "network-config.json")
        cnc = SnuddaInit(struct_def={}, config_file=config_name)

        # In a 1x1x0.15 mm slice there are 12000 neurons normally
        # We want 10% of MS population only, since those are the ones being
        # stimulated (Chuhma et al 2011)
        #
        # 47.5% dSPN normally, now 4.75% of normal density = 570 dSPN, 570 iSPN
        # We assume we measure only monosynaptic connections.
        #
        # Adding 10 FS, 10 ChIN, 10 dSPN, 10 iSPN to measure from
        #

        if False:
            # Small debug version
            # cnc.defineStriatum(nMSD1=20,nMSD2=20,nFS=0,nLTS=0,nChIN=0,
            #                   volumeType="slice",sideLen=200e-6)
            # cnc.defineStriatum(nMSD1=20,nMSD2=20,nFS=10,nLTS=0,nChIN=10,
            #                   volumeType="slice",sideLen=200e-6)
            cnc.define_striatum(num_dSPN=153, num_iSPN=153, num_FS=10, num_LTS=0, num_ChIN=10,
                                volume_type="slice", side_len=500e-6)

        if sim_type == "Chuhma2011":
            cnc.define_striatum(num_dSPN=1140 + self.n_nrns,
                                num_iSPN=1140 + self.n_nrns,
                                num_FS=5,
                                num_LTS=0,
                                num_ChIN=self.n_nrns,
                                volume_type="slice",
                                side_len=1000e-6,
                                slice_depth=300e-6)  # 400mum, assume 100 mum dead

        elif sim_type == "Straub2016FS":
            # nFS must be correct density, but readout neurons can be any density
            cnc.define_striatum(num_dSPN=self.n_nrns,
                                num_iSPN=self.n_nrns,
                                num_FS=182, num_LTS=0,
                                num_ChIN=self.n_nrns,
                                volume_type="slice",
                                side_len=1000e-6,
                                slice_depth=175e-6)  # 275e-6 m slice, assume 100e-6 dead

        elif sim_type == "Straub2016LTS":
            cnc.define_striatum(num_dSPN=self.n_nrns,
                                num_iSPN=self.n_nrns,
                                num_FS=0, num_LTS=98,
                                num_ChIN=self.n_nrns,
                                volume_type="slice",
                                side_len=1000e-6,
                                slice_depth=175e-6)
        elif sim_type == "Szydlowski2013":
            cnc.define_striatum(num_dSPN=0,
                                num_iSPN=0,
                                num_FS=156,
                                num_LTS=self.n_nrns,
                                num_ChIN=0,
                                volume_type="slice",
                                side_len=1000e-6,
                                slice_depth=150e-6)
        else:
            print("setup : Unkown simType: " + str(sim_type))
            sys.exit(-1)

        dir_name = os.path.dirname(config_name)

        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        cnc.write_json(config_name)

    ############################################################################

    def simulate_network(self, sim_name, sim_type=None):

        if sim_type is None:
            sim_type = self.sim_type

        if sim_type == "Chuhma2011":
            self.simulate_network_chuhma2011(sim_name)

        elif sim_type == "Straub2016FS" or sim_type == "Straub2016LTS":
            self.simulateNetworkStraub2016(sim_name, sim_type)

        elif sim_type == "Szydlowski2013":
            self.simulate_network_szydlowski_2013(sim_name)

        else:
            print("simulateNetwork: unknown simType = " + str(sim_type))
            sys.exit(-1)

    ############################################################################

    def simulate_network_chuhma2011(self, sim_name):

        sim_type = "Chuhma2011"

        if self.snudda_sim is None:
            log_file = sim_name + "/log/simlog.txt"

            cut_file = sim_name + "/network-cut-slice.hdf5"
            if os.path.exists(cut_file):
                self.network_file = cut_file
            else:
                self.network_file = sim_name + "/network-synapses.hdf5"

            print("Using network file: " + str(self.network_file))

            self.snudda_sim = SnuddaSimulate(network_file=self.network_file,
                                             input_file=None,
                                             log_file=log_file,
                                             disable_gap_junctions=True)

        # Get neuronID of neurons that will get artificial stimulation
        stim_id = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"] if "SPN" in x["type"]]

        # Pick neurons that we will measure from
        measure_fs = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"] if x["type"] == "FS"]

        measure_chin = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"] if x["type"] == "ChIN"]

        # For future: maybe pick the ones more centrally
        measure_dspn = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"]
                        if x["type"] == "dSPN"][:self.n_nrns]
        measure_ispn = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"]
                        if x["type"] == "iSPN"][:self.n_nrns]

        # Remove the overlap, ie dont stimulate the neurons we measure from
        stim_id = np.setdiff1d(stim_id, np.union1d(measure_dspn, measure_ispn))

        measure_id = np.union1d(np.union1d(measure_dspn, measure_ispn), np.union1d(measure_fs, measure_chin))

        self._simulate_network_helper(sim_name, sim_type, stim_id, measure_id)

    ############################################################################

    def simulateNetworkStraub2016(self, sim_name, sim_type):

        if self.snudda_sim is None:
            log_file = os.path.join(sim_name, "log", "simlog.txt")
            self.network_file = sim_name + "/network-synapses.hdf5"

            self.snudda_sim = SnuddaSimulate(network_file=self.network_file,
                                             input_file=None,
                                             log_file=log_file,
                                             disable_gap_junctions=True)

        # Get neuronID of neurons that will get artificial stimulation
        if sim_type == "Straub2016FS":
            stim_id = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"] if "FS" in x["type"]]
        elif sim_type == "Straub2016LTS":
            stim_id = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"] if "LTS" in x["type"]]
        else:
            print(f"simulateNetworkStraub2016: Unknown simType : {sim_type}")
            sys.exit(-1)

        measure_id = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"]
                      if x["type"] == "ChIN" or x["type"] == "dSPN" or x["type"] == "iSPN"]

        self._simulate_network_helper(sim_name, sim_type, stim_id, measure_id)

    ############################################################################

    def simulate_network_szydlowski_2013(self, sim_name):

        if self.snudda_sim is None:
            log_file = os.path.join(sim_name, "log", "simlog.txt")
            self.network_file = os.path.join(sim_name, "network-synapses.hdf5")

            self.snudda_sim = SnuddaSimulate(network_file=self.network_file,
                                             input_file=None,
                                             log_file=log_file,
                                             disable_gap_junctions=True)

        stim_id = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"] if "FS" in x["type"]]
        measure_id = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"] if x["type"] == "LTS"]

        self._simulate_network_helper(sim_name, sim_type, stim_id, measure_id)

    ############################################################################

    def _simulate_network_helper(self, sim_name, sim_type, stim_id, measure_id):

        if self.snudda_sim is None:
            log_file = os.path.join(sim_name, "log", "simlog.txt")
            self.network_file = os.path.join(sim_name, "network-synapses.hdf5")

            self.snudda_sim = SnuddaSimulate(network_file=self.network_file,
                                             input_file=None,
                                             log_file=log_file,
                                             disable_gap_junctions=True)

        # Set up stimulation protocol
        for n_id in stim_id:
            self.snudda_sim.add_current_injection(neuron_id=n_id,
                                                  start_time=self.t_inj,
                                                  end_time=self.t_inj + self.inj_duration,
                                                  amplitude=self.cur_inj)

        # Add recordings
        self.snudda_sim.add_voltage_clamp(cell_id=measure_id,
                                          voltage=self.hold_v,
                                          duration=self.sim_end,
                                          save_i_flag=True)

        # Also add voltage recording for debugging reasons
        save_voltage = True  # False #True
        if save_voltage:
            self.snudda_sim.add_recording(cell_id=stim_id)

        self.set_gaba_rev(self.gaba_rev)

        self.snudda_sim.run(self.sim_end * 1e3)

        self.current_file = os.path.join(sim_name, f"{sim_type}-network-stimulation-current.txt")

        self.snudda_sim.write_current(self.current_file)

        if save_voltage:
            voltage_file = os.path.join(sim_name,  f"{sim_type}-network-stimulation-voltage.txt")

            self.snudda_sim.write_voltage_OLD(voltage_file)

    ############################################################################

    def create_network(self, sim_name):

        sn = Snudda(sim_name)

        class FakeArgs(object):
            def __init__(self):
                setattr(self, "h5legacy", "latest")
                setattr(self, "volumeID", None)
                setattr(self, "hvsize", None)  # Use default value
                setattr(self, "cont", None)
                setattr(self, "mergeonly", False)

        args = FakeArgs()

        sn.place_neurons(args)
        sn.touch_detection(args)
        sn.prune_synapses(args)

    ############################################################################

    def setup_exp_data_dict(self):

        self.exp_data_dict = dict()

        # Straub et al 2016
        lts_2_spn = np.array([0.0316, 0.0433, 0.0474, 0.1253, 0.1839,
                              0.1860, 0.1946, 0.1968, 0.2082, 0.2203,
                              0.2384, 0.2439, 0.2793, 0.3091, 0.3234,
                              0.3271, 0.3383, 0.3500, 0.3540, 0.3662,
                              0.3662, 0.3831, 0.4053, 0.4099, 0.4288,
                              0.4337, 0.4966, 0.5023, 0.5196, 0.5314,
                              0.5436, 0.5560, 0.5817, 0.6017, 0.6736,
                              0.6968, 0.7047, 0.7047, 0.7127, 0.7979,
                              0.9034, 1.0461])

        lts_2_chin = np.array([0.2466, 0.5080, 0.5196, 0.6017, 0.6660,
                               0.7541, 0.7713, 0.8442, 1.1069, 1.2391,
                               1.2818, 1.4030, 2.3315])

        fs_2_spn = np.array([0.3091, 0.5137, 0.5255, 0.5687, 0.6890,
                             0.8161, 0.8832, 0.8932, 0.9667, 1.0228,
                             1.0228, 1.0822, 1.1844, 1.2391, 1.2964,
                             1.3111, 1.4189, 1.4350, 1.5530, 1.6247,
                             1.7385, 1.7984, 1.9028, 2.1063, 2.2539,
                             2.3580, 2.4669, 2.4949, 2.5232, 2.7307,
                             2.7930, 2.8247, 3.2711, 3.3458, 3.4222,
                             4.2648, 4.4617, 4.9668, 5.3148])

        fs_2_chin = np.array([0.0233, 0.0378, 0.0419, 0.0428, 0.0666,
                              0.0762])

        self.exp_data_dict[("Straub2016LTS", "dSPN")] = lts_2_spn
        self.exp_data_dict[("Straub2016LTS", "iSPN")] = lts_2_spn
        self.exp_data_dict[("Straub2016LTS", "ChIN")] = lts_2_chin
        self.exp_data_dict[("Straub2016FS", "dSPN")] = fs_2_spn
        self.exp_data_dict[("Straub2016FS", "iSPN")] = fs_2_spn
        self.exp_data_dict[("Straub2016FS", "ChIN")] = fs_2_chin

        if self.plot_exp_trace:
            self.exp_trace_dict = dict()

            lts_2_spn = np.genfromtxt("DATA/Straub2016/LTSItoSPN_Straub2.txt")
            lts_2_chin = np.genfromtxt("DATA/Straub2016/LTSItoChIN_Straub2.txt")
            fs_2_spn = np.genfromtxt("DATA/Straub2016/FSItoSPN_Straub2_shorter.txt")

            # Convert current from pA to nA
            lts_2_spn[:, 1:] = 1e-3 * lts_2_spn[:, 1:]
            lts_2_chin[:, 1:] = 1e-3 * lts_2_chin[:, 1:]
            fs_2_spn[:, 1:] = 1e-3 * fs_2_spn[:, 1:]

            self.exp_trace_dict[("Straub2016LTS", "dSPN")] = lts_2_spn
            self.exp_trace_dict[("Straub2016LTS", "iSPN")] = lts_2_spn
            self.exp_trace_dict[("Straub2016LTS", "ChIN")] = lts_2_chin
            self.exp_trace_dict[("Straub2016FS", "dSPN")] = fs_2_spn
            self.exp_trace_dict[("Straub2016FS", "iSPN")] = fs_2_spn

            spn_2_spn = np.genfromtxt("DATA/Chuhma2011/SPNtoSPN_Chuhma.txt")
            spn_2_chin = np.genfromtxt("DATA/Chuhma2011/SPNtoChIN_Chuhma.txt")

            # Convert current from pA to nA
            spn_2_spn[:, 1:] = 1e-3 * spn_2_spn[:, 1:]
            spn_2_chin[:, 1:] = 1e-3 * spn_2_chin[:, 1:]

            self.exp_trace_dict[("Chuhma2011", "dSPN")] = spn_2_spn
            self.exp_trace_dict[("Chuhma2011", "iSPN")] = spn_2_spn
            self.exp_trace_dict[("Chuhma2011", "ChIN")] = spn_2_chin

    ############################################################################

    def analyse_network(self, sim_name, sim_type=None, n_plot_max=10):

        fig_dir = os.path.join(sim_name, "figures")
        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)

        if sim_type is None:
            sim_type = self.sim_type

        if sim_type == "Straub2016LTS":
            pre_type = "LTS"
            self.setup_exp_data_dict()
        elif sim_type == "Straub2016FS":
            pre_type = "FS"
            self.setup_exp_data_dict()
        elif sim_type == "Chuhma2011":
            pre_type = "SPN"
            self.setup_exp_data_dict()
        elif sim_type == "Szydlowski2013":
            pre_type = "FS"
        else:
            print(f"Unknown simType : {sim_type}")
            sys.exit(-1)

        print(f"Analysing data in {sim_name}")
        voltFile = os.path.join(sim_name, f"{sim_type}-network-stimulation-current.txt")

        # Read data from file
        data = np.genfromtxt(voltFile, delimiter=",")

        assert (data[0, 0] == -1)  # First column should be time
        time = data[0, 1:] * 1e-3

        current = dict()

        for rows in data[1:, :]:
            c_id = int(rows[0])
            current[c_id] = rows[1:] * 1e-9

        # Data in time, current now

        # Read the network info
        network_file = os.path.join(sim_name, "network-synapses.hdf5")
        self.snudda_load = SnuddaLoad(network_file)
        self.data = self.snudda_load.data

        recorded_neurons = [x for x in current]

        # Group the neurons by type

        if sim_type == "Chuhma2011":
            neuron_type_list = ["dSPN", "iSPN", "FS", "ChIN"]
        elif sim_type == "Straub2016FS" or sim_type == "Straub2016LTS":
            neuron_type_list = ["dSPN", "iSPN", "ChIN"]
        elif sim_type == "Szydlowski2013":
            neuron_type_list = ["LTS"]
        else:
            print(f"simulate: Unknown simType: {sim_type}")
            sys.exit(-1)

        neuron_plot_list = []

        min_time_idx = np.where(time > self.t_inj)[0][0]
        max_time_idx = np.where(time > self.t_inj + self.t_window)[0][0]

        for nt in neuron_type_list:
            id_list = [x for x in current if self.data["neurons"][x]["type"] == nt]
            max_idx = [np.argmax(np.abs(current[x][min_time_idx:max_time_idx]
                                        - current[x][min_time_idx])) + min_time_idx
                       for x in id_list]

            neuron_plot_list.append((id_list, max_idx))

        matplotlib.rcParams.update({'font.size': 22})

        for plot_id, max_idx in neuron_plot_list:

            if len(plot_id) == 0:
                continue

            plot_type = self.data["neurons"][plot_id[0]]["type"]
            fig_name = os.path.join(fig_dir, f"{sim_type}-{plot_type}-current-traces.pdf")
            fig_name_hist = os.path.join(fig_dir, f"{sim_type}-{plot_type}-current-histogram.pdf")

            good_max = []

            plt.figure()

            peak_amp = []
            peak_time = []
            volt_curve = []

            for p_id, m_idx in zip(plot_id, max_idx):
                t_idx = np.where(np.logical_and(time > self.t_inj, time < self.t_inj + self.t_window))[0]

                cur_amp = current[p_id][t_idx] - current[p_id][t_idx[0] - 1]
                max_amp = current[p_id][m_idx] - current[p_id][t_idx[0] - 1]

                if (m_idx < min_time_idx or
                        m_idx > max_time_idx or
                        abs(max_amp) < 1e-12):
                    # No peaks
                    continue

                good_max.append(max_amp * 1e9)
                peak_amp.append(max_amp * 1e9)
                peak_time.append((time[m_idx] - time[t_idx[0]]) * 1e3)
                volt_curve.append(((time[t_idx] - time[t_idx[0]]) * 1e3, cur_amp * 1e9))

            # Pick which curves to plot
            sort_idx = np.argsort(peak_amp)
            if len(sort_idx) < n_plot_max:
                keep_idx = sort_idx
            else:
                keep_idx = [sort_idx[int(np.round(x))] for x in np.linspace(0, len(sort_idx) - 1, n_plot_max)]

            for x in keep_idx:
                plt.plot(volt_curve[x][0], volt_curve[x][1], 'k-')

            plt.scatter(peak_time, peak_amp, marker=".", c="blue", s=100)

            n_type = self.data["neurons"][plot_id[0]]["type"]
            if (sim_type, n_type) in self.exp_data_dict:
                exp_data = self.exp_data_dict[(sim_type, n_type)]
                t = self.t_window * 1e3 * (1 + 0.03 * np.random.rand(exp_data.shape[0]))
                plt.scatter(t, -exp_data, marker=".", c="red", s=100)

            if self.plot_exp_trace and (sim_type, n_type) in self.exp_trace_dict:
                data = self.exp_trace_dict[(sim_type, n_type)]
                t_exp = data[:, 0]
                v_exp = data[:, 1:]
                t_idx = np.where(t_exp < self.t_window * 1e3)[0]
                plt.plot(t_exp[t_idx], v_exp[t_idx, :], c="red")

            plt.title(f"{pre_type} to {plot_type}")
            plt.xlabel("Time (ms)")
            plt.ylabel("Current (nA)")

            # Remove part of the frame
            plt.gca().spines["right"].set_visible(False)
            plt.gca().spines["top"].set_visible(False)

            plt.tight_layout()
            plt.ion()
            plt.show()
            plt.savefig(fig_name, dpi=300)

            # Also plot histogram
            plt.figure()
            plt.hist(good_max)
            plt.xlabel("Current (nA)")
            plt.title(f"{pre_type} to {plot_type}")

            # Remove part of the frame
            plt.gca().spines["right"].set_visible(False)
            plt.gca().spines["top"].set_visible(False)

            plt.tight_layout()
            plt.ion()
            plt.show()
            plt.savefig(fig_name_hist, dpi=300)

        import pdb
        pdb.set_trace()

    ############################################################################

    def set_gaba_rev(self, v_rev_cl):

        print("Setting GABA reversal potential to " + str(v_rev_cl * 1e3) + " mV")

        for s in self.snudda_sim.synapse_list:
            assert s.e == -65, "It should be GABA synapses only that we modify!"
            s.e = v_rev_cl * 1e3

    ############################################################################


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser("Simulate Chuhma 2011 and Straub2016 experiments")
    parser.add_argument("action", choices=["setup", "run", "analyse"])
    parser.add_argument("simType", help="Experiment we want to perform",
                        choices=["Chuhma2011",
                                 "Straub2016FS",
                                 "Straub2016LTS",
                                 "Szydlowski2013"])

    parser.add_argument("simName",
                        help="Simulation name, eg. networks/Chuhma2011-v1",
                        type=str)

    args = parser.parse_args()

    sim_name = args.simName
    sim_type = args.simType

    sm = SnuddaModelCurrentInjections(sim_name, sim_type)

    if args.action == "setup":
        print(f"Setup {sim_name}")
        # simName = "networks/Chuhma2011-v1"

        sm.define_network(sim_name)
        sm.create_network(sim_name)

    if args.action == "run":
        print(f"Running {sim_name}")
        sm.simulate_network(sim_name)

    if args.action == "analyse":
        print(f"Analyse {sim_name}")
        sm.analyse_network(sim_type=sim_type, sim_name=sim_name)
