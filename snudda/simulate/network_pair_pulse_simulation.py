# This code assumes you have created a small network of neurons, it will
# setup current injections
#
#
# OBS! We set a holding current to keep neuron at around -80mV
# We also change GABA reversal potential to -40mV, since internal Cl is 30mM
# and external Cl is 135 mM, temperature is 306K
#
#
# How to:
#
# * Edit snudda_init_custom.py to have the neurons you want.
#
# * Generate network 
#
#  python3 snudda/simulate/network_pair_pulse_simulation.py setup Planert2010 networks/Planert2010-v1
#  snudda place networks/Planert2010-v1
#  snudda detect networks/Planert2010-v1
#  snudda prune networks/Planert2010-v1

# * Figure out where to put the slcie cut (with plotOnly equation is ignored)
# 
#  python3 snudda/utils/cut.py networks/Planert2010-v1/network-synapses.hdf5 "z>0" --plotOnly
#
# * Look at networks/Planert2010-v1/network-cut-slice.hdf5.pdf to decide cut plane
#
# * Compile mod files (we now have failure rates for GABA)
#
# nrnivmodl data/neurons/mechanisms/
#
# * Cut the slice, so z > 0.00504 is kept
#
#  python3 cut.py networks/Planert2010-v1/network-synapses.hdf5 "abs(z)<100e-6"
#
# * Look at networks/Planert2010-v1/network-cut-slice.hdf5.pdf to verify cut plane
#
# * Run dSPN -> iSPN calibration (you get dSPN -> dSPN data for free then)
#
#  mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda/simulate/network_pair_pulse_simulation.py run Planert2010 networks/Planert2010-v1/network-cut-slice.hdf5 --pre dSPN --post iSPN
#
# *  Analyse
#
#  python3 snudda/simulate/network_pair_pulse_simulation.py analyse networks/Planert2010-v1/network-cut-slice.hdf5 dSPN iSPN
#  python3 snudda/simulate/network_pair_pulse_simulation.py analyse Planert2010 networks/Planert2010-v1/network-cut-slice.hdf5 --pre dSPN --post dSPN
#
# * Look at plot with traces overlayed and histogram of voltage amplitudes
# (When you do preType to postType, you also get preType to preType for free
# since both voltages are recorded

import os
import glob
import sys

import numpy as np
from snudda.simulate.simulate import SnuddaSimulate
from snudda.utils.load import SnuddaLoad
from snudda.utils import SnuddaLoadNetworkSimulation

import matplotlib
import matplotlib.pyplot as plt
import neuron


# We want to match Taverna 2008 data:

# The slices were 300 μm thick.  MSNs sampled were 25–100 μm from the
# surface of the slice to facilitate paired recordings.  Pairs always
# had cell bodies that were similar in depth (near the same focal
# plane) in an attempt to minimize the probability that the local axon
# collateral of one or the other cell was cut. After choosing a given
# pair of cells (lateral distance between somata, ≤50 μm)
#
# Assume loss of 100 micrometer in depth, at KI they start with 250 micrometers
# and get 150 micrometers after.


class SnuddaNetworkPairPulseSimulation:

    # TODO: Allow hold_voltage to be a dictionary, with neuron type as lookup
    #       to allow different holding voltages

    def __init__(self, network_path,
                 pre_type, post_type=None,
                 exp_type=None,
                 current_injection=10e-9,
                 hold_voltage=-80e-3,
                 max_dist=50e-6,
                 log_file=None):

        if os.path.isfile(network_path):
            self.network_file = network_path
            self.network_path = os.path.dirname(network_path)
        else:
            self.network_file = os.path.join(network_path, "network-synapses.hdf5")
            self.network_path = network_path

        self.exp_type = exp_type

        self.pre_type = pre_type

        if post_type:
            self.post_type = post_type
        else:
            self.post_type = "ALL"

        self.cur_inj = current_injection
        self.hold_v = hold_voltage

        if log_file:
            self.log_file = log_file
        else:
            self.log_file = os.path.join(network_path, "log", "pair-pulse.log")

        print(f"Using log file {self.log_file}")

        self.max_dist = max_dist

        print(f"Checking depolarisation/hyperpolarisation of {pre_type} to {post_type} synapses")

        self.inj_spacing = 0.5  # Tried with 0.2 before, too close
        self.inj_duration = 1e-3

        # Voltage file
        # self.volt_file = os.path.join(self.network_path,
        #                               f"synapse-calibration-volt-{self.pre_type}-{self.post_type}.txt")
        # self.volt_file_alt_mask = os.path.join(self.network_path,
        #                                        f"synapse-calibration-volt-{self.pre_type}-*.txt")

        self.snudda_sim = None   # Defined in run_sim
        self.snudda_load = None  # Defined in analyse
        self.data = None         # Defind in analyse
        self.holding_i_clamp_list = []
        self.pre_id = []
        self.possible_post_id = []
        self.inj_info = []
        self.exp_data = dict()

    ############################################################################

    def setup(self, n_dSPN=120, n_iSPN=120, n_FS=20, n_LTS=0, n_ChIN=0,
              volume_type=None,
              neuron_density=80500,
              side_len=200e-6,
              slice_depth=150e-6):

        """ Setup network for pair pulse simulation. If volume_type is 'slice', then side_len and slice_depth are used.
            If volume_type is 'cube' then side_len is used. If side_len is set to None then neuron_density is used.

        """

        from snudda.init.init import SnuddaInit

        if volume_type is None:
            volume_type = "slice"

        config_name = os.path.join(self.network_path, "network-config.json")
        cnc = SnuddaInit(struct_def={}, config_file=config_name)
        cnc.define_striatum(num_dSPN=n_dSPN, num_iSPN=n_iSPN, num_FS=n_FS, num_LTS=n_LTS, num_ChIN=n_ChIN,
                            volume_type=volume_type, side_len=side_len, slice_depth=slice_depth,
                            neuron_density=neuron_density)

        dir_name = os.path.dirname(config_name)

        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        cnc.write_json(config_name)

        print(f"\n\nsnudda place {self.network_path}")
        print(f"snudda detect {self.network_path}")
        print(f"snudda prune {self.network_path}")
        print(f"python3 snudda/utils/cut.py {self.network_path}/network-synapses.hdf5 abs(z)<100e-6")

        print("\nThe last command will pop up a figure and enter debug mode,"
              " press ctrl+D in the terminal window after inspecting the plot to continue")

        print("\n!!! Remember to compile the mod files: nrnivmodl data/neurons/mechanisms")

        print("\nTo run for example dSPN -> iSPN (and dSPN->dSPN) calibration:")
        print(f"mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_network_pair_pulse_simulation.py "
              f"run {self.exp_type} {self.network_path}/network-cut-slice.hdf5 dSPN iSPN")

        print(f"\npython3 snudda/simulate/network_pair_pulse_simulation.py analyse {self.exp_type} "
              f"{self.network_path}/network-cut-slice.hdf5 --pre dSPN --post iSPN"
              f"\npython3 snudda_network_pair_pulse_simulation.py analyse {self.network_path}/network-cut-slice.hdf5 "
              f"--pre iSPN --post dSPN")

    ############################################################################

    def setup_holding_volt(self, hold_v=None, sim_end=None):

        assert sim_end is not None, \
            "setup_holding_volt: Please set sim_end, for holding current"

        if hold_v is None:
            hold_v = self.hold_v

        if hold_v is None:
            print("Not using holding voltage, skipping.")
            return

        # Setup vClamps to calculate what holding current will be needed
        soma_v_clamp = []

        soma_list = [self.snudda_sim.neurons[x].icell.soma[0] for x in self.snudda_sim.neurons]

        for s in soma_list:
            vc = neuron.h.SEClamp(s(0.5))
            vc.rs = 1e-9
            vc.amp1 = hold_v * 1e3
            vc.dur1 = 100

            soma_v_clamp.append((s, vc))

        neuron.h.finitialize(hold_v * 1e3)
        neuron.h.tstop = 100
        neuron.h.run()

        self.holding_i_clamp_list = []

        # Setup iClamps
        for s, vc in soma_v_clamp:
            cur = float(vc.i)
            ic = neuron.h.IClamp(s(0.5))
            ic.amp = cur
            ic.dur = 2 * sim_end * 1e3
            self.holding_i_clamp_list.append(ic)

        # Remove vClamps
        v_clamps = None
        vc = None

        ############################################################################

    def set_gaba_rev(self, v_rev_cl):

        print(f"Setting GABA reversal potential to {v_rev_cl * 1e3} mV")

        for s in self.snudda_sim.synapse_list:
            assert s.e == -65, "It should be GABA synapses only that we modify!"
            s.e = v_rev_cl * 1e3

    ############################################################################

    def run_sim(self, gaba_rev, pre_id=None):

        self.snudda_sim = SnuddaSimulate(network_file=self.network_file,
                                         input_file=None,
                                         log_file=self.log_file,
                                         disable_gap_junctions=True)

        self.snudda_sim.setup()

        # A current pulse to all pre synaptic neurons, one at a time
        if pre_id:
            print(f"Using user defined pre_id: {pre_id}")
            self.pre_id = pre_id
        else:
            self.pre_id = [x["neuronID"] for x in self.snudda_sim.network_info["neurons"] if x["type"] == self.pre_type]

        # injInfo contains (preID,injStartTime)
        self.inj_info = list(zip(self.pre_id, self.inj_spacing + self.inj_spacing * np.arange(0, len(self.pre_id))))

        sim_end = self.inj_info[-1][1] + self.inj_spacing

        # Set the holding voltage
        self.setup_holding_volt(hold_v=self.hold_v, sim_end=sim_end)

        self.set_gaba_rev(gaba_rev)

        # Add current injections defined in init
        for (nid, t) in self.inj_info:
            print(f"Current injection to {nid} at {t} s")
            self.snudda_sim.add_current_injection(neuron_id=nid,
                                                  start_time=t,
                                                  end_time=t + self.inj_duration,
                                                  amplitude=self.cur_inj)

        # !!! We could maybe update code so that for postType == "ALL" we
        # record voltage from all neurons

        if self.post_type == "ALL":
            self.snudda_sim.add_recording()
        else:
            # Record from all the potential post synaptic neurons
            self.snudda_sim.add_recording_of_type(self.post_type)

            # Also save the presynaptic traces for debugging, to make sure they spike
            self.snudda_sim.add_recording_of_type(self.pre_type)

        # Run simulation
        self.snudda_sim.run(sim_end * 1e3, hold_v=self.hold_v)

        # Write results to disk
        # self.snudda_sim.write_voltage_OLD(self.volt_file)
        self.snudda_sim.write_output()

    ############################################################################

    # This extracts all the voltage deflections, to see how strong they are

    def analyse(self, max_dist=None, n_max_show=10, pre_id=None, post_type=None):

        self.setup_exp_data()

        if max_dist is None:
            max_dist = self.max_dist

        # Read the data
        self.snudda_load = SnuddaLoad(self.network_file)
        self.data = self.snudda_load.data

        ssd = SnuddaLoadNetworkSimulation(network_path=self.network_path)
        voltage, time = ssd.get_voltage()

        check_width = 0.05

        # Generate current info structure
        # A current pulse to all pre synaptic neurons, one at a time
        if pre_id:
            print(f"Using user provided pre_id = {pre_id}\n"
                  f"This must match what was used for simulation! BE CAREFUL!")
            self.pre_id = pre_id
        else:
            self.pre_id = [x["neuronID"] for x in self.data["neurons"] if x["type"] == self.pre_type]

        if post_type is None:
            post_type = self.post_type

        assert post_type != "ALL", "You need to specify a neuron type as post_type, e.g. FS"

        self.possible_post_id = [x["neuronID"] for x in self.data["neurons"] if x["type"] == post_type]

        # injInfo contains (preID,injStartTime)
        self.inj_info = zip(self.pre_id, self.inj_spacing + self.inj_spacing * np.arange(0, len(self.pre_id)))

        # For each pre synaptic neuron, find the voltage deflection in each
        # of its post synaptic neurons

        synapse_data = []
        too_far_away = 0

        for (pre_id, t) in self.inj_info:
            # Post synaptic neuron to preID
            synapses, coords = self.snudda_load.find_synapses(pre_id=pre_id)

            post_id_set = set(synapses[:, 1]).intersection(self.possible_post_id)
            pre_pos = self.snudda_load.data["neuronPositions"][pre_id, :]

            for post_id in post_id_set:

                if max_dist is not None:
                    post_pos = self.snudda_load.data["neuronPositions"][post_id, :]
                    if np.linalg.norm(pre_pos - post_pos) > max_dist:
                        too_far_away += 1
                        continue

                # There is a bit of synaptic delay, so we can take voltage
                # at first timestep as baseline
                t_idx = np.where(np.logical_and(t <= time, time <= t + check_width))[0]
                synapse_data.append((time[t_idx], voltage[post_id][t_idx]))

        if max_dist is not None:
            print(f"Number of pairs excluded, distance > {max_dist * 1e6} mum : {too_far_away}")

        # Fig names:
        trace_fig = os.path.join(os.path.dirname(self.network_file),
                                 "figures",
                                 f"{self.exp_type}-synapse-calibration-volt-traces-{self.pre_type}-{post_type}.pdf")

        hist_fig = os.path.join(os.path.dirname(self.network_file),
                                "figures",
                                f"{self.exp_type}-synapse-calibration-volt-histogram-{self.pre_type}-{post_type}.pdf")

        fig_dir = os.path.join(os.path.dirname(self.network_file), "figures")

        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)

        # Extract the amplitude of all voltage pulses
        amp = np.zeros((len(synapse_data),))
        idx_max = np.zeros((len(synapse_data),), dtype=int)
        t_max = np.zeros((len(synapse_data),))

        for i, (t, v) in enumerate(synapse_data):
            # Save the largest deflection -- with sign
            idx_max[i] = np.argmax(np.abs(v - v[0]))
            t_max[i] = t[idx_max[i]] - t[0]
            amp[i] = v[idx_max[i]] - v[0]

        assert len(amp) > 0, "No responses... too short distance!"

        print(f"Min amp: {np.min(amp)}")
        print(f"Max amp: {np.max(amp)}")
        print(f"Mean amp: {np.mean(amp)} +/- {np.std(amp)}")
        print(f"Amps: {amp}")

        # Now we have all synapse deflections in synapseData
        matplotlib.rcParams.update({'font.size': 22})

        sort_idx = np.argsort(amp)
        if len(sort_idx) > n_max_show:
            keep_idx = [sort_idx[int(np.round(x))] for x in np.linspace(0, len(sort_idx) - 1, n_max_show)]
        else:
            keep_idx = sort_idx

        plt.figure()
        for x in keep_idx:
            t, v = synapse_data[x]

            plt.plot((t - t[0]) * 1e3, (v - v[0]) * 1e3, color="black")

        plt.scatter(t_max * 1e3, amp * 1e3, color="blue", marker=".", s=100)

        if (self.exp_type, self.pre_type, post_type) in self.exp_data:
            exp_mean, exp_std = self.exp_data[(self.exp_type, self.pre_type, post_type)]

            t_end = (t[-1] - t[0]) * 1e3

            axes = plt.gca()
            ay = axes.get_ylim()
            # Plot SD or 1.96 SD?
            plt.errorbar(t_end, exp_mean, exp_std, ecolor="red",
                         marker='o', color="red")

            model_mean = np.mean(amp) * 1e3
            model_std = np.std(amp) * 1e3
            plt.errorbar(t_end - 2, model_mean, model_std, ecolor="blue",
                         marker="o", color="blue")

            axes.set_ylim(ay)

        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")
        # plt.title(str(len(synapseData)) + " traces")
        plt.title(f"{self.pre_type} to {post_type}")

        # Remove part of the frame
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)

        plt.tight_layout()
        plt.ion()
        plt.show()
        plt.savefig(trace_fig, dpi=300)

        plt.figure()
        plt.hist(amp * 1e3, bins=20)
        plt.title(f"{self.pre_type} to {post_type}")
        plt.xlabel("Voltage deflection (mV)")

        # Remove part of the frame
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)

        plt.tight_layout()
        plt.show()
        plt.savefig(hist_fig, dpi=300)

        plt.pause(10)

    ############################################################################

    def setup_exp_data(self):

        self.exp_data = dict()

        planert_d1_d1 = (0.24, 0.15)
        planert_d1_d2 = (0.33, 0.15)
        planert_d2_d1 = (0.27, 0.09)
        planert_d2_d2 = (0.45, 0.44)
        planert_fs_d1 = (4.8, 4.9)
        planert_fs_d2 = (3.1, 4.1)

        self.exp_data[("Planert2010", "dSPN", "dSPN")] = planert_d1_d1
        self.exp_data[("Planert2010", "dSPN", "iSPN")] = planert_d1_d2
        self.exp_data[("Planert2010", "iSPN", "dSPN")] = planert_d2_d1
        self.exp_data[("Planert2010", "iSPN", "iSPN")] = planert_d2_d2
        self.exp_data[("Planert2010", "FS", "dSPN")] = planert_fs_d1
        self.exp_data[("Planert2010", "FS", "iSPN")] = planert_fs_d2

    ############################################################################


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Calibrate synapse conductances")
    parser.add_argument("task", choices=["setup", "run", "analyse"])
    parser.add_argument("expType", help="Experiment we replicate", choices=["Planert2010", "Szydlowski2013"])
    parser.add_argument("networkFile", help="Network file (hdf5) or network directory")
    parser.add_argument("--preType", "--pre", help="Pre synaptic neuron type", default="dSPN")
    parser.add_argument("--postType", "--post", default="ALL",
                        help="Post synaptic neuron type (for run task, "
                             "postType can be 'ALL' to record from all neuron)")
    parser.add_argument("--maxDist", help="Only check neuron pairs within (mum)", type=float, default=None)
    args = parser.parse_args()

    if args.maxDist is None:
        max_dist = 50e-6
    elif args.maxDist == "None":
        max_dist = None
    else:
        max_dist = float(args.maxDist)

    print(f"Using maxDist = {max_dist}")

    if args.expType == "Planert2010":
        n_dSPN = 120
        n_iSPN = 120
        n_FS = 20
        n_LTS = 0
        n_ChIN = 0

        hold_v = -80e-3
        max_dist = 100e-6 if args.maxDist is None else args.maxDist
        GABA_rev = -40e-3

    elif args.expType == "Szydlowski2013":
        n_dSPN = 10
        n_iSPN = 10
        n_FS = 20
        n_LTS = 20
        n_ChIN = 0

        hold_v = -76e-3
        max_dist = 200e-6 if args.maxDist is None else args.maxDist
        GABA_rev = -39e-3

    else:
        print(f"Unknown expType = {args.expType}")
        sys.exit(-1)

    pps = SnuddaNetworkPairPulseSimulation(network_path=args.networkFile,
                                           exp_type=args.expType,
                                           pre_type=args.preType,
                                           post_type=args.postType,
                                           max_dist=max_dist,
                                           hold_voltage=hold_v)

    if args.task == "setup":
        pps.setup(args.networkFile,
                  n_dSPN=n_dSPN, n_iSPN=n_iSPN,
                  n_FS=n_FS, n_LTS=n_LTS, n_ChIN=n_ChIN)

    elif args.task == "run":
        pps.run_sim(gaba_rev=GABA_rev)

    elif args.task == "analyse":
        pps.analyse()
