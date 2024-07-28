# TODO:
#
# What cells to record from
# What cells to give current injection to


# For each cell:
#  --

#   Pre synaptic cells (list)
#   List of current injections for the pre-synaptic cells
#
#   Post synaptic cells (list)
#   List of current injections
#
#
#

import json
import os
from collections import OrderedDict

import neuron
import numpy as np

from snudda.simulate import SnuddaSimulate

neuron.h.load_file("stdrun.hoc")


class PairRecording(SnuddaSimulate):
    """ Runs simulation with current injections to one or more neurons. This code can be used to simulate
        a subset of the neurons, ie the neurons receiving current injections and their post synaptic targets.

        First create your own network using Snudda init/place/detect/prune. Then use PairRecording to setup
        and run simulation.

        For example experiment configuration files, see Snudda/snudda/data/experiment-config/pair_recording directory

    """

    def __init__(self, network_path, experiment_config_file,
                 output_file=None, network_file=None,
                 snudda_data=None,
                 disable_gap_junctions=False,
                 disable_synapses=False,
                 verbose=False):

        # Do stuff with experiment_config
        self.experiment_config_file = experiment_config_file
        self.experiment_config = self.read_experiment_config(experiment_config_file=experiment_config_file)

        if output_file is None:
            if "pair_recording_output_file" in self.experiment_config["meta"].keys():
                output_file = os.path.join(network_path, "simulation",
                                           self.experiment_config["meta"]["pair_recording_output_file"])
            else:
                output_file = os.path.join(network_path, "simulation", "output.hdf5")

        print(f"Pair recording output file: {output_file}")

        super().__init__(network_path=network_path,
                         network_file=network_file,
                         output_file=output_file,
                         disable_synapses=disable_synapses,
                         disable_gap_junctions=disable_gap_junctions,
                         snudda_data=snudda_data,
                         verbose=verbose)

        self.simulate_neuron_ids = None

        self.sim_duration = None
        self.holding_i_clamp_list = dict()
        self.v_init_saved = dict()
        self.v_hold_saved = dict()

        # Variables for synapse current recording
        # self.synapse_currents = []
        self.record_from_pair = []

        self.parse_experiment_config()

    # def __del__(self):
    #     super().__del__()

    @staticmethod
    def read_experiment_config(experiment_config_file):

        """ Loads the experimental config from JSON file. """

        with open(experiment_config_file, "r") as f:
            return json.load(f, object_pairs_hook=OrderedDict)

    def parse_experiment_config(self):

        """ Parses the experimental config file, updating the simulation as needed. """

        if "neuron_subset" in self.experiment_config["meta"]:
            neuron_subset = self.experiment_config["meta"]["neuron_subset"]
        else:
            neuron_subset = "all"  # "prepost"

        self.set_neurons_to_simulate(neuron_subset)

        if "record_synaptic_current" in self.experiment_config:
            syn_cur_record = self.experiment_config["record_synaptic_current"]

            for pre_id, post_id in syn_cur_record:
                self.write_log(f"Marking synapses for recording: {pre_id} -> {post_id}")
                self.mark_synapses_for_recording(pre_neuron_id=pre_id, post_neuron_id=post_id)

        # Setup the network given in network_config
        super().setup()

        self.sim_duration = self.experiment_config["meta"]["simulation_duration"]

        if self.output_file is None:
            if "pair_recording_output_file" in self.experiment_config["meta"]:
                self.output_file = os.path.join(self.network_path, "simulation",
                                                self.experiment_config["meta"]["pair_recording_output_file"])
                self.record.set_new_output_file(self.output_file)

        # Setup v_init for each neuron_id specified
        if "v_init" in self.experiment_config["meta"]:

            if type(self.experiment_config["meta"]["v_init"]) == list:
                neuron_id, v_init = zip(*self.experiment_config["meta"]["v_init"])
            else:
                neuron_id = None
                v_init = self.experiment_config["meta"]["v_init"]

            self.set_v_init(neuron_id=neuron_id, v_init=v_init)

        if "v_hold" in self.experiment_config["meta"]:
            if type(self.experiment_config["meta"]["v_hold"]) == list:
                neuron_id, v_hold = zip(*self.experiment_config["meta"]["v_hold"])
            else:
                neuron_id = None
                v_hold = self.experiment_config["meta"]["v_hold"]

            self.set_v_hold(neuron_id=neuron_id, v_hold=v_hold)

        if "reversal_potential" in self.experiment_config["meta"]:
            for channel_name, v_rev in self.experiment_config["meta"]["reversal_potential"]:
                self.set_channel_rev(channel_name=channel_name, v_rev=v_rev)

        # Iterate over current injections
        for cur_info in self.experiment_config["current_injection"]:
            stim_neuron_id = self.to_list(cur_info["neuron_id"])
            stim_start_time = self.to_list(cur_info["start"])
            stim_end_time = self.to_list(cur_info["end"])

            if "amp_spread" in cur_info:
                amp_spread = cur_info["amp_spread"]
            else:
                amp_spread = None

            for nid in stim_neuron_id:

                if "amplitude" in cur_info:
                    stim_amplitude = self.to_list(cur_info["amplitude"])
                elif "requested_frequency" in cur_info:
                    requested_freq = self.to_list(cur_info["requested_frequency"])
                    stim_amplitude = self.get_corresponding_current_injection(neuron_id=nid,
                                                                              frequency=requested_freq)

                    if nid in self.holding_i_clamp_list:
                        _, cur = self.holding_i_clamp_list[nid]
                        self.write_log(f"Compensating for holding current, reducing {stim_amplitude} A by {cur} A")
                        stim_amplitude -= cur

                    v_holder = self.v_hold_saved[nid]
                else:
                    assert False, f"You need to specify 'amplitude' or 'requestedFrequency' for neuron_id {stim_neuron_id}"

                self.add_current_pulses(neuron_id=nid, start_times=stim_start_time,
                                        end_times=stim_end_time, amplitudes=stim_amplitude,
                                        amp_spread=amp_spread)

                if "noise_amplitude" in cur_info:

                    if "noise_dt" in cur_info:
                        noise_dt = cur_info["noise_dt"]  # dt = how often does noise current change
                    else:
                        noise_dt = None

                    self.add_noise(neuron_id=nid, noise_std=cur_info["noise_amplitude"],
                                   duration=self.sim_duration,
                                   dt=noise_dt)

        if "record" in self.experiment_config["meta"]:
            if self.experiment_config["meta"]["record"].lower() == "soma":
                self.add_volt_recording_soma()
            elif self.experiment_config["meta"]["record"].lower() == "all":
                self.add_volt_recording_all()
            else:
                assert False, (f"Unknown record option {self.experiment_config['meta']['record']}. "
                               f"Valid options are 'soma' or 'all'")
        else:
            # Add voltage recordings to neurons
            self.add_volt_recording_soma()

        self.setup_synaptic_recordings()

    def get_corresponding_current_injection(self, neuron_id, frequency):

        """ Extracts what current injection is needed to get requested frequency from if_info.json """

        if_file = os.path.join(self.network_info["neurons"][neuron_id]["neuron_path"], "if_info.json")
        parameter_key = self.network_info["neurons"][neuron_id]["parameter_key"]
        morphology_key = self.network_info["neurons"][neuron_id]["morphology_key"]

        assert os.path.isfile(if_file), f"Unable determine current to get frequency {frequency}, needs {if_file}"

        with open(if_file, "r") as f:
            if_data = json.load(f)

        current_list = if_data[parameter_key][morphology_key]["current"]
        frequency_list = if_data[parameter_key][morphology_key]["frequency"]

        assert (np.diff(current_list) > 0).all(), f"Injected currents are not increasing in file: {if_file} " \
                                                  f"for parameter_key {parameter_key}, morphology_key {morphology_key}"

        current = np.interp(frequency, frequency_list, current_list)

        if (np.array(frequency) < frequency_list[0]).any():
            self.write_log(f"WARNING: Neuron {neuron_id} requested frequency {frequency}, "
                           f"but lowest frequency in {if_file} is {frequency_list[0]}", force_print=True)

        if (np.array(frequency) > frequency_list[-1]).any():
            self.write_log(f"WARNING: Neuron {neuron_id} requested frequency {frequency}, "
                           f"but highest frequency in {if_file} is {frequency_list[-1]}", force_print=True)

        print(f'neuron_id: {neuron_id}, current:{current}')

        return current

    @staticmethod
    def to_list(val, new_list_len=1):

        """ If val is not a list, returns a list with new_list_len elements, each with value val

        Args:
            val : variable to upgrade to list (or leave unchanged if val is a list already)
            new_list_len : length of new list created (if val is not already a list)

        """

        if type(val) not in [list, np.ndarray, range]:
            val = [val] * new_list_len

        return val

    def set_v_hold(self, v_hold, neuron_id=None):

        """ This function is currently not called by default. If you choose to use it beware that
            the holding current may change the excitability of your neuron.

            Set holding voltage v_hold for neurons speciefied by neuron_id.
            If neuron_id = None (default), all neurons get v_hold set.

            Args:
                v_hold = Holding voltage (list or int) in volt (SI-units)
                neuron_id = Neuron ID of neurons affected (list, int or None)

            """

        if neuron_id is None:
            neuron_id = self.neuron_id
            neuron_id = self.to_list(neuron_id)
            v_hold = self.to_list(v_hold, new_list_len=len(neuron_id))
        else:
            neuron_id = [x for x in neuron_id]
            v_hold = [x for x in v_hold]

        assert self.sim_duration is not None, f"setup_holding_volt: Please set self.end_time, for holding current"

        # Setup vClamps to calculate what holding current will be needed
        soma_v_clamp = []
        soma_list = [(self.neurons[x].icell.soma[0], x) for x in neuron_id if x in self.neuron_id]

        missing_id = [x for x in self.neurons.keys() if x not in neuron_id and x in self.neuron_id]

        if len(missing_id) > 0:
            print(f"Warning: v_hold is not specified for neurons {missing_id}")

        for (s, nid), vi in zip(soma_list, v_hold):
            vc = neuron.h.SEClamp(s(0.5))
            vc.rs = 1e-9
            vc.amp1 = vi * 1e3
            vc.dur1 = 200
            soma_v_clamp.append((s, vc, nid))

        self.set_v_hold_helper(neuron_id=neuron_id, v_hold=v_hold)
        neuron.h.finitialize()

        neuron.h.tstop = 200
        neuron.h.run()

        # Dont override
        # self.holding_i_clamp_list = dict()

        assert self.sim_duration is not None, ("set_v_hold: self.end_time must be set before calling, "
                                               "IClamps need to know their duration.")

        # Setup iClamps
        for s, vc, nid in soma_v_clamp:
            cur = float(vc.i)
            ic = neuron.h.IClamp(s(0.5))
            ic.amp = cur
            ic.dur = 2 * self.sim_duration * 1e3
            self.holding_i_clamp_list[nid] = ic, cur * 1e-9

        # Remove vClamps
        v_clamps = None
        vc = None

        # Set voltage also
        self.set_v_hold_helper(neuron_id=neuron_id, v_hold=v_hold)

    def set_v_hold_helper(self, neuron_id, v_hold):
        for nid, v in zip(neuron_id, v_hold):
            self.set_resting_voltage(neuron_id=nid, rest_volt=v)
            self.v_hold_saved[nid] = v
            self.v_init_saved[nid] = v

    def set_v_init(self, neuron_id, v_init):
        """ Set initial voltage v_init for neurons speciefied by neuron_id.
            If neuron_id = None (default), all neurons get v_init set.

            Args:
                v_init = Initial voltage (list or int) in volt (SI-units)
                neuron_id = Neuron ID of neurons affected (list, int or None)

            """
        if neuron_id is None:
            neuron_id = self.neuron_id

        neuron_id = self.to_list(neuron_id)
        v_init = self.to_list(v_init, new_list_len=len(neuron_id))
        for nid, v in zip(neuron_id, v_init):
            self.set_resting_voltage(neuron_id=nid, rest_volt=v)
            self.v_init_saved[nid] = v

    def set_neurons_to_simulate(self, selection=None):

        """ Sets subset of neurons to simulate. If selection "prepost" neurons receiving current injection and
            their post synaptic targets are included. If selection is "all" or None, then all neurons are included.

            If you simulate with gap junctions, you need to make sure that each neuron simulated also have all
            gap junction coupled neuron simulated.

        Args:
            selection (string or list of int) : What neurons to include in simulation ("prepost", "all"),
                                                or a list of neuron_id to simulate

        """
        if type(selection) == list:
            sim_id = selection
            self.write_log(f"Simulating user selected neuron_ids: {sim_id}")

        if selection is None or selection.lower() == "all":
            # No restrictions, use all
            self.simulate_neuron_ids = None
            self.write_log("Simulating all neurons in network.")

        elif selection.lower() == "prepost":

            pre_id = set(sum([c["neuron_id"] if type(c["neuron_id"]) == list else [c["neuron_id"]]
                              for c in self.experiment_config["current_injection"]], []))
            sim_id = pre_id

            for pid in pre_id:
                found_syn = self.snudda_loader.find_synapses(pre_id=pid)[0]
                found_gj = self.snudda_loader.find_gap_junctions(neuron_id=pid)[0]

                if found_syn is not None:
                    post_id = set(found_syn[:, 1])
                    sim_id = sim_id.union(post_id)

                if found_gj is not None:
                    gj_id = set(found_gj[:, 1])
                    sim_id = sim_id.union(gj_id)

            self.simulate_neuron_ids = sorted(list(sim_id))

            self.write_log(f"Simulating neuron IDs {sim_id}")

        else:
            assert False, f"Unsupported selection {selection}, use 'all' or 'prepost'"

    def distribute_neurons(self):

        """ Overloading the normal distribute_neurons to allow only a subset of neurons to be simulated. """

        if self.simulate_neuron_ids is None:
            super().distribute_neurons()
        else:
            self.write_log("Only simulating some of the neurons.")

            # Only include the neurons that are in simulate_neurons_id
            assert len(self.simulate_neuron_ids) >= int(self.pc.nhost()), \
                (f"Do not allocate more workers ({int(self.pc.nhost())}) than there "
                 f"are neurons included ({len(self.simulate_neuron_ids)}).")

            idx = np.arange(int(self.pc.id()), len(self.simulate_neuron_ids), int(self.pc.nhost()))
            self.neuron_id = np.array(sorted([self.simulate_neuron_ids[x] for x in idx]))

            self.neuron_id_on_node = np.zeros((self.num_neurons,), dtype=bool)
            self.neuron_id_on_node[self.neuron_id] = True

            self.write_log(f"Node {self.pc.id()} processing neurons {self.neuron_id}")

            # Clear this variable, since old value is not valid
            self.neuron_nodes = None  # TODO, set it correctly!

            # We also need to remove gap junctions that belong to neurons not simulated
            keep_gj_flag = np.ones((self.gap_junctions.shape[0], ), dtype=bool)

            for idx, gj_row in enumerate(self.gap_junctions):
                if gj_row[0] not in self.simulate_neuron_ids or gj_row[1] not in self.simulate_neuron_ids:
                    keep_gj_flag[idx] = False

            self.gap_junctions = self.gap_junctions[keep_gj_flag, :]

    def run(self):

        """ Run simulation. """

        for nid, v in self.v_init_saved.items():
            self.set_resting_voltage(neuron_id=nid, rest_volt=v)

        # Run simulation
        super().run(self.sim_duration * 1e3, hold_v=None)

        # Write results to disk
        try:
            self.record.output_file = self.output_file
            self.write_output()

        except:
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str)
            self.write_log(f"Saving failed, whoops. Entering debug mode.\n"
                           f"You might have had {self.output_file} opened elsewhere. Try closing it, then type:\n"
                           f"> self.write_output()\n"
                           f"If you are lucky, this will work and you wont loose any data.", is_error=True)
            import pdb
            pdb.set_trace()

    def setup_synaptic_recordings(self):
        for src_id, dest_id in self.record_from_pair:
            self.add_synapse_current_recording(src_id, dest_id)

    def plot_trace_overview(self, experiment_name=None):
        from snudda.plotting import PlotTraces

        if experiment_name is None:
            if "experiment_name" in self.experiment_config["meta"]:
                experiment_name = self.experiment_config["meta"]["experiment_name"]
            else:
                experiment_name = None

        pt = PlotTraces(output_file=self.output_file,
                        network_file=self.network_file,
                        experiment_name=experiment_name)

        pt.plot_traces([x for x in pt.voltage])

    def plot_traces(self, mark_current_y=-80.05e-3, trace_id=None, offset=150e-3):

        """
           Plot traces of post synaptic neuron activation, this assumes only one neuron was stimulated

           Args:
               mark_current_y (float) : Y-axis value to mark current injections at
               trace_id (list, optional) : Trace ID to show (default None, meaning all post synaptic traces)
           """

        for cur_info in self.experiment_config["current_injection"]:
            pre_id = cur_info["neuron_id"]
            cur_start = self.to_list(cur_info["start"])
            cur_end = self.to_list(cur_info["end"])
            cur_times = list(zip(cur_start, cur_end))

            skip_time = cur_start[0] / 2

            assert type(pre_id) == int, f"Plot traces assumes one pre-synaptic neuron stimulated: {pre_id}"
            if self.snudda_loader.find_synapses(pre_id=pre_id)[0] is None:
                post_id = []
            else:
                post_id = set(self.snudda_loader.find_synapses(pre_id=pre_id)[0][:, 1])

            experiment_name = self.get_experiment_name(empty_is_none=False)

            for pid in post_id:

                if trace_id is not None and pid not in trace_id:
                    print(f"Skipping trace {pid}, not in trace_id={trace_id}")
                    continue

                fig_name = f"Current-injection-{experiment_name}-pre-{pre_id}-post-{pid}.pdf"

                self.plot_trace(pre_id=pre_id, post_id=pid, fig_name=fig_name,
                                mark_current=cur_times, mark_current_y=mark_current_y,
                                skip_time=skip_time, offset=offset)

    def get_experiment_name(self, empty_is_none=True):
        if "experiment_name" in self.experiment_config["meta"]:
            experiment_name = self.experiment_config["meta"]["experiment_name"]
        elif empty_is_none:
            experiment_name = None
        else:
            experiment_name = ""

        return experiment_name

    def plot_trace(self, pre_id, post_id, offset=0, title=None, fig_name=None, skip_time=0,
                   mark_current=None, mark_current_y=None):

        from snudda.plotting import PlotTraces

        if not title:
            synapses, _ = self.snudda_loader.find_synapses(pre_id=pre_id, post_id=post_id)
            if synapses is not None:
                n_synapses = synapses.shape[0]
            else:
                n_synapses = 0

            title = f"{self.neurons[pre_id].name} ({pre_id}) -> " \
                    f"{self.neurons[post_id].name} ({post_id}) ({n_synapses} synapses)"

        experiment_name = self.get_experiment_name()

        pt = PlotTraces(output_file=self.output_file, network_file=self.network_file, experiment_name=experiment_name)
        pt.plot_traces(trace_id=post_id, offset=offset, title=title, fig_name=fig_name, skip_time=skip_time,
                       mark_current=mark_current, mark_current_y=mark_current_y)

    def plot_synaptic_currents(self, post_id):

        from snudda.plotting import PlotTraces

        experiment_name = self.get_experiment_name()
        pt = PlotTraces(output_file=self.output_file, network_file=self.network_file, experiment_name=experiment_name)
        pt.plot_synaptic_currents(post_id=post_id)

    def mark_synapses_for_recording(self, pre_neuron_id, post_neuron_id):

        """ What neuron pair synapses should we record current from? Add them to self.record_from_pair

            Args:
                pre_neuron_id : Neuron ID of presynaptic neuron
                post_neuron_id : Neuron ID of postsynaptic neuron

        """

        self.record_from_pair.append((pre_neuron_id, post_neuron_id))

    def set_channel_rev(self, channel_name, v_rev):

        print(f"Setting {channel_name} reversal potential to {v_rev * 1e3} mV")

        for syn in self.synapse_dict.values():
            if channel_name == syn.hname().split("[")[0]:
                syn.e = v_rev * 1e3


if __name__ == "__main__":
    import sys
    if '-python' in sys.argv:
        print("Network_simulate.py called through nrniv, fixing arguments")
        pythonidx = sys.argv.index('-python')
        if len(sys.argv) > pythonidx:
            sys.argv = sys.argv[pythonidx + 1:]

    from argparse import ArgumentParser, RawTextHelpFormatter
    parser = ArgumentParser("Pair recording", formatter_class=RawTextHelpFormatter)
    parser.add_argument("network_path")
    parser.add_argument("experiment_config_file")
    args = parser.parse_args()
    pr = PairRecording(network_path=args.network_path, experiment_config_file=args.experiment_config_file)
    pr.run()
