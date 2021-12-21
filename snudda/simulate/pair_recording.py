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

import os
from collections import OrderedDict

import numpy as np
import json
from snudda.simulate import SnuddaSimulate
from snudda.utils.save_network_activity import SnuddaSaveNetworkActivity

import neuron
neuron.h.load_file("stdrun.hoc")


class PairRecording(SnuddaSimulate):

    """ Runs simulation with current injections to one or more neurons. This code can be used to simulate
        a subset of the neurons, ie the neurons receiving current injections and their post synaptic targets.

        First create your own network using Snudda init/place/detect/prune. Then use PairRecording to setup
        and run simulation.

        For example experiment configuration files, see Snudda/snudda/data/experiment-config/ directory

    """

    def __init__(self, network_path, experiment_config_file):

        super().__init__(network_path=network_path)

        self.simulate_neuron_ids = None

        self.sim_duration = None
        self.holding_i_clamp_list = []
        self.v_init_saved = dict()

        # Do stuff with experiment_config
        self.experiment_config_file = experiment_config_file
        self.experiment_config = self.read_experiment_config(experiment_config_file=experiment_config_file)

        self.output_file_name = None

        # Variables for synapse current recording
        self.synapse_currents = []
        self.record_from_pair = []

        self.parse_experiment_config()

    @staticmethod
    def read_experiment_config(experiment_config_file):

        """ Loads the experimental config from JSON file. """

        with open(experiment_config_file, "r") as f:
            return json.load(f, object_pairs_hook=OrderedDict)

    def parse_experiment_config(self):

        """ Parses the experimental config file, updating the simulation as needed. """

        if "neuronSubset" in self.experiment_config["meta"]:
            neuron_subset = self.experiment_config["meta"]["neuronSubset"]
        else:
            neuron_subset = "prepost"

        self.set_neurons_to_simulate(neuron_subset)

        if "recordSynapticCurrent" in self.experiment_config:
            syn_cur_record = self.experiment_config["recordSynapticCurrent"]

            for pre_id, post_id in syn_cur_record:
                self.write_log(f"Marking synapses for recording: {pre_id} -> {post_id}")
                self.mark_synapses_for_recording(pre_neuron_id=pre_id, post_neuron_id=post_id)

        # Setup the network given in network_config
        super().setup()

        self.sim_duration = self.experiment_config["meta"]["simulationDuration"]

        if "pairRecordingOutputFile" in self.experiment_config["meta"]:
            self.output_file_name = os.path.join(self.network_path, "simulation",
                                                 self.experiment_config["meta"]["pairRecordingOutputFile"])

        # Setup v_init for each neuron_id specified
        if "vInit" in self.experiment_config["meta"]:

            if type(self.experiment_config["meta"]["vInit"]) == list:
                neuron_id, v_init = zip(*self.experiment_config["meta"]["vInit"])
            else:
                neuron_id = None
                v_init = self.experiment_config["meta"]["vInit"]

            self.set_v_init(neuron_id=neuron_id, v_init=v_init)

        if "reversal_potential" in self.experiment_config["meta"]:
            for channel_name, v_rev in self.experiment_config["meta"]["reversal_potential"]:
                self.set_channel_rev(channel_name=channel_name, v_rev=v_rev)

        # Iterate over current injections
        for cur_info in self.experiment_config["currentInjection"]:
            stim_neuron_id = self.to_list(cur_info["neuronID"])
            stim_start_time = self.to_list(cur_info["start"])
            stim_end_time = self.to_list(cur_info["end"])
            stim_amplitude = self.to_list(cur_info["amplitude"])

            for nid in stim_neuron_id:
                self.add_current_pulses(neuron_id=nid, start_times=stim_start_time,
                                        end_times=stim_end_time, amplitudes=stim_amplitude)

        # Add voltage recordings to neurons
        self.add_recording()

    @staticmethod
    def to_list(val, new_list_len=1):

        """ If val is not a list, returns a list with new_list_len elements, each with value val

        Args:
            val : variable to upgrade to list (or leave unchanged if val is a list already)
            new_list_len : length of new list created (if val is not already a list)

        """

        if type(val) not in [list, np.ndarray, range]:
            val = [val]*new_list_len

        return val

    def set_v_init(self, v_init, neuron_id=None):

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

        assert self.sim_duration is not None, \
            f"setup_holding_volt: Please set self.end_time, for holding current"

        # Setup vClamps to calculate what holding current will be needed
        soma_v_clamp = []
        soma_list = [self.neurons[x].icell.soma[0] for x in neuron_id if x in self.neuron_id]

        missing_id = [x for x in self.neurons.keys() if x not in neuron_id and x in self.neuron_id]

        if len(missing_id) > 0:
            print(f"Warning: v_init is not specified for neurons {missing_id}")

        for s, vi in zip(soma_list, v_init):
            vc = neuron.h.SEClamp(s(0.5))
            vc.rs = 1e-9
            vc.amp1 = vi * 1e3
            vc.dur1 = 200
            soma_v_clamp.append((s, vc))

        self.set_v_init_helper(neuron_id=neuron_id, v_init=v_init)
        neuron.h.finitialize()

        neuron.h.tstop = 200
        neuron.h.run()

        self.holding_i_clamp_list = []

        assert self.sim_duration is not None, ("set_v_init: self.end_time must be set before calling, "
                                               "IClamps need to know their duration.")

        # Setup iClamps
        for s, vc in soma_v_clamp:
            cur = float(vc.i)
            ic = neuron.h.IClamp(s(0.5))
            ic.amp = cur
            ic.dur = 2 * self.sim_duration * 1e3
            self.holding_i_clamp_list.append(ic)

        # Remove vClamps
        v_clamps = None
        vc = None

        # Set voltage also
        self.set_v_init_helper(neuron_id=neuron_id, v_init=v_init)

    def set_v_init_helper(self, neuron_id, v_init):
        for nid, v in zip(neuron_id, v_init):
            self.set_resting_voltage(neuron_id=nid, rest_volt=v)
            self.v_init_saved[nid] = v

    def set_neurons_to_simulate(self, selection=None):

        """ Sets subset of neurons to simulate. If selection "prepost" neurons receiving current injection and
            their post synaptic targets are included. If selection is "all" or None, then all neurons are included.

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

            pre_id = set(sum([c["neuronID"] if type(c["neuronID"]) == list else [c["neuronID"]]
                              for c in self.experiment_config["currentInjection"]], []))
            sim_id = pre_id

            for pid in pre_id:
                post_id = set(self.snudda_loader.find_synapses(pre_id=pid)[0][:, 1])
                sim_id = sim_id.union(post_id)

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
            self.neuron_id = np.array([self.simulate_neuron_ids[x] for x in idx])

            self.write_log(f"Node {self.pc.id()} processing neurons {self.neuron_id}")

            # Clear this variable, since old value is not valid
            self.neuron_nodes = None  # TODO, set it correctly!

    def run(self):

        """ Run simulation. """

        for nid, v in self.v_init_saved.items():
            self.set_resting_voltage(neuron_id=nid, rest_volt=v)

        # Run simulation
        super().run(self.sim_duration * 1e3, hold_v=None)

        # Write results to disk
        try:
            save = SnuddaSaveNetworkActivity(output_file=self.output_file_name)
            save.write(t_save=self.t_save, v_save=self.v_save, v_key=self.v_key,
                       t_spikes=self.t_spikes, id_spikes=self.id_spikes)

            pre_id = np.array([x[0] for x in self.synapse_currents])
            post_id = np.array([x[1] for x in self.synapse_currents])
            cur = [np.array(x[2]) for x in self.synapse_currents]
            save.write_currents(t_save=self.t_save, i_save=cur, pre_id=pre_id, post_id=post_id)

        except:
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str)
            self.write_log(f"Saving failed, whoops. Entering debug mode.")
            import pdb
            pdb.set_trace()

    def connect_neuron_synapses(self, start_row, end_row):

        """ Connects the synapses present in the synapse matrix between start_row and end_row-1. """

        source_id_list, dest_id, dend_sections, sec_id, sec_x, synapse_type_id, \
            axon_distance, conductance, parameter_id = self.get_synapse_info(start_row=start_row, end_row=end_row)

        for (src_id, section, section_x, s_type_id, axon_dist, cond, p_id) \
                in zip(source_id_list, dend_sections, sec_x, synapse_type_id,
                       axon_distance, conductance, parameter_id):

            try:
                syn = self.add_synapse(cell_id_source=src_id,
                                       dend_compartment=section,
                                       section_dist=section_x,
                                       synapse_type_id=s_type_id,
                                       axon_dist=axon_dist,
                                       conductance=cond,
                                       parameter_id=p_id)

                if (src_id, dest_id) in self.record_from_pair:
                    self.write_log(f"Adding synaptic recording {src_id} -> {dest_id} (NEURON)")
                    syn_i = neuron.h.Vector().record(syn._ref_i)
                    self.synapse_currents.append((src_id, dest_id, syn_i))

            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr, is_error=True)
                import pdb
                pdb.set_trace()

    def plot_trace_overview(self):
        from snudda.plotting import PlotTraces
        pt = PlotTraces(output_file=self.output_file_name,
                        network_file=self.network_file)

        pt.plot_traces([x for x in pt.voltage])

    def plot_traces(self, mark_current_y=-80.05e-3):

        for cur_info in self.experiment_config["currentInjection"]:
            pre_id = cur_info["neuronID"]
            cur_start = self.to_list(cur_info["start"])
            cur_end = self.to_list(cur_info["end"])
            cur_times = list(zip(cur_start, cur_end))

            skip_time = cur_start[0]/2

            assert type(pre_id) == int, f"Plot traces assumes one pre-synaptic neuron stimulated: {pre_id}"

            post_id = set(self.snudda_loader.find_synapses(pre_id=pre_id)[0][:, 1])

            for pid in post_id:
                fig_name = f"Current-injection-pre-{pre_id}-post-{pid}.pdf"
                self.plot_trace(pre_id=pre_id, post_id=pid, fig_name=fig_name,
                                mark_current=cur_times, mark_current_y=mark_current_y,
                                skip_time=skip_time)

    def plot_trace(self, pre_id, post_id, offset=0, title=None, fig_name=None, skip_time=0,
                   mark_current=None, mark_current_y=None):

        from snudda.plotting import PlotTraces

        if not title:
            n_synapses = self.snudda_loader.find_synapses(pre_id=pre_id, post_id=post_id)[0].shape[0]
            title = f"{self.neurons[pre_id].name} -> {self.neurons[post_id].name} ({n_synapses} synapses)"

        pt = PlotTraces(output_file=self.output_file_name, network_file=self.network_file)
        fig = pt.plot_traces(trace_id=post_id, offset=offset, title=title, fig_name=fig_name, skip_time=skip_time)

        if mark_current:
            import matplotlib.pyplot as plt
            for t_start, t_end in mark_current:
                plt.plot([t_start-skip_time, t_end-skip_time], [mark_current_y, mark_current_y], 'r-', linewidth=5)

    def mark_synapses_for_recording(self, pre_neuron_id, post_neuron_id):

        """ What neuron pair synapses should we record current from? Add them to self.record_from_pair

            Args:
                pre_neuron_id : Neuron ID of presynaptic neuron
                post_neuron_id : Neuron ID of postsynaptic neuron

        """

        self.record_from_pair.append((pre_neuron_id, post_neuron_id))

    def set_channel_rev(self, channel_name, v_rev):

        print(f"Setting {channel_name} reversal potential to {v_rev * 1e3} mV")

        for syn in self.synapse_list:
            if channel_name == syn.hname().split("[")[0]:
                syn.e = v_rev * 1e3

