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
import numpy as np
import json
from snudda.simulate import SnuddaSimulate

import neuron
neuron.h.load_file("stdrun.hoc")


class PairRecording(SnuddaSimulate):

    def __init__(self, network_path, experiment_config_file):

        super().__init__(network_path=network_path)

        self.simulate_neuron_ids = None

        self.sim_duration = None
        self.holding_i_clamp_list = []
        self.v_init_saved = dict()

        # Do stuff with experiment_config
        self.experiment_config_file = experiment_config_file
        self.experiment_config = self.read_experiment_config(experiment_config_file=experiment_config_file)

        self.output_voltage_file_name = None

        # Variables for synapse current recording
        self.synapse_currents = []
        self.record_from_pair = []

        self.parse_experiment_config()

    def read_experiment_config(self, experiment_config_file):

        with open(experiment_config_file, "r") as f:
            return json.load(f)

    def parse_experiment_config(self):

        self.set_neurons_to_simulate("prepost")

        # Setup the network given in network_config
        super().setup()

        self.sim_duration = self.experiment_config["meta"]["simulationDuration"]

        if "pair_recording_voltage_file" in self.experiment_config["meta"]:
            self.output_voltage_file_name = os.path.join(self.network_path, "simulation",
                                                         self.experiment_config["meta"]["pair_recording_voltage_file"])

        # Setup v_init for each neuron_id specified
        neuron_id, v_init = zip(*self.experiment_config["meta"]["vInit"])
        self.set_v_init(neuron_id=neuron_id, v_init=v_init)

        # Iterate over current injections
        for cur_info in self.experiment_config["currentInjection"]:
            stim_neuron_id = self.to_list(cur_info["neuronID"])
            stim_start_time = self.to_list(cur_info["start"])
            stim_end_time = self.to_list(cur_info["end"])
            stim_amplitude = self.to_list(cur_info["amplitude"])

            for nid, st, et, amp in zip(stim_neuron_id, stim_start_time, stim_end_time, stim_amplitude):
                self.add_current_injection(neuron_id=nid, start_time=st, end_time=et, amplitude=amp)

        # Add voltage recordings to neurons
        self.add_recording()

    @staticmethod
    def to_list(val):
        if type(val) != list:
            val = [val]

        return val

    def set_v_init(self, neuron_id, v_init):

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
            vc.dur1 = 100
            soma_v_clamp.append((s, vc))

        self.set_v_init_helper(neuron_id=neuron_id, v_init=v_init)
        neuron.h.finitialize()

        neuron.h.tstop = 100
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

        if selection is None or selection == "ALL":
            # No restrictions, use all
            self.simulate_neuron_ids = None

        elif selection == "prepost":

            pre_id = set(sum([c["neuronID"] if type(c["neuronID"]) == list else [c["neuronID"]]
                              for c in self.experiment_config["currentInjection"]], []))
            sim_id = pre_id

            for pid in pre_id:
                post_id = set(self.snudda_loader.find_synapses(pre_id=pid)[0][:, 1])
                sim_id = sim_id.union(post_id)

            self.simulate_neuron_ids = sorted(list(sim_id))

        else:
            assert False, f"Unsupported selection {selection}, use 'ALL' or 'prepost'"

    def distribute_neurons(self):

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

        for nid, v in self.v_init_saved.items():
            self.set_resting_voltage(neuron_id=nid, rest_volt=v)

        # Run simulation
        super().run(self.sim_duration * 1e3, hold_v=None)

        # Write results to disk
        try:
            self.write_voltage(output_file=self.output_voltage_file_name)
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
                    syn_i = neuron.h.Vector().record(syn._ref_i)
                    self.synapse_currents.append((src_id, dest_id, syn_i))

            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr, is_error=True)
                import pdb
                pdb.set_trace()

    def mark_synapses_for_recording(self, pre_neuron_id, post_neuron_id):

        self.record_from_pair.append((pre_neuron_id, post_neuron_id))

    def add_synapse_current_recording(self, pre_neuron_id, post_neuron_id):

        # 1. Find out which synapses connect pre and post
        synapses = self.snudda_loader.find_synapses(pre_id=pre_neuron_id,
                                                    post_id=post_neuron_id)

        # 2. How do we access that synapse, what type is it?
        channel_model_id = synapses[:, 6]
        segment_id = synapses[:, 9]
        segment_x = synapses[:, 10] * 1e-3

        segment_x[segment_x == 0] = 0.01
        segment_x[segment_x == 1] = 0.99

        # TODO: WE NEED TO HANDLE SYNAPSES ON THE SOMA ALSO!!

    def write_meta_data(self):