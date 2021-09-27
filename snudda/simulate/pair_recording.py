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
from snudda.simulate import SnuddaSimulate

import neuron
neuron.h.load_file("stdrun.hoc")


class PairRecording(SnuddaSimulate):

    def __init__(self, network_path, experiment_config_file):

        super().__init__(network_path=network_path)

        self.sim_duration = None
        self.holding_i_clamp_list = []
        self.v_init_saved = dict()

        # Do stuff with experiment_config
        self.experiment_config_file = experiment_config_file
        self.experiment_config = self.read_experiment_config(experiment_config_file=experiment_config_file)

        self.parse_experiment_config()

    def read_experiment_config(self, experiment_config_file):

        with open(experiment_config_file, "r") as f:
            return json.load(f)

    def parse_experiment_config(self):

        # Setup the network given in network_config
        super().setup()

        self.sim_duration = self.experiment_config["meta"]["simulationDuration"]

        # Setup v_init for each neuron_id specified
        neuron_id, v_init = zip(*self.experiment_config["meta"]["vInit"])
        self.set_v_init(neuron_id=neuron_id, v_init=v_init)

        # Iterate over current injections
        for cur_info in self.experiment_config["currentInjection"]:
            stim_neuron_id = self.to_list(cur_info["neuronID"])
            stim_start_time = self.to_list(cur_info["start"])
            stim_end_time = self.to_list(cur_info["end"])
            stim_amplitude = self.to_list(cur_info["amplitude"])

            self.add_current_injection(neuron_id=stim_neuron_id,
                                                  start_time=stim_start_time,
                                                  end_time=stim_end_time,
                                                  amplitude=stim_amplitude)
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

        soma_list = [self.neurons[x].icell.soma[0] for x in neuron_id]

        missing_id = [x for x in self.neurons.keys() if x not in neuron_id]

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

    def run(self):

        for nid, v in self.v_init_saved.items():
            self.set_resting_voltage(neuron_id=nid, rest_volt=v)

        # Run simulation
        super().run(self.sim_duration * 1e3, hold_v=None)

        # Write results to disk
        self.write_voltage(self.volt_file)
