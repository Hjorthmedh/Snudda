from snudda.synaptic_fitting.run_synapse_run import RunSynapseRun
from snudda.simulate import SnuddaSimulate
import timeit
import neuron
import numpy as np

class SingleCellNetwork(RunSynapseRun):

    def __init__(self, neuron_morphology,
                 neuron_mechanisms,
                 neuron_parameters,
                 neuron_modulation,
                 stim_times,
                 synapse_density,
                 num_synapses,
                 synapse_section_id=None,  # if given, nSynapses is ignored
                 synapse_section_x=None,  # same # of elements as synapseSectionID
                 neuron_parameter_id=0,  # Which param set in parameter file to use
                 neuron_modulation_id=0,
                 holding_voltage=-70e-3,
                 holding_current=None,
                 synapse_type='glut',
                 params={},
                 time=2,
                 random_seed=None,
                 log_file=None,
                 verbose=True):

        super(SingleCellNetwork, self).__init__(
                 neuron_morphology,
                 neuron_mechanisms,
                 neuron_parameters,
                 neuron_modulation,
                 stim_times,
                 synapse_density,
                 num_synapses,
                 synapse_section_id, # if given, nSynapses is ignored
                 synapse_section_x,  # same # of elements as synapseSectionID
                 neuron_parameter_id,  # Which param set in parameter file to use
                 neuron_modulation_id,
                 holding_voltage,
                 holding_current,
                 synapse_type,
                 params,
                 time,
                 random_seed,
                 log_file,
                 verbose)

    def run(self, time=None, cond=1e-8):

        if time is None:
            time = self.time
        else:
            self.time = time

        neuron.h.finitialize(self.holding_voltage * 1e3)
        for ncs in self.nc_syn:
            ncs.weight[0] = cond

        self.set_resting_voltage(self.holding_voltage * 1e3)

        neuron.h.v_init = self.holding_voltage * 1e3
        neuron.h.tstop = time * 1e3
        self.write_log("About to start NEURON... stay safe")
        neuron.h.run()
        self.write_log("NEURON actually completed?!")

        # Convert results back to SI units
        return (np.array(self.t_save) * 1e-3,
                np.array(self.v_save) * 1e-3,
                np.array(self.i_save) * 1e-9)
