import numpy as np
from snudda.input.time_varying_input import TimeVaryingInput


class VirtualInput:

    def __init__(self, spike_file, mapping_file):

        self.spike_file = spike_file
        self.mapping_file = mapping_file

        self.data = dict()

    def add_input(self, neuron_id, spike_times):

        self.data[neuron_id] = spike_times

    def write_data(self):

        data = []
        mapping = []

        for idx, (neuron_id, spike_data) in enumerate(self.data.items()):
            mapping.append([neuron_id, idx])
            data.append(spike_data)

        np.savetxt(self.mapping_file, X=np.array(mapping, dtype=int), fmt="%i")

        with open(self.spike_file, "wt") as f:
            for row in data:
                s = " ".join([f"{x:.6f}" for x in row])
                f.write(f"{s}\n")

    def poisson_spikes(self, frequency, max_time, rng=None):

        if rng is None:
            rng = np.random.default_rng()

        spike_times = TimeVaryingInput._poisson_helper(end_time=max_time*frequency, rng=rng) / frequency

        return spike_times
