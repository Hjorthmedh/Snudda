import os
import numpy as np

from snudda import SnuddaLoad

# TODO: We now have kept, removed and added synapses.
#       Plot the neuron morphology of the original and the degenerated neuron
#       plot synapses depending on type

class PlotDegenerationAndGrowth:

    def __init__(self, original_network_path, degenerated_network_path, neuron_id):

        self.neuron_id = neuron_id
        self.original_file = SnuddaLoad(original_network_path)
        self.degenerated_file = SnuddaLoad(degenerated_network_path)

        self.check_same_network()

        self.fig_path = os.path.join(degenerated_network_path, "figures")
        if not os.path.exists(self.fig_path):
            os.mkdir(self.fig_path)

    def check_same_network(self):

        assert len(self.original_file.data["neurons"]) == len(self.degenerated_file.data["neurons"]), \
            f"NOT THE SAME NETWORK!! Different number of neurons in the two networks, aborting."

        for orig_neuron, degen_neuron in zip(self.original_file.data["neurons"], self.degenerated_file.data["neurons"]):
            for var in ["name", "position", "rotation"]:
                assert orig_neuron[var] == degen_neuron[var], \
                    f"NOT SAME NETWORK!! Mismatch in {var} for {orig_neuron} and {degen_neuron}"

    def categorising_synapses(self):

        orig_pre_synapses, orig_pre_coords = self.original_file.find_synapses(pre_id=self.neuron_id)
        orig_post_synapses, orig_post_coords = self.original_file.find_synapses(post_id=self.neuron_id)

        degen_pre_synapses, degen_pre_coords = self.degenerated_file.find_synapses(pre_id=self.neuron_id)
        degen_post_synapses, degen_post_coords = self.degenerated_file.find_synapses(post_id=self.neuron_id)

        self.num_neurons = self.original_file.data["nNeurons"]

        done = False
        orig_idx = 0
        degen_idx = 0
        orig_priority =
        degen_priority =

        while not done:


    def categorise_helper(self, orig_synapses, degen_synapses):

        done = False
        orig_idx = 0
        degen_idx = 0

        removed_synapses = []
        added_synapses = []
        kept_synapses = []

        orig_iter = self.iter_synapses(orig_synapses)
        degen_iter = self.iter_synapses(degen_synapses)

        orig_syn, orig_prio = next(orig_iter)
        degen_syn, degen_prio = next(degen_iter)

        while not done:
            if orig_prio == degen_prio:
                kept_synapses.append((orig_syn, degen_syn))
                orig_syn, orig_prio = next(orig_iter)
                degen_syn, degen_prio = next(degen_iter)

            elif degen_syn is None or orig_prio < degen_prio:
                removed_synapses.append(orig_syn)
                orig_syn, orig_prio = next(orig_iter)

            elif orig_syn is None or degen_prio < orig_prio:
                added_synapses.append(degen_synapses)
                degen_syn, degen_prio = next(degen_iter)

        return removed_synapses, added_synapses, kept_synapses

    def iter_synapses(self, synapses):
        for syn_row in synapses:
            yield syn_row, self.calculate_order(syn_row)

    def calculate_order(self, synapse_row):
        if synapse_row is None:
            return None
        else:
            return (synapse_row[1] * self.num_neurons + synapse_row[0])*max(100, self.num_neurons) + synapse_row[6]

