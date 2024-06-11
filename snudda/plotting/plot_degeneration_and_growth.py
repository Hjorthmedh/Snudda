import os
import numpy as np

from snudda import SnuddaLoad
import matplotlib.pyplot as plt


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

    def plot_synapses(self):

        added_synapse_coords, removed_synapse_coords, kept_synapse_coords = self.categorising_synapses()

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        if removed_synapse_coords.size > 0:
            ax.scatter(removed_synapse_coords[:, 0],
                       removed_synapse_coords[:, 1],
                       removed_synapse_coords[:, 2],
                       marker='.', color="grey")

        if kept_synapse_coords.size > 0:
            ax.scatter(kept_synapse_coords[:, 0],
                       kept_synapse_coords[:, 1],
                       kept_synapse_coords[:, 2],
                       marker=".", color="black", s=50)

        if added_synapse_coords.size > 0:
            ax.scatter(added_synapse_coords[:, 0],
                       added_synapse_coords[:, 1],
                       added_synapse_coords[:, 2],
                       marker=".", color="red")

        plt.ion()
        plt.show()

        # import pdb
        # pdb.set_trace()

        return fig, ax

    def check_same_network(self):

        assert len(self.original_file.data["neurons"]) == len(self.degenerated_file.data["neurons"]), \
            f"NOT THE SAME NETWORK!! Different number of neurons in the two networks, aborting."

        for orig_neuron, degen_neuron in zip(self.original_file.data["neurons"], self.degenerated_file.data["neurons"]):
            for var in ["name", "position", "rotation"]:
                comp = orig_neuron[var] == degen_neuron[var]
                assert (type(comp) == bool and orig_neuron[var] == degen_neuron[var]) \
                    or (not type(comp) == bool and (orig_neuron[var] == degen_neuron[var]).all()), \
                    f"NOT SAME NETWORK!! Mismatch in {var} for {orig_neuron} and {degen_neuron}"

    def categorising_synapses(self):

        orig_pre_synapses, orig_pre_coords = self.original_file.find_synapses(pre_id=self.neuron_id)
        orig_post_synapses, orig_post_coords = self.original_file.find_synapses(post_id=self.neuron_id)

        degen_pre_synapses, degen_pre_coords = self.degenerated_file.find_synapses(pre_id=self.neuron_id)
        degen_post_synapses, degen_post_coords = self.degenerated_file.find_synapses(post_id=self.neuron_id)

        added_coords, removed_coords, kept_coords = self.get_synapse_locations(orig_post_coords, degen_post_coords)

        return added_coords, removed_coords, kept_coords


#        removed_synapses = np.stack(removed_synapses_list)
#        added_synapses = np.stack(added_synapses_list)
#        kept_synapses_A = np.stack([x for x, y in kept_synapses_list])
#        kept_synapses_B = np.stack([y for x, y in kept_synapses_list])
#
#        removed_synapse_coords = removed_synapses[:, 2:5] * self.original_file.data["voxelSize"] + self.original_file.data["simulationOrigo"]
#        kept_synapse_coords_A = kept_synapses_A[:, 2:5] * self.original_file.data["voxelSize"] + self.original_file.data["simulationOrigo"]
#        kept_synapse_coords_B = kept_synapses_B[:, 2:5] * self.degenerated_file.data["voxelSize"] + self.degenerated_file.data["simulationOrigo"]
#        added_synapse_coords = removed_synapses[:, 2:5] * self.degenerated_file.data["voxelSize"] + self.degenerated_file.data["simulationOrigo"]
#
#        return removed_synapses, added_synapses, kept_synapses_A, kept_synapses_B, \
#                removed_synapse_coords, added_synapse_coords, kept_synapse_coords_A, kept_synapse_coords_B


    def get_synapse_locations(self, coords_1, coords_2):

        if coords_1 is None:
            coords_1 = np.zeros((3, 0))

        if coords_2 is None:
            coords_2 = np.zeros((3, 0))

        voxel_size = self.original_file.data["voxel_size"]

        rcoords_1 = np.round(coords_1/voxel_size, decimals=0).astype(int)
        rcoords_2 = np.round(coords_2/voxel_size, decimals=0).astype(int)

        coord_set_1 = set([tuple(x) for x in rcoords_1])
        coord_set_2 = set([tuple(x) for x in rcoords_2])

        removed = coord_set_1.difference(coord_set_2)
        added = coord_set_2.difference(coord_set_1)
        kept = coord_set_2.intersection(coord_set_1)

        print(f"Added: {len(added)}, removed: {len(removed)}, kept: {len(kept)}")

        added_coords = np.array(list(added), dtype=float)*voxel_size
        removed_coords = np.array(list(removed), dtype=float)*voxel_size
        kept_coords = np.array(list(kept), dtype=float)*voxel_size

        # import pdb
        # pdb.set_trace()

        return added_coords, removed_coords, kept_coords

    def categorise_helper(self, orig_synapses, orig_coords, degen_synapses, degen_coords):

        done = False

        removed_synapses = []
        added_synapses = []
        kept_synapses = []

        orig_iter = self.iter_synapses(orig_synapses, orig_coords)
        degen_iter = self.iter_synapses(degen_synapses, degen_coords)

        orig_syn, orig_prio = next(orig_iter, (None, None))
        degen_syn, degen_prio = next(degen_iter, (None, None))

        while not done:
            if orig_prio == degen_prio:
                kept_synapses.append((orig_syn, degen_syn))
                orig_syn, orig_prio = next(orig_iter, (None, None))
                degen_syn, degen_prio = next(degen_iter, (None, None))

            elif degen_syn is None or orig_prio < degen_prio:
                removed_synapses.append(orig_syn)
                orig_syn, orig_prio = next(orig_iter, (None, None))

            elif orig_syn is None or degen_prio < orig_prio:
                added_synapses.append(degen_syn)
                degen_syn, degen_prio = next(degen_iter, (None, None))

            if orig_syn is None and degen_syn is None:
                done = True

        return removed_synapses, added_synapses, kept_synapses

    def iter_synapses(self, synapses, coords):
        for syn_row in synapses:
            yield syn_row, self.calculate_order(syn_row)

    def calculate_order(self, synapse_row):
        if synapse_row is None:
            return None
        else:
            return (synapse_row[1] * self.num_neurons + synapse_row[0])*max(100, self.num_neurons) + synapse_row[6]

