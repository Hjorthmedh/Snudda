import numpy as np
from snudda import SnuddaLoad
from matplotlib.pyplot import plt


class PlotDensitySlice:

    def __init__(self, network_file):

        sl = SnuddaLoad(network_file)

        self.neuron_positions = sl.data["neuronPositions"]
        self.neuron_type_lookup = dict()

        for nt in set([n["type"] for n in sl.data["neurons"]]):
            self.neuron_type_lookup[nt] = self.sl.get_neuron_id_of_type(neuron_type=nt)

    def plot_slice(self, neuron_type,
                   x_min=None, x_max=None,
                   y_min=None, y_max=None,
                   z_min=None, z_max=None):

        cell_id = self.neuron_type_lookup[neuron_type]
        ok_idx = np.zeros((len(cell_id),), dtype=int)
        ok_idx[cell_id] = 1

        if x_min is not None:
            ok_idx = np.logical_and(ok_idx, x_min < self.neuron_positions[:, 0])

        if x_max is not None:
            ok_idx = np.logical_and(ok_idx, self.neuron_positions[:, 0] < x_max)

        if y_min is not None:
            ok_idx = np.logical_and(ok_idx, y_min < self.neuron_positions[:, 1])

        if y_max is not None:
            ok_idx = np.logical_and(ok_idx, self.neuron_positions[:, 1] < y_max)

        if z_min is not None:
            ok_idx = np.logical_and(ok_idx, z_min < self.neuron_positions[:, 2])

        if z_max is not None:
            ok_idx = np.logical_and(ok_idx, self.neuron_positions[:, 2] < z_max)

        cell_pos = self.neuron_positions[ok_idx, :]

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        plt.scatter(x=cell_pos[:, 0], y =cell_pos[:, 1], z=cell_pos[:, 2])

        plt.show()