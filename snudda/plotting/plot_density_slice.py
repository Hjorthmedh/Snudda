import numpy as np
from snudda import SnuddaLoad
import matplotlib.pyplot as plt


class PlotDensitySlice:

    def __init__(self, network_file):

        self.sl = SnuddaLoad(network_file)

        self.neuron_positions = self.sl.data["neuronPositions"]
        self.neuron_type_lookup = dict()

        for nt in set([n["type"] for n in self.sl.data["neurons"]]):
            self.neuron_type_lookup[nt] = self.sl.get_cell_id_of_type(neuron_type=nt)

    def plot_slice(self, neuron_type,
                   x_min=None, x_max=None,
                   y_min=None, y_max=None,
                   z_min=None, z_max=None):

        n_neurons = self.sl.data["nNeurons"]

        cell_id = self.neuron_type_lookup[neuron_type]
        ok_idx = np.zeros((n_neurons,), dtype=int)
        ok_idx[cell_id] = 1

        if x_min is not None:
            ok_idx = np.logical_and(ok_idx, x_min <= self.neuron_positions[:, 0])

        if x_max is not None:
            ok_idx = np.logical_and(ok_idx, self.neuron_positions[:, 0] <= x_max)

        if y_min is not None:
            ok_idx = np.logical_and(ok_idx, y_min <= self.neuron_positions[:, 1])

        if y_max is not None:
            ok_idx = np.logical_and(ok_idx, self.neuron_positions[:, 1] <= y_max)

        if z_min is not None:
            ok_idx = np.logical_and(ok_idx, z_min <= self.neuron_positions[:, 2])

        if z_max is not None:
            ok_idx = np.logical_and(ok_idx, self.neuron_positions[:, 2] <= z_max)

        cell_pos = self.neuron_positions[ok_idx, :].copy()

        #import pdb
        #pdb.set_trace()

        fig = plt.figure()

        if (x_min is not None or x_max is not None) \
                and y_min is None and y_max is None and z_min is None and z_max is None:
            plt.scatter(cell_pos[:, 1], cell_pos[:, 2])
            plt.axis("equal")
        elif (y_min is not None or y_max is not None) \
                and x_min is None and x_max is None and z_min is None and z_max is None:
            plt.scatter(cell_pos[:, 0], cell_pos[:, 2])
            plt.axis("equal")

        elif (z_min is not None or z_max is not None) \
                and y_min is None and y_max is None and x_min is None and x_max is None:
            plt.scatter(cell_pos[:, 0], cell_pos[:, 1])
            plt.axis("equal")

        else:
            ax = fig.add_subplot(projection='3d')
            ax.scatter(xs=cell_pos[:, 0], ys=cell_pos[:, 1], zs=cell_pos[:, 2])
            # elev, azim = 0, np.pi/2
            # ax.view_init(elev, azim)
            plt.ion()

        fig.tight_layout()
        plt.show()
