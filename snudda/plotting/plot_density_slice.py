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
                   z_min=None, z_max=None,
                   projection=None):

        if projection is None:
            projection = "3d"

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

        fig = plt.figure()

        if projection == "yz":
            plt.scatter(cell_pos[:, 1], cell_pos[:, 2])
            plt.axis("equal")
            fig.tight_layout()

        elif projection == "xz":
            plt.scatter(cell_pos[:, 0], cell_pos[:, 2])
            plt.axis("equal")
            fig.tight_layout()

        elif projection == "xz":
            plt.scatter(cell_pos[:, 0], cell_pos[:, 1])
            plt.axis("equal")
            fig.tight_layout()

        elif projection == "3d":
            ax = fig.add_subplot(projection='3d')
            ax.scatter(xs=cell_pos[:, 0], ys=cell_pos[:, 1], zs=cell_pos[:, 2])
            # elev, azim = 0, np.pi/2
            # ax.view_init(elev, azim)
            plt.ion()
        else:
            print(f"Unknown projection: {projection} (use 'xy', 'xz', 'yz' or '3d')")
            return

        plt.show()
