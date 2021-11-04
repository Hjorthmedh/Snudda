import numpy as np
from snudda import SnuddaLoad
import matplotlib.pyplot as plt


class PlotDensitySlice:

    def __init__(self, network_file):

        self.sl = SnuddaLoad(network_file)

        self.neuron_positions = self.sl.data["neuronPositions"]
        self.neuron_type_lookup = dict()

        for nt in set([n["type"] for n in self.sl.data["neurons"]]):
            self.neuron_type_lookup[nt] = self.sl.get_neuron_id_of_type(neuron_type=nt)

    def plot_slice(self, neuron_type,
                   x_min=None, x_max=None,
                   y_min=None, y_max=None,
                   z_min=None, z_max=None,
                   projection=None):

        if projection is None:
            projection = "3d"

        n_neurons = self.sl.data["nNeurons"]

        cell_id = self.neuron_type_lookup[neuron_type]
        ok_mask = np.zeros((n_neurons,), dtype=int)
        ok_mask[cell_id] = 1

        if x_min is not None:
            ok_mask = np.logical_and(ok_mask, x_min <= self.neuron_positions[:, 0])

        if x_max is not None:
            ok_mask = np.logical_and(ok_mask, self.neuron_positions[:, 0] <= x_max)

        if y_min is not None:
            ok_mask = np.logical_and(ok_mask, y_min <= self.neuron_positions[:, 1])

        if y_max is not None:
            ok_mask = np.logical_and(ok_mask, self.neuron_positions[:, 1] <= y_max)

        if z_min is not None:
            ok_mask = np.logical_and(ok_mask, z_min <= self.neuron_positions[:, 2])

        if z_max is not None:
            ok_mask = np.logical_and(ok_mask, self.neuron_positions[:, 2] <= z_max)

        # cell_pos = self.neuron_positions[ok_idx, :].copy()
        ok_idx = np.where(ok_mask)[0]
        cell_pos = self.neuron_positions[ok_idx, :].copy()

        print(f"Plotting {np.count_nonzero(ok_mask)} {neuron_type} neurons")

        # for pos in cell_pos:
        #     print(f"xyz: {pos}")

        fig = plt.figure()

        if projection == "yz":
            plt.scatter(cell_pos[:, 1], cell_pos[:, 2])
            plt.axis("equal")
            fig.tight_layout()

        elif projection == "xz":
            plt.scatter(cell_pos[:, 0], cell_pos[:, 2])
            plt.axis("equal")
            fig.tight_layout()

        elif projection == "xy":
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
