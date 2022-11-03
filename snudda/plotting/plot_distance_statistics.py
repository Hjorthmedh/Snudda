import numpy as np
from snudda.utils.load import SnuddaLoad
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt


class PlotDistanceStatistics:

    def __init__(self, network_path):

        self.snudda_load = SnuddaLoad(network_file=network_path)

    def calculate_distance_matrix(self, neuron_type, neuron_type2=None):

        if neuron_type2 is None:
            neuron_type2 = neuron_type

        neuron_id = self.snudda_load.get_neuron_id_of_type(neuron_type=neuron_type)
        neuron_id2 = self.snudda_load.get_neuron_id_of_type(neuron_type=neuron_type2)

        positions = self.snudda_load.data["neuronPositions"][neuron_id, :]
        positions2 = self.snudda_load.data["neuronPositions"][neuron_id2, :]

        dist_matrix = distance_matrix(positions, positions2)

        return dist_matrix

    def plot_distance_histogram(self, neuron_type, neuron_type2=None):

        dist_matrix = self.calculate_distance_matrix(neuron_type=neuron_type, neuron_type2=neuron_type2)

        # Nan the diagonal elements
        dist_matrix += np.diag(np.full((dist_matrix.shape[0],), np.nan))

        plt.figure()
        plt.hist(dist_matrix.flatten()*1e6, bins=30)
        plt.xlabel("Distance to neighbour ($\mu$m)")
        plt.ylabel("Count")
        if neuron_type2 is None:
            plt.title(f"Distance between {neuron_type} (mean: {np.nanmean(dist_matrix)*1e6:.1f})")
        else:
            plt.title(f"Distance between {neuron_type} and {neuron_type2} (mean: {np.nanmean(dist_matrix)*1e6:.1f})")

        # plt.ion()
        plt.show()

    def plot_distance_histogram_closest(self, neuron_type, neuron_type2=None, n_closest=5):

        dist_matrix = self.calculate_distance_matrix(neuron_type=neuron_type, neuron_type2=neuron_type2)

        # Nan the diagonal elements
        dist_matrix += np.diag(np.full((dist_matrix.shape[0],), np.nan))
        partial_dist_matrix = np.sort(dist_matrix, axis=1)[:, :n_closest]

        plt.figure()
        plt.hist(partial_dist_matrix.flatten()*1e6, bins=30)
        plt.xlabel("Distance to neighbour ($\mu$m)")
        plt.ylabel("Count")
        if neuron_type2 is None:
            plt.title(f"Distance between {neuron_type} ({n_closest} neighbours) (mean: {np.nanmean(partial_dist_matrix)*1e6:.1f})")
        else:
            plt.title(f"Distance between {neuron_type} and {neuron_type2}  ({n_closest} neighbours) (mean: {np.nanmean(partial_dist_matrix)*1e6:.1f})")

        # plt.ion()
        plt.show()


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser("Plot distance statistics")
    parser.add_argument("network_path", help="Path to neuron folder (or the network-synapse.hdf5 file)")
    parser.add_argument("neuron_type", help="Neuron type to do statistics for", type=str)
    parser.add_argument("--neuron_type2", help="Second neuron type", default=None, type=str)
    parser.add_argument("--closest", help="Only include distance to N closest neighbours", type=int)

    args = parser.parse_args()

    pds = PlotDistanceStatistics(network_path=args.network_path)

    if args.closest is None:
        pds.plot_distance_histogram(neuron_type=args.neuron_type, neuron_type2=args.neuron_type2)
    else:
        pds.plot_distance_histogram_closest(neuron_type=args.neuron_type, neuron_type2=args.neuron_type2,
                                            n_closest=args.closest)


