import numpy as np
from snudda import SnuddaLoad
import matplotlib.pyplot as plt


class PlotConnectivity:

    def __init__(self, network_path=None, network_name="default"):

        self.network_data = dict()
        self.connection_matrix = dict()
        self.distance_matrix = dict()

        if network_path is not None:
            self.add_network(network_path=network_path, network_name=network_name)

    def add_network(self, network_path, network_name):

        """ Load network.

        Args:
            network_path (str) : Path to network
            network_name (str) : Name of network (default: default)"""

        self.network_data[network_name] = SnuddaLoad(network_file=network_path)
        self.connection_matrix[network_name] = self.network_data[network_name].create_connection_matrix()
        self.distance_matrix[network_name] = self.network_data[network_name].create_distance_matrix()

    def add_plot(self,
                 pre_neuron_type, post_neuron_type,
                 network_name="default",
                 plot_label=None,
                 colour="black",
                 error_colour="grey",
                 bin_width=20e-6, num_bins=10,
                 pre_neuron_id=None, post_neuron_id=None,
                 axis=None):

        if pre_neuron_id is None:
            pre_neuron_id = self.network_data[network_name].get_neuron_id_of_type(neuron_type=pre_neuron_type)

        if post_neuron_id is None:
            post_neuron_id = self.network_data[network_name].get_neuron_id_of_type(neuron_type=post_neuron_type)

        sub_con_mat = self.connection_matrix[network_name][pre_neuron_id, :][:, post_neuron_id]
        sub_dist_mat = self.network_data[network_name].create_distance_matrix(pre_id=pre_neuron_id,
                                                                              post_id=post_neuron_id)

        n_connected, n_total, bin_edges = self.matrix_to_bins(con_mat=sub_con_mat, dist_mat=sub_dist_mat,
                                                              bin_width=bin_width, num_bins=num_bins)

        if plot_label is None:
            if network_name == "default":
                plot_label = f"{pre_neuron_type} to {post_neuron_type}"
            else:
                plot_label = f"{network_name}: {pre_neuron_type} to {post_neuron_type}"

        if axis is None:
            fig, axis = plt.subplots()

        z = 1.96  # 95 % confidence interval

        p = np.divide(n_connected + (z ** 2) / 2, n_total + z ** 2)
        p_error = np.multiply(np.divide(z, n_total + z ** 2),
                              np.sqrt(np.divide(np.multiply(n_connected, (n_total - n_connected)), n_total)
                                      + (z ** 2)/4))

        dist = (bin_edges[0:-1] + bin_edges[1:]) / 2 * 1e6

        axis.fill_between(dist, 100*(p-p_error), 100*(p+p_error), color=error_colour, step=None, alpha=0.4)
        axis.plot(dist, 100*p, color=colour, linewidth=2, label=plot_label)

        axis.legend(fontsize=18)

        axis.set_xlabel("Distance between somas ($\mu$m)", fontsize=22)
        axis.set_ylabel("Connection probability (%)", fontsize=22)

        return axis, (n_connected, n_total, dist)

    def draw_plot(self, axis, figure_path=None):

        plt.sca(axis)

        if figure_path is not None:
            plt.savefig(fname=figure_path)

        plt.show()

    @staticmethod
    def matrix_to_bins(con_mat, dist_mat, bin_width, num_bins):

        bin_edges = bin_width * np.arange(0, num_bins+1)
        n_connected = np.zeros((num_bins,), dtype=int)
        n_total = np.zeros((num_bins,), dtype=int)

        for is_connected, distance in zip(con_mat.todense().flatten().T, dist_mat.flatten().T):
            bin_idx = int(distance / bin_width)

            if bin_idx < num_bins:
                n_total[bin_idx] += 1

                if is_connected:
                    n_connected[bin_idx] += 1

        return n_connected, n_total, bin_edges


