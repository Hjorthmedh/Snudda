import numpy as np
import pandas as pd
from snudda import SnuddaLoad
import matplotlib.pyplot as plt


class PlotConnectivity:

    def __init__(self, network_path=None, network_name="default"):

        self.network_data = dict()
        self.connection_matrix = dict()
        self.distance_matrix = dict()
        self.exp_and_model_data = []
        self.exp_data = None

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

        p, p_error = self.estimate_p(n_connected=n_connected, n_total=n_total)

        dist = (bin_edges[0:-1] + bin_edges[1:]) / 2 * 1e6

        axis.fill_between(dist, 100*(p-p_error), 100*(p+p_error), color=error_colour, step=None, alpha=0.4)
        axis.plot(dist, 100*p, color=colour, linewidth=2, label=plot_label)

        axis.legend(fontsize=18)

        axis.set_xlabel("Distance between somas ($\mu$m)", fontsize=22)
        axis.set_ylabel("Connection probability (%)", fontsize=22)

        return axis, (n_connected, n_total, dist)

    def estimate_p(self, n_connected, n_total, z=1.96):  # 95 % confidence interval

        p = np.divide(n_connected + (z ** 2) / 2, n_total + z ** 2)
        p_error = np.multiply(np.divide(z, n_total + z ** 2),
                              np.sqrt(np.divide(np.multiply(n_connected, (n_total - n_connected)), n_total)
                                      + (z ** 2)/4))

        return p, p_error

    def draw_plot(self, axis, figure_path=None):

        plt.sca(axis)

        if figure_path is not None:
            plt.savefig(fname=figure_path)

        plt.show()

    def add_experimental_data(self, network_name,
                              pre_neuron_type, post_neuron_type,
                              distance_range, n_connected, n_total,
                              exp_colour="red",
                              model_colour="black",
                              label=None):

        # Calculate the corresponding model data
        pre_neuron_id = self.network_data[network_name].get_neuron_id_of_type(neuron_type=pre_neuron_type)
        post_neuron_id = self.network_data[network_name].get_neuron_id_of_type(neuron_type=post_neuron_type)

        sub_con_mat = self.connection_matrix[network_name][pre_neuron_id, :][:, post_neuron_id]
        sub_dist_mat = self.network_data[network_name].create_distance_matrix(pre_id=pre_neuron_id,
                                                                              post_id=post_neuron_id)

        if len(distance_range) != 2 or distance_range[1] < distance_range[0]:
            raise ValueError("distance_range must be a tuple (start_dist, end_dist) in meters")

        bin_width = distance_range[1] - distance_range[0]
        num_bins = 1

        n_connected_model, n_total_model, bin_edges = self.matrix_to_bins(con_mat=sub_con_mat, dist_mat=sub_dist_mat,
                                                              bin_width=bin_width, num_bins=num_bins)

        if label is None:
            if network_name == "default":
                label = f"{pre_neuron_type} to {post_neuron_type}\n{distance_range[0]*1e6:.0f}-{distance_range[1]*1e6:.0f} $\mu$m"
            else:
                label = f"{pre_neuron_type} to {post_neuron_type} ({network_name})\n{distance_range[0] * 1e6:.0f}-{distance_range[1] * 1e6:.0f} $\mu$m"

        model_label = fr"{label} (model)"
        exp_label = fr"{label}"

        self.exp_and_model_data.append(((n_connected, n_total, exp_label, exp_colour),
                                        (n_connected_model[0], n_total_model[0], model_label, model_colour)))

    def clear_plot_data(self):
        self.exp_and_model_data = []

    def draw_exp_model_data(self, figure_path=None, label_rotation=90):

        fig, axis = plt.subplots()

        n_data_sets = len(self.exp_and_model_data)
        data_set_id = np.arange(n_data_sets)

        label_text = []
        label_pos = []

        for ctr, ((n_con_exp, n_total_exp, exp_label, exp_colour),
                  (n_con_model, n_total_model, model_label, model_colour)) in \
                enumerate(self.exp_and_model_data):

            p_exp, p_exp_error = self.estimate_p(n_connected=n_con_exp, n_total=n_total_exp)
            p_model, p_model_error = self.estimate_p(n_connected=n_con_model, n_total=n_total_model)

            exp_alpha = 1  # 1 - 2*p_exp_error / p_exp
            model_alpha = 1  # 1 - 2*p_model_error / p_model

            axis.bar(ctr-0.15, height=p_exp, yerr=p_exp_error, width=0.3, color=exp_colour, align="center", alpha=exp_alpha, capsize=5)
            axis.bar(ctr+0.15, height=p_model, yerr=p_model_error, width=0.3, color=model_colour, align="center", alpha=model_alpha, capsize=5)

            label_text.append(exp_label)
            label_pos.append(ctr)

        axis.set_xticks(label_pos, label_text, rotation=label_rotation)

        if figure_path is not None:
            plt.savefig(figure_path)

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

    def load_exp_data(self, data_file):

        self.exp_data = pd.read_csv(filepath_or_buffer=data_file)

    def add_loaded_exp_data_to_plot(self, phenotype, data_source, network_name,
                                    pre_neuron_type, post_neuron_type,
                                    exp_colour="red", model_colour="black"):

        idx = (self.exp_data["phenotype"] == phenotype) & (self.exp_data["data_source"] == data_source) \
            & (self.exp_data["pre_neuron"] == pre_neuron_type) & (self.exp_data["post_neuron"] == post_neuron_type)

        mean_dist = self.exp_data[idx]["mean_dist"].iloc[0]
        std_dist = self.exp_data[idx]["std_dist"].iloc[0]
        min_dist = self.exp_data[idx]["min_dist"].iloc[0]
        max_dist = self.exp_data[idx]["max_dist"].iloc[0]
        n_connected = self.exp_data[idx]["n_connected"].iloc[0]
        n_total = self.exp_data[idx]["n_total"].iloc[0]

        if not np.isnan(mean_dist) and not np.isnan(std_dist):
            distance_range = np.array([mean_dist - 1.96*std_dist, mean_dist + 1.96*std_dist])
        elif not np.isnan(min_dist) and not np.isnan(max_dist):
            distance_range = np.array([min_dist, max_dist])
        else:
            raise ValueError(f"No usable distances for {phenotype}, {data_source}, {pre_neuron_type}, {post_neuron_type}")

        self.add_experimental_data(network_name=network_name,
                                   pre_neuron_type=pre_neuron_type,
                                   post_neuron_type=post_neuron_type,
                                   distance_range=distance_range,
                                   n_connected=n_connected,
                                   n_total=n_total,
                                   exp_colour=exp_colour,
                                   model_colour=model_colour)
