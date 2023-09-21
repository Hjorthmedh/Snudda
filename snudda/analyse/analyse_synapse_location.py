import os
import numpy as np
from snudda.utils import SnuddaLoad
import matplotlib.pyplot as plt


class AnalyseSynapseLocation:

    def __init__(self, network_file):

        self.snudda_loader = SnuddaLoad(network_file=network_file, load_synapses=True)
        self.network_path = os.path.dirname(self.snudda_loader.network_file)
        self.figure_path = os.path.join(self.network_path, "figures")

        self.connection_matrix = None

    def setup_figure_directory(self):

        if not os.path.isdir(self.figure_path):
            os.mkdir(self.figure_path)

    def count_synapses(self, pre_type, post_type):

        count = 0
        for syn_row in self.iterate_synapses(pre_type=pre_type, post_type=post_type):
            count += 1

    def synapse_distance_to_soma(self, pre_type=None, post_type=None):

        dist = np.zeros((self.snudda_loader.data["synapses"].shape[0],), dtype=np.float32)
        ctr = 0

        for idx, syn_row in enumerate(self.iterate_synapses(pre_type=pre_type, post_type=post_type)):
            dist[idx] = syn_row[8]*1e-6
            ctr += 1

        return dist[:ctr]

    def synapses_per_pair(self, pre_type=None, post_type=None):

        if self.connection_matrix is None:
            self.connection_matrix = self.snudda_loader.create_connection_matrix()

        pre_id = self.snudda_loader.get_neuron_id_of_type(neuron_type=pre_type)
        post_id = self.snudda_loader.get_neuron_id_of_type(neuron_type=post_type)

        sub_mat = self.connection_matrix[pre_id, :][:, post_id]

        return sub_mat[sub_mat > 0].todense().flatten().T

    def plot_synapses_per_pair(self, pre_type=None, post_type=None, fig_path=None,
                               figure=None, fig_size=None,
                               label=None, title=None, show_plot=True, colour=None, linestyle=None):

        self.setup_figure_directory()

        if pre_type is None:
            pre_text = "ALL"
        else:
            pre_text = pre_type

        if post_type is None:
            post_text = "ALL"
        else:
            post_text = post_type

        if fig_path is None:
            fig_path = os.path.join(self.figure_path, f"Distance-to-soma-{pre_text}-to-{post_text}.pdf")

        num_synapses = self.synapses_per_pair(pre_type=pre_type, post_type=post_type)

        if len(num_synapses) == 0:
            print(f"No synapses found between {pre_type} and {post_type}")
            return None

        if figure is None:
            figure = plt.figure(figsize=fig_size)
        else:
            plt.figure(figure.number)

        bins = np.arange(0, self.connection_matrix.tocsr().max()+1)
        try:
            plt.hist(num_synapses, bins=bins, label=label, histtype="step", color=colour, linestyle=linestyle)
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

        if title is None:
            plt.title(f"Number of synapses between {pre_text} and {post_text} pair")
        else:
            plt.title(title)

        plt.xlabel("Number of synapses")
        plt.ylabel("Count")

        if show_plot:
            plt.ion()
            plt.show()

        if fig_path:
            plt.savefig(fig_path)
            print(f"Writing figure to {fig_path}")

        return figure

    def plot_synapse_distance_to_soma(self, pre_type=None, post_type=None, fig_path=None, bin_width=20,
                                      figure=None, fig_size=None,
                                      label=None, title=None, show_plot=True, colour=None, linestyle=None):

        self.setup_figure_directory()

        if pre_type is None:
            pre_text = "ALL"
        else:
            pre_text = pre_type

        if post_type is None:
            post_text = "ALL"
        else:
            post_text = post_type

        if fig_path is None:
            fig_path = os.path.join(self.figure_path, f"Distance-to-soma-{pre_text}-to-{post_text}.pdf")

        dist = self.synapse_distance_to_soma(pre_type=pre_type, post_type=post_type)

        if len(dist) == 0:
            print(f"No synapses found between {pre_type} and {post_type}")
            return None

        if figure is None:
            figure = plt.figure(figsize=fig_size)
        else:
            plt.figure(figure.number)

        dist_scaled = dist * 1e6
        bins = np.arange(0, np.ceil(np.max(dist_scaled) / bin_width + 1) * bin_width, bin_width)
        plt.hist(dist_scaled, bins=bins, label=label, histtype="step", color=colour, linestyle=linestyle)
        if title is None:
            plt.title(f"Distance to soma for {len(dist)} synapses from {pre_text} to {post_text}")
        else:
            plt.title(title)
        plt.xlabel("Distance ($\mu$m)")
        plt.ylabel("Count")

        if show_plot:
            plt.ion()
            plt.show()

        if fig_path:
            plt.savefig(fig_path)
            print(f"Writing figure to {fig_path}")

        return figure

    def iterate_synapses(self, pre_type=None, post_type=None):

        neuron_types = self.snudda_loader.get_neuron_types()

        for syn_row in self.snudda_loader.data["synapses"]:

            if pre_type is not None:
                pret = neuron_types[syn_row[0]]
                if pre_type != pret:
                    continue

            if post_type is not None:
                postt = neuron_types[syn_row[1]]
                if post_type != postt:
                    continue

            yield syn_row


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser("Analyse synapse distribution")
    parser.add_argument("network_file")
    parser.add_argument("--pre", type=str, default=None)
    parser.add_argument("--post", type=str, default=None)

    args = parser.parse_args()

    asl = AnalyseSynapseLocation(network_file=args.network_file)

    asl.plot_synapse_distance_to_soma(pre_type=args.pre, post_type=args.post)
