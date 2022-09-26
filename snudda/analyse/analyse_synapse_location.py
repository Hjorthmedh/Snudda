import os
import numpy as np
from snudda.utils import SnuddaLoad
import matplotlib.pyplot as plt


class AnalyseSynapseLocation:

    def __init__(self, network_file):

        self.snudda_loader = SnuddaLoad(network_file=network_file, load_synapses=True)
        self.network_path = os.path.dirname(network_file)
        self.figure_path = os.path.join(self.network_path, "figures")

    def setup_figure_directory(self):

        if not os.path.isdir(self.figure_path):
            os.mkdir(self.figure_path)

    def count_synapses(self, pre_type, post_type):

        count = 0
        for syn_row in self.iterate_synapses(pre_type=pre_type, post_type=post_type):
            count += 1

    def synapse_distance_to_soma(self, pre_type=None, post_type=None):

        dist = []

        for syn_row in self.iterate_synapses(pre_type=pre_type, post_type=post_type):
            dist.append(syn_row[8]*1e-6)

        return dist

    def plot_synapse_distance_to_soma(self, pre_type=None, post_type=None, fig_path=None, bin_width=20):

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
            fig_path = os.path.join(self.figure_path, f"Distance-to-soma-{pre_text}-to-{post_text}.png")

        dist = self.synapse_distance_to_soma(pre_type=pre_type, post_type=post_type)

        plt.figure()
        dist_scaled = [d*1e6 for d in dist]
        bins = np.arange(0, np.ceil(np.max(dist_scaled) / bin_width + 1) * bin_width, bin_width)
        plt.hist(dist_scaled, bins=bins)
        plt.title(f"Distance to soma for {len(dist)} synapses from {pre_text} to {post_text}")
        plt.xlabel("Distance ($\mu$m)")
        plt.ylabel("Count")
        plt.ion()
        plt.show()

        if fig_path:
            plt.savefig(fig_path)
            print(f"Writing figure to {fig_path}")

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
