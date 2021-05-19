import numpy as np
from matplotlib import pyplot as plt
import sys
import os
import json
from collections import OrderedDict

from snudda.plotting.plot_node_benchmark import PlotNodeBenchmark


class PlotSizeBenchmark(PlotNodeBenchmark):

    def __init__(self, network_path_list=None, network_size=None, size_prefix=True):

        assert len(network_path_list) == len(network_size)

        self.data = self.merge_data(map(self.load_benchmark, network_path_list), sort=False)
        self.network_size = np.array(network_size)
        self.fig_dir = os.path.join(network_path_list[0], "figures")
        self.size_prefix = size_prefix

        assert (self.network_size == np.sort(self.network_size)).all(), "Networks must be sorted by size"

    def plot_data(self):

        stage = ["Place", "Detect", "Prune"]

        nodes = self.data["place"][:, 1]
        assert (nodes == nodes[0]).all(), "All sizes should have same number of nodes"

        duration = np.zeros((len(stage), len(nodes)))

        for idx, s in enumerate(stage):
            duration[idx, :] = self.data[s.lower()][:, 0] / 3600
            assert (self.data[s.lower()][:, 1] == nodes).all(), f"Stage {s} is missing some node values that Place has"

        small_size = 15
        medium_size = 22
        bigger_size = 30

        plt.rc('font', size=small_size)  # controls default text sizes
        plt.rc('axes', titlesize=small_size)  # fontsize of the axes title
        plt.rc('axes', labelsize=medium_size)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=small_size)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=small_size)  # fontsize of the tick labels
        plt.rc('legend', fontsize=small_size)  # legend fontsize
        plt.rc('figure', titlesize=bigger_size)  # fontsize of the figure title

        fig = plt.figure()
        ax = plt.subplot()

        ax.stackplot(self.network_size, duration, labels=stage)
        ax.legend(loc="upper left")
        ax.set_xlabel("Network size (log scale)")
        ax.set_ylabel("Runtime (h)")
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.set_xticks(self.network_size, minor=False)

        if self.size_prefix:
            ax.set_xticklabels(self.convert_labels(self.network_size))
        else:
            ax.set_xticklabels(self.network_size.astype(int))

        # Prevent cropping of x-label
        fig.subplots_adjust(bottom=0.15)

        if not os.path.isdir(self.fig_dir):
            os.mkdir(self.fig_dir)

        fig_name = os.path.join(self.fig_dir, "benchmark-size.pdf")
        plt.savefig(fig_name, dpi=300)
        plt.show()

    @staticmethod
    def convert_labels(values):

        prefix_lookup = [(1000000, "M"), (1000, "k")]
        converted = []

        for v in values:
            cv = f"{v}"
            for pval, prefix in prefix_lookup:
                if v % pval == 0:
                    cv = f"{int(v/pval)}{prefix}"
                    break
            converted.append(cv)

        return converted


if __name__ == "__main__":

    if len(sys.argv) > 2:
        network_path_list = sys.argv[1::2]
        network_size_list = [int(x) for x in sys.argv[2::2]]

        pb = PlotSizeBenchmark(network_path_list, network_size_list)
        pb.plot_data()

    else:
        print(f"Usage python3 {sys.argv[0]} network_path1 size1 network_path2 size2 ...")
