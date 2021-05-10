import numpy as np
from matplotlib import pyplot as plt
import sys
import os
import json
from collections import OrderedDict


class PlotNodeBenchmark:

    def __init__(self, network_path_list=None):

        self.data = self.merge_data(map(self.load_benchmark, network_path_list))
        self.fig_dir = os.path.join(network_path_list[0], "figures")

    def plot_data(self):

        # stage = ["Init", "Place", "Detect", "Prune"]
        stage = ["Place", "Detect", "Prune"]

        colour = np.array([[55, 126, 184], [77, 175, 74], [152, 78, 163], [255, 127, 0]]) / 255

        nodes = self.data["place"][:, 1]

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
        ax.stackplot(nodes, duration, labels=stage)
        ax.legend(loc="upper right")
        ax.set_xlabel("CPU cores (log scale)")
        ax.set_ylabel("Runtime (h)")
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.set_xticks(nodes, minor=False)
        ax.set_xticklabels(nodes.astype(int))

        # Prevent cropping of x-label
        fig.subplots_adjust(bottom=0.15)

        if not os.path.isdir(self.fig_dir):
            os.mkdir(self.fig_dir)

        fig_name = os.path.join(self.fig_dir, "benchmark-n-workers.pdf")
        plt.savefig(fig_name, dpi=300)
        plt.show()

    @staticmethod
    def load_benchmark(network_path):

        file_path = os.path.join(network_path, "benchmark_log.json")
        with open(file_path, "r") as f:
            data = json.load(f, object_pairs_hook=OrderedDict)

        new_data = OrderedDict()

        for main_keys in data.keys():

            for item, value in data[main_keys].items():
                new_data[item] = np.array(value)

        return new_data

    @staticmethod
    def merge_data(data_list, sort=True):

        new_data = OrderedDict()

        for data in data_list:
            for item, value in data.items():
                if item not in new_data:
                    new_data[item] = value
                else:
                    new_data[item] = np.concatenate([new_data[item], value])

        if sort:
            # Sort the items in node order
            for item, value in new_data.items():
                num_workers = np.unique(new_data[item][:, 1])
                d = []
                for nw in num_workers:
                    idx = np.where(new_data[item][:, 1] == nw)
                    mean_duration = np.mean(new_data[item][idx, 0])
                    d.append([mean_duration, nw])

                new_data[item] = np.array(d)

        return new_data


if __name__ == "__main__":

    if len(sys.argv) > 1:
        directories = sys.argv[1:]

        pb = PlotNodeBenchmark(network_path_list=directories)
        pb.plot_data()

    else:
        print("You must specify network folder(s) with benchmark data.")