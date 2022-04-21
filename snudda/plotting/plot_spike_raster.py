# Example:
#
# python3 plot_spike_raster.py save/traces/network-output-spikes-0.txt save/network-connect-synapse-file-0.hdf5


import sys
import os
import numpy as np
import re
import ntpath
from snudda.utils.load import SnuddaLoad
import time


class PlotSpikeRaster(object):

    def __init__(self, spike_file_name, network_file=None, skip_time=0.0, type_order=None, end_time=2.0,
                 figsize=None,type_division=None):
        #type_division: divides plot in two parts based on neuron type.
        #   Example 1: [["dspn","ispn"],["chin", "fsn","lts"]]
        #   Example 2: [["dspn","ispn","chin", "fsn","lts"],[]]
        self.spike_file_name = spike_file_name

        self.time = []
        self.spike_id = []
        self.end_time = end_time

        try:
            self.ID = int(re.findall('\d+', ntpath.basename(spike_file_name))[0])
        except:
            self.ID = 0

        self.neuron_name_remap = {"FSN": "FS"}

        self.read_csv()

        if network_file is not None:
            self.network_info = SnuddaLoad(network_file)
            self.network_file = network_file

            # assert(int(self.ID) == int(self.networkInfo.data["SlurmID"]))
        else:
            self.network_info = None
            self.network_file = None

        if self.network_info is None:
            print("If you also give network file, then the plot shows neuron types")
            self.plot_raster(skip_time=skip_time)
            time.sleep(1)
        else:
            if type_division is None:
                type_division = [["dspn", "ispn"], ["chin", "fsn", "lts"]]
            self.sort_traces()
            self.plot_colour_raster(skip_time=skip_time, type_order=type_order,type_division=type_division)
            time.sleep(1)

    ############################################################################

    def neuron_name(self, neuron_type):

        if neuron_type in self.neuron_name_remap:
            return self.neuron_name_remap[neuron_type]
        else:
            return neuron_type

    ############################################################################

    def read_csv(self):

        data = np.genfromtxt(self.spike_file_name, delimiter='\t')
        self.time = data[:, 0] * 1e-3
        self.spike_id = data[:, 1].astype(int)

    ############################################################################

    def plot_raster(self, skip_time=0, figsize=None):

        if not figsize:
            figsize = (6, 4)

        import matplotlib.pyplot as plt
        plt.figure(figsize=figsize)
        plt.scatter(self.time - skip_time, self.spike_id, color='black', s=1)
        plt.xlabel('Time (s)')
        plt.ylabel('Neurons')
        xl = plt.xlim()
        newx = (0, np.max(self.time) - skip_time)
        plt.xlim(newx)
        # plt.savefig('figures/Network-spike-raster-' + str(self.ID) + ".pdf")
        plt.savefig(f"figures/Network-spike-raster-{self.ID}.png", dpi=600)

        plt.ion()
        plt.show()
        plt.draw()
        plt.pause(0.001)

        print("Figure done")

    ############################################################################

    def plot_colour_raster(self, skip_time, plot_idx="sort", type_order=None, figsize=None, type_division=None):

        if plot_idx == "sort":
            plot_idx, tick_pos, tick_text, type_dict = self.sort_traces(type_order)
        else:
            tick_pos, tick_text = None, None

        plot_lookup = self.make_plot_lookup(plot_idx)

        import matplotlib.pyplot as plt

        cell_types = [n["name"].split("_")[0].lower() for n in self.network_info.data["neurons"]]

        colours = {"dSPN".lower(): (77. / 255, 151. / 255, 1.0),
                   "iSPN".lower(): (67. / 255, 55. / 255, 181. / 255),
                   "FS".lower(): (6. / 255, 31. / 255, 85. / 255),
                   "FSN".lower(): (6. / 255, 31. / 255, 85. / 255),
                   "ChIN".lower(): (252. / 255, 102. / 255, 0.0),
                   "LTS".lower(): (150. / 255, 63. / 255, 212. / 255)}

        cols = [colours[c] for c in cell_types]

        if not figsize:
            figsize = (6, 4)

        fig = plt.figure(figsize=figsize)
        r = 4
        grid = plt.GridSpec(r, r, hspace=0, wspace=0)
        ax = fig.add_subplot(grid[2:, :])
        if (len(type_division[0]) > 0) & (len(type_division[1]) > 0):
            atop = fig.add_subplot(grid[0, :])
            atop2 = fig.add_subplot(grid[1, :])
        elif (len(type_division[0]) == 0) & (len(type_division[1]) == 0):
            print("No neuron type to show in histogram due to empty variable type_division.")
        elif (len(type_division[0]) == 0) | (len(type_division[1]) == 0):
            atop = fig.add_subplot(grid[0:2, :])
            atop2 = atop

        t_idx = np.where(self.time >= skip_time)[0]

        cols2 = [colours[cell_types[int(s)]] for s in self.spike_id]

        ax.scatter(self.time[t_idx] - skip_time,
                   [plot_lookup[x] for x in self.spike_id[t_idx]],
                   color=[cols2[t] for t in t_idx], s=1,
                   linewidths=0.1)

        # histogram
        spike_count = []

        max0=0
        max1=0
    
        for t in type_order:
            pruned_spikes = self.time[type_dict[t]] - skip_time
            pruned_spikes = pruned_spikes[pruned_spikes>0]
            num_of_type = len([x["type"] for x in self.network_info.data["neurons"] if x["type"].lower() == t])
            bin_width = 0.01 #10  # ms
            bin_range = np.arange(0, self.end_time + skip_time + bin_width, bin_width)

            if t in type_division[0]:
                counts0, bins0, bars0 = atop.hist(pruned_spikes,
                                                  bins=bin_range,
                                                  range=(skip_time, self.time[-1]),
                                                  density=0,
                                                  color=colours[t],
                                                  alpha=1.0,
                                                  histtype='step',
                                                  weights=np.ones_like(pruned_spikes) / bin_width / num_of_type)
                spike_count.append((t, len(pruned_spikes), num_of_type))
                max0 = max(np.append(counts0, max0))
            elif t in type_division[1]:
                counts1, bins1, bars1 = atop2.hist(pruned_spikes,
                                                   bins=bin_range,
                                                   range=(skip_time, self.time[-1]),
                                                   density=0,
                                                   color=colours[t],
                                                   alpha=1.0,
                                                   histtype='step',
                                                   weights=np.ones_like(pruned_spikes) / bin_width / num_of_type)
                spike_count.append((t, len(pruned_spikes), num_of_type))
                max1 = max(np.append(counts1, max1))

        ax.invert_yaxis()

        ax.set_xlabel('Time (s)')
        if tick_pos is not None:
            ax.set_yticks(tick_pos)
            ax.set_yticklabels(tick_text)
        else:
            ax.ylabel('Neurons')

        # set axes ---------------------------------------------------
        atop.set_xticklabels([])
        atop.set_ylabel('Mean spikes/s')
        # UPDATE here to set specific range for plot window!!!!
        end_time = np.max([self.end_time, np.max(self.time)]) - skip_time
        for st, sc, n in spike_count:
            print(f"{st} ({n}): {sc/(n*end_time)} Hz, total spikes {sc}")
        if len(type_division[0]) > 0:
            atop.set_xlim([-0.01, end_time + 0.01])
            atop.set_yticks([0, round(max0*0.8)])  #Make sure ticks don't overlap between subplots
        if len(type_division[1]) > 0:
            atop2.set_xticklabels([])
            atop2.set_xlim([-0.01, end_time + 0.01])
            atop2.set_yticks([0, round(max1*0.8)])
        ax.set_xlim([-0.01, end_time + 0.01])
        m = len(self.network_info.data["neurons"])
        offset = m * 0.08  # 5%
        ax.set_ylim([-offset, m + offset])
        # -----------------------------------------------------------
        # plt.savefig('figures/Network-spike-raster-' + str(self.ID) + "-colour.pdf")

        # TODO: this is not working if run from the same folder as the networkFile
        # if so -> fig_path = "/figs"
        fig_path = os.path.join(os.path.dirname(self.network_file), "figures")
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)

        # have updated the name of the saved file to be the same as the fileName
        fn = os.path.basename(self.spike_file_name)
        #fig_name = '{}/{}{}'.format(fig_path, fn.split('.')[0], '-colour3.pdf')
        fig_name = '{}/{}{}'.format(fig_path, fn.split('.')[0], '-colour3.png')
        print(f"Saving {fig_name}")
        plt.savefig(fig_name, dpi=600)

        plt.ion()
        plt.show()
        plt.draw()
        plt.pause(0.001)

        ############################################################################

    def sort_traces(self, type_order=None):

        print("Sort the traces")

        all_types = [x["type"].lower() for x in self.network_info.data["neurons"]]

        if type_order is None:
            type_order = np.unique(all_types)

        # This one was fun to write. For every type t in typeOrder we want
        # to find all the indexes associated with it (here using enumerate and if)
        idx = [i for t in type_order for i, x in enumerate(all_types) if x == t]

        # new dict with cell type specific spikes
        type_dict = {'order': type_order}
        for t in type_order:
            type_dict[t] = [i for i, x in enumerate(self.spike_id) if all_types[x] == t]

        prev_pos = 0
        tick_pos = []
        tick_text = []

        for t in type_order:
            num_element = np.sum([x == t for x in all_types])

            if num_element == 0:  # No labels for missing types
                continue

            tick_pos.append(prev_pos + num_element / 2)
            tick_text.append(self.neuron_name(t))
            prev_pos += num_element

        return idx, tick_pos, tick_text, type_dict

    ############################################################################

    def make_plot_lookup(self, plot_idx):

        plot_lookup = dict()  # np.nan * np.zeros(len(plot_idx))

        for i, p in enumerate(plot_idx):
            plot_lookup[p] = i

        return plot_lookup

############################################################################


if __name__ == "__main__":

    # TODO: Update to use argparser
    print("Usage: " + sys.argv[0] + " network-output-spikes-XXX.txt")

    if len(sys.argv) > 1:
        file_name = sys.argv[1]
    else:
        file_name = None

    if len(sys.argv) > 2:
        network_file = sys.argv[2]
    else:
        network_file = None

    if len(sys.argv) > 3:
        end_time = float(sys.argv[3])
    else:
        end_time = 2.0

    if file_name is not None:
        # type_order = ["FS", "dSPN", "LTS", "iSPN", "ChIN"]
        type_order = ["fs", "fsn", "dspn", "lts", "ispn", "chin"]

        npsr = PlotSpikeRaster(file_name, network_file, skip_time=0.0,
                               end_time=end_time,
                               type_order=type_order)

    # import pdb
    # pdb.set_trace()
