# This script uses the functions defined in Network_analyse.py
#
# Performs analysis of the Striatum
#
#

import numpy as np
import sys
import os

import matplotlib.pyplot as plt

from snudda.analyse import SnuddaAnalyse


class SnuddaAnalyseStriatum(SnuddaAnalyse):

    def __init__(self, sim_dir, volume_type="cube", side_len=300e-6):

        if os.path.isfile(sim_dir):
            # We allow the user to also send in a hdf5 file as simDir...
            hdf5_file = sim_dir
            self.simDir = os.path.basename(sim_dir)
        else:
            self.simDir = sim_dir
            hdf5_file = sim_dir + "/network-pruned-synapses.hdf5"

            if not os.path.exists(hdf5_file):
                alt_hdf5_file = sim_dir + "/network-connect-voxel-pruned-synapse-file.hdf5"

                if os.path.exists(alt_hdf5_file):
                    hfd5_file = alt_hdf5_file

        print("Loading " + str(hdf5_file))

        super().__init__(hdf5_file=hdf5_file, load_cache=True,
                         volume_type=volume_type,
                         side_len=side_len)

    ############################################################################

    # Validation. How well do the synapse location for dSPN and iSPN match
    # the experimental data from Straub,..., Sabatini 2016

    def plot_fs_lts_cum_dist(self, plot_fs=True, plot_lts=True):

        pair_list_list = [[("FSN", "dSPN"), ("LTS", "dSPN")],
                          [("FSN", "iSPN"), ("LTS", "iSPN")]]
        figure_name_list = ["synapseCumulativeDistance-FSN-and-LTS-to-dSPN.png",
                            "synapseCumulativeDistance-FSN-and-LTS-to-iSPN.png"]
        figure_colour_list = [(6. / 255, 31. / 255, 85. / 255),
                              (150. / 255, 63. / 255, 212. / 255)]
        fill_range = [[0, 100e-6], [50e-6, 250e-6]]

        plot_flag = (plot_fs, plot_lts)

        assert plot_fs or plot_lts, "You must plot either FS or LTS, or both"

        if not plot_fs:
            figure_name_list = [x.replace("FSN-and-", "") for x in figure_name_list]
        if not plot_lts:
            figure_name_list = [x.replace("and-LTS-", "") for x in figure_name_list]

        for pairList, figName \
                in zip(pair_list_list, figure_name_list):

            plt.rcParams.update({'font.size': 22})
            fig = plt.figure()
            ax = plt.subplot(111)
            # fig.tight_layout()
            fig.subplots_adjust(bottom=0.15, left=0.15)

            for pair, figCol, fillR, plotMeFlag \
                    in zip(pairList, figure_colour_list, fill_range, plot_flag):

                if not plotMeFlag:
                    continue

                try:
                    pair_id = tuple([self.allTypes.index(x) for x in pair])
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    continue

                cum_dist = np.cumsum(self.dend_position_bin[pair_id]) \
                           / np.sum(self.dend_position_bin[pair_id])

                # Dont plot the full range
                end_idx = np.where(self.dend_position_edges <= 400e-6)[0][-1]

                ax.plot(self.dend_position_edges[:end_idx] * 1e6, cum_dist[:end_idx],
                        color=figCol, label=pair[0], linewidth=3)

                fill_idx = np.where(np.logical_and(fillR[0] <= self.dend_position_edges,
                                                   self.dend_position_edges <= fillR[1]))[0]
                fill_start = fill_idx[0]
                fill_end = fill_idx[-1]

                # Add the area marking
                ax.fill_between(self.dend_position_edges[fill_idx] * 1e6,
                                np.zeros((len(fill_idx),)),
                                cum_dist[fill_idx], alpha=0.95, color=figCol,
                                label=None)

                ax.set_xlabel('Distance from soma ($\mu$m)')
                ax.set_ylabel('Cumulative distrib.')

                if plot_fs and plot_lts:
                    ax.set_title("Synapse locations onto " + pair[1])
                else:
                    ax.set_title("Synapses " + pair[0] + " to " + pair[1])

            if plot_fs and plot_lts:
                # Only do legend if both are in figure
                ax.legend(loc="lower right")

            plt.ion()
            plt.show()
            plt.draw()
            plt.pause(0.0001)

            self.save_figure(plt, figName)

    ############################################################################


if __name__ == "__main__":

    sim_dir = None

    if len(sys.argv) > 1:
        sim_dir = sys.argv[1]
        print("Reading network from " + str(sim_dir))
    else:
        print("Please specify which directory the striatum network files is in")
        exit(-1)

    nas = SnuddaAnalyseStriatum(sim_dir, volume_type="cube")

    dist3D = False
    yMaxH = None  # 0.5

    nas.plot_incoming_connections(neuron_type="dSPN", pre_type="iSPN", num_bins=20)

    nas.plot_connection_probability("iSPN", "dSPN",
                                    dist_3d=dist3D,
                                    exp_max_dist=[50e-6, 100e-6],
                                    exp_data=[13 / 47.0, 10 / 80.0],
                                    exp_data_detailed=[(13, 47), (10, 80)],
                                    y_max=yMaxH)

    nas.plot_num_synapses_per_pair("iSPN", "dSPN")

    # nas.plotSynapseCumDist()
