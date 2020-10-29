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
                    hdf5_file = alt_hdf5_file

        print("Loading " + str(hdf5_file))

        super().__init__(hdf5_file=hdf5_file, load_cache=True,
                         volume_type=volume_type,
                         side_len=side_len)

    ############################################################################

    # Validation. How well do the synapse location for dSPN and iSPN match
    # the experimental data from Straub,..., Sabatini 2016

    def plot_fslt_scum_dist(self, plot_fs=True, plot_lts=True):

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

                # Don't plot the full range
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

    if len(sys.argv) > 1:
        sim_dir = sys.argv[1]
        print("Reading network from " + str(sim_dir))
    else:
        print("Please specify which directory the striatum network files is in")
        exit(-1)

    nas = SnuddaAnalyseStriatum(sim_dir, volume_type="cube")

    # import pdb
    # pdb.set_trace()
    #
    # nas.plotNeurons(0,showSynapses=True)

    plotHenrike = True
    plotChIN = True
    plotLTS = True

    dist3D = False
    # dist3D = True

    # import pdb
    # pdb.set_trace()

    nas.plot_synapse_cum_dist_summary(pair_list=[("dSPN", "ChIN"),
                                                 ("iSPN", "ChIN"),
                                                 ("LTS", "ChIN")])

    nas.plot_synapse_cum_dist_summary(pair_list=[("dSPN", "dSPN"),
                                                 ("iSPN", "dSPN"),
                                                 ("FSN", "dSPN"),
                                                 ("LTS", "dSPN"),
                                                 ("ChIN", "dSPN")])

    nas.plot_synapse_cum_dist_summary(pair_list=[("dSPN", "iSPN"),
                                                 ("iSPN", "iSPN"),
                                                 ("FSN", "iSPN"),
                                                 ("LTS", "iSPN"),
                                                 ("ChIN", "iSPN")])

    if plotHenrike:
        yMaxH = None  # 0.5

        nas.plot_connection_probability("dSPN", "iSPN", 
                                        dist_3d=dist3D, 
                                        exp_max_dist=[50e-6, 100e-6], 
                                        exp_data=[3 / 47.0, 3 / 66.0],
                                        exp_data_detailed=[(3, 47), (3, 66)],
                                        y_max=yMaxH)
        nas.plot_connection_probability("dSPN", "dSPN", 
                                        dist_3d=dist3D, 
                                        exp_max_dist=[50e-6, 100e-6], 
                                        exp_data=[5 / 19.0, 3 / 43.0],
                                        exp_data_detailed=[(5, 19), (3, 43)],
                                        y_max=yMaxH)
        nas.plot_connection_probability("iSPN", "dSPN", 
                                        dist_3d=dist3D, 
                                        exp_max_dist=[50e-6, 100e-6], 
                                        exp_data=[13 / 47.0, 10 / 80.0],
                                        exp_data_detailed=[(13, 47), (10, 80)],
                                        y_max=yMaxH)
        nas.plot_connection_probability("iSPN", "iSPN",
                                        dist_3d=dist3D,
                                        exp_max_dist=[50e-6, 100e-6],
                                        exp_data=[14 / 39.0, 7 / 31.0],
                                        exp_data_detailed=[(14, 39), (7, 31)],
                                        y_max=yMaxH)

    nas.plot_num_synapses_per_pair("dSPN", "dSPN")
    nas.plot_num_synapses_per_pair("dSPN", "iSPN")
    nas.plot_num_synapses_per_pair("iSPN", "dSPN")
    nas.plot_num_synapses_per_pair("iSPN", "iSPN")

    # !!! Check edge effects

    nas.plot_incoming_connections(neuron_type="iSPN", pre_type="FSN")
    nas.plot_incoming_connections(neuron_type="iSPN", pre_type="ChIN")
    nas.plot_incoming_connections(neuron_type="iSPN", pre_type="LTS")

    nas.plot_incoming_connections(neuron_type="dSPN", pre_type="dSPN")
    nas.plot_incoming_connections(neuron_type="dSPN", pre_type="iSPN")
    nas.plot_incoming_connections(neuron_type="iSPN", pre_type="dSPN")
    nas.plot_incoming_connections(neuron_type="iSPN", pre_type="iSPN")

    if True:
        # 2-5 ChIN should connect to each MS (approx)
        nas.plot_incoming_connections(neuron_type="dSPN", pre_type="ChIN")
        nas.plot_incoming_connections(neuron_type="iSPN", pre_type="ChIN")

    if True:
        nas.plot_connection_probability("FSN", "iSPN",
                                        dist_3d=dist3D,
                                        exp_max_dist=[100e-6, 150e-6, 250e-6],
                                        exp_data=[6 / 9.0, 21 / 54.0, 27 / 77.0],
                                        exp_data_detailed=[(6, 9), (21, 54), (27, 77)],
                                        y_max=None)

        nas.plot_connection_probability("FSN", "dSPN",
                                        dist_3d=dist3D,
                                        exp_max_dist=[100e-6, 150e-6, 250e-6],
                                        exp_data=[8 / 9.0, 29 / 48.0, 48 / 90.0],
                                        exp_data_detailed=[(8, 9), (29, 48), (48, 90)],
                                        y_max=None)

        nas.plot_num_synapses_per_pair("FSN", "dSPN")
        nas.plot_num_synapses_per_pair("FSN", "iSPN")

        #  Gittis,...,Kreitzer 2010 (p2228) -- 7/12 (and 3/4 reciprocal) -- distance?
        # FS->FS synapses weaker, 1.1 +/- 1.5nS

        nas.plot_connection_probability("FSN", "FSN",
                                        dist_3d=dist3D,
                                        exp_max_dist=[250e-6],
                                        exp_data=[7 / 12.0],
                                        exp_data_detailed=[(7, 12)])

        nas.plot_num_synapses_per_pair("FSN", "FSN")

        # Koos & Tepper 1999, 2/6
        nas.plot_connection_probability("FSN", "FSN",
                                        dist_3d=dist3D,
                                        connection_type="gapjunctions",
                                        exp_max_dist=[250e-6, 250e-6],
                                        exp_data=[2 / 6.0, 3 / 7.0],
                                        exp_data_detailed=[(2, 6), (3, 7)], )

        nas.plot_num_synapses_per_pair("FSN", "FSN", connection_type="gapjunctions")

        nas.plot_incoming_connections(neuron_type="FSN", pre_type="FSN",
                                      connection_type="gapjunctions")

    nas.plot_fslt_scum_dist()
    nas.plot_fslt_scum_dist(plot_fs=False)
    nas.plot_fslt_scum_dist(plot_lts=False)

    nas.plot_num_synapses_per_pair("dSPN", "ChIN")
    nas.plot_num_synapses_per_pair("iSPN", "ChIN")

    nas.plot_incoming_connections(neuron_type="FSN", pre_type="FSN")
    nas.plot_incoming_connections(neuron_type="FSN", pre_type="FSN")

    nas.plot_num_synapses_per_pair("ChIN", "FSN")

    nas.plot_synapse_cum_dist()

    nas.plot_synapse_dist(density_flag=True)

    if plotLTS:
        # 3/21 LTS->MS, Basal Ganglia book --- distance??
        # Ibanez-Sandoval, ..., Tepper  2011 3/21 -- if patching around visual axon
        # but 2/60 when patching blind
        nas.plot_connection_probability("LTS", "dSPN",
                                        dist_3d=dist3D,
                                        exp_max_dist=[250e-6],
                                        exp_data=[2 / 60.0],
                                        exp_data_detailed=[(2, 60)],
                                        x_max=500)

        nas.plot_connection_probability("LTS", "iSPN",
                                        dist_3d=dist3D,
                                        exp_max_dist=[250e-6],
                                        exp_data=[2 / 60.0],
                                        exp_data_detailed=[(2, 60)],
                                        x_max=500)

        # Silberberg et al 2013, 2/12 FS-> LTS connected --- distance??
        # Voltage deflection... 0.5mV and 0.8mV
        # (check Szydlowski et al 2013, what Cl rev)
        #
        nas.plot_connection_probability("FSN", "LTS",
                                        dist_3d=dist3D,
                                        exp_max_dist=[250e-6],
                                        exp_data=[2.0 / 12],
                                        exp_data_detailed=[(2, 12)])

        nas.plot_num_synapses_per_pair("LTS", "dSPN")
        nas.plot_num_synapses_per_pair("LTS", "iSPN")
        nas.plot_num_synapses_per_pair("LTS", "ChIN")
        nas.plot_num_synapses_per_pair("FSN", "LTS")

        # This plots figures for the article

    # 100e-6 from Planert 2010, and 250e-6 data from Gittis 2010
    # 150e-6 from Gittis 2011 (actually 100 +/- 50 micrometers)

    # MS <-> MS

    # FS -> MS

    nas.plot_num_synapses_per_pair("ChIN", "dSPN")
    nas.plot_num_synapses_per_pair("ChIN", "iSPN")
    nas.plot_num_synapses_per_pair("ChIN", "LTS")

    nas.plot_connection_probability("ChIN", "LTS",
                                    dist_3d=dist3D)

    # Janicova 2015?? --- distance??!
    nas.plot_connection_probability("ChIN", "iSPN",
                                    dist_3d=dist3D,
                                    exp_max_dist=[250e-6],
                                    exp_data=[0.05])

    nas.plot_connection_probability("ChIN", "dSPN",
                                    dist_3d=dist3D,
                                    exp_max_dist=[250e-6],
                                    exp_data=[0.05])

    if True:
        nas.plot_connection_probability("LTS", "ChIN",
                                        dist_3d=dist3D)

        # ALSO ADD GAP JUNCTIONS PLOT!!!
        # No exp data for this -- Gittis,...,Kreitzer 2010 (p2228) -- 7/12 (and 3/4 reciprocal) -- distance?
        # FS->FS synapses weaker, 1.1 +/- 1.5nS

    if plotChIN:
        # "I Janickova et al. 2017 så har de 2018 varicosities i en area på 655 um²,
        # deras slices är 70 um tjocka och om man antar att det inte är några
        # varicositites som täcker varandra så är volym-densiteten/mm³: 4.4*10⁷/mm3"
        # 1.7e6/24*0.01 = 708 ChIN per mm3
        # 4.4e7 / 708 = 62000 varicosities per ChIN
        #
        # 325 ChIN synapser per MS
        # 2-5 ChIN per MS
        # --> 65-160 synapser between a ChIN-MS pair
        # --> Each ChIN connect to 400 - 950 MS
        #
        # Number of MS within 350 micrometer radius 4*pi*(350e-6)^3/3*1.76e6/24e-9
        # --> 13100 MS reachable by ChIN at most (or rather number of MS somas
        # within radius of axonal arbour)
        # -->  3-7% connectivity probability??

        # nas.plotConnectionProbability("ChIN","iSPN", \
        #                              dist3D=dist3D,
        #                              expMaxDist=[200e-6],
        #                              expData=[62/89.0],
        #                              expDataDetailed=[(62,89)],
        #                              yMax=1.0)
        # This is from a targeted experiment, when they looked at where axon were?
        #
        # nas.plotConnectionProbability("ChIN","dSPN", \
        #                              dist3D=dist3D,
        #                              expMaxDist=[200e-6],
        #                              expData=[62/89.0],
        #                              expDataDetailed=[(62,89)],
        #                              yMax=1.0)

        nas.plot_connection_probability("ChIN", "FSN",
                                        dist_3d=dist3D,
                                        y_max=None)

        # A MS neuron receives 1e4 assymetrical synapses (Kincaid et al 1998),
        # and 2500 symmetrical synapses (Ingham et al 1998). Symmetrical synapses
        # can be dopaminergic, cholinergic or GABAergic, with dopaminergic
        # being 13% (Roberts et al 2002). Assuming that cholinergic inputs are
        # a similar percentage, 650 symmetrical synapses per MS are not GABAergic.
        #
        # --> 0.13*2500 = 325 ChIN inputs to MS
        nas.plot_incoming_connections(neuron_type="ChIN", pre_type="dSPN")
        nas.plot_incoming_connections(neuron_type="ChIN", pre_type="iSPN")
        nas.plot_incoming_connections(neuron_type="ChIN", pre_type="LTS")

        # 2-5 ChIN should connect to each MS (approx) --- ref? ?!?!?!
        nas.plot_incoming_connections(neuron_type="dSPN", pre_type="ChIN")
        nas.plot_incoming_connections(neuron_type="iSPN", pre_type="ChIN")

        # Om vi antar 2000 MS skulle kunna nå varje ChIN, när 10% aktiverade
        # så är 200 MS aktiva, om 75% av ChIN känner av MS input
        # (1-p)^200 = 0.25 --> 0.7 %
        nas.plot_connection_probability("dSPN", "ChIN",
                                        dist_3d=dist3D)
        nas.plot_connection_probability("iSPN", "ChIN",
                                        dist_3d=dist3D)

        # nas.nearestPreNeighbourDistance("LTS","dSPN")
        # nas.nearestPreNeighbourDistance("LTS","iSPN")

    if True:
        nas.plot_incoming_connections(neuron_type="dSPN", pre_type="iSPN")
        nas.plot_incoming_connections(neuron_type="dSPN", pre_type="dSPN")
        nas.plot_incoming_connections(neuron_type="dSPN", pre_type="FSN")

        nas.plot_incoming_connections(neuron_type="iSPN", pre_type="iSPN")
        nas.plot_incoming_connections(neuron_type="iSPN", pre_type="dSPN")
        nas.plot_incoming_connections(neuron_type="iSPN", pre_type="FSN")

        nas.plot_incoming_connections(neuron_type="dSPN", pre_type="LTS")
        nas.plot_incoming_connections(neuron_type="iSPN", pre_type="LTS")
        nas.plot_incoming_connections(neuron_type="ChIN", pre_type="LTS")

        nas.plot_incoming_connections(neuron_type="LTS", pre_type="ChIN")
        nas.plot_incoming_connections(neuron_type="LTS", pre_type="FSN")

        nas.plot_incoming_connections(neuron_type="ChIN", pre_type="dSPN")
        nas.plot_incoming_connections(neuron_type="ChIN", pre_type="iSPN")
