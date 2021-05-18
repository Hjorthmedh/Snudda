import os
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

from snudda.analyse import SnuddaAnalyse
from snudda import SnuddaInit


class MethodsPaperFigure2:

    def __init__(self, network_file_list):

        self.analysis = []
        self.network_file_list = network_file_list

        for nf in self.network_file_list:
            print(f"Loading {nf}")
            self.analysis.append(SnuddaAnalyse(hdf5_file=nf,
                                               volume_type="cube"))

        self.side_len = 250e-6

        self.legend = ['No pruning',
                       'DP',
                       'DP, f1',
                       'DP, f1, SM',
                       'DP, f1, SM, mu2',
                       'DP, f1, SM, mu2, a3']

        self.plot_style = ['-', '-', '--', '--', ':', '-']
        self.plot_colour = [(0.5, 0.5, 0.5),
                            (0.75, 0.75, 0.75),
                            (0, 0, 0),
                            (0.5, 0.5, 0.5),
                            (0, 0, 0),
                            (0, 0, 0)]

    ############################################################################

    def make_connection_probability_summary(self,
                                            pre_type,
                                            post_type):

        fig, ax = plt.subplots(1)
        matplotlib.rcParams.update({'font.size': 24})

        x_max = 0.0
        y_max = 0.0
        n_bins = 86
        dist_3d = False

        connection_type = "synapses"

        for idx, a in enumerate(self.analysis):

            pre_id = a.populations[pre_type]
            post_id = a.populations[post_type]

            # We can in principle use all pairs, but here we restrict to just the
            # pairs who's post partner are in the centre
            post_id = a.get_sub_pop(volume_type=a.volume_type,
                                    volume_part="centre",
                                    side_len=self.side_len,
                                    neuron_id=post_id)

            (dist, p_con, count_con, count_all) = \
                a.connection_probability(pre_id, post_id, n_bins, dist_3d=dist_3d,
                                         connection_type=connection_type)

            print(f"!!!!    !!!!  count_con = {count_con}")

            d_half_step = (dist[1] - dist[0]) / 2
            plt.plot((dist + d_half_step) * 1e6, p_con,
                     color=self.plot_colour[idx],
                     linestyle=self.plot_style[idx],
                     linewidth=2)

            try:
                x_max = np.maximum((dist[-1] + d_half_step) * 1e6, x_max)
                y_max = np.maximum(np.nanmax(p_con), y_max)
            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)

                import pdb
                pdb.set_trace()

        plt.xticks(fontsize=14, rotation=0)
        plt.yticks(fontsize=14, rotation=0)

        label_size = 22

        plt.xlabel("Distance ($\mu$m)", fontsize=label_size)
        plt.ylabel("Con Prob (%)", fontsize=label_size)

        # Override max
        x_max = 400

        plt.xlim([0, x_max])
        plt.ylim([0, y_max])

        locs, labels = plt.yticks()
        new_labels = ["{0:g}".format(yy * 100) for yy in locs]

        if locs[0] < 0:
            locs = locs[1:]
            new_labels = new_labels[1:]

        if locs[-1] > plt.ylim()[1]:
            locs = locs[:-1]
            new_labels = new_labels[:-1]

        plt.yticks(locs, new_labels)

        font_p = matplotlib.font_manager.FontProperties()
        font_p.set_size('xx-small')
        plt.legend(self.legend, prop=font_p)

        plt.title(f"{self.analysis[0].neuron_name(pre_type)} to {self.analysis[0].neuron_name(post_type)}")
        plt.tight_layout()

        print("Debugging figure")
        import pdb
        pdb.set_trace()

        fig_name = 'Summary-pruning-dist-dep-connection-probability-' \
                   + str(pre_type) + "-to-" + str(post_type) \
                   + "-" + str(connection_type)

        self.analysis[0].save_figure(plt, fig_name, fig_type="png")

        plt.show()
        plt.pause(0.001)

    ############################################################################

    def make_num_synapses_summary_figure(self, pre_type, post_type):

        fig, ax = plt.subplots(1)
        matplotlib.rcParams.update({'font.size': 24})

        connection_type = "synapses"

        for idx, a in enumerate(self.analysis):

            pre_id = a.populations[pre_type]
            post_id = a.populations[post_type]

            if connection_type == "synapses":
                con_mat = a.connection_matrix
            elif connection_type == "gapjunctions":
                con_mat = a.connection_matrix_gj
            else:
                con_mat = None
                print(f"Unknown connection_type: {connection_type}")
                print("Please use 'synapses' or 'gapjunctions'")
                import pdb
                pdb.set_trace()

            # We can in principle use all pairs, but here we restrict to just the
            # pairs who's post partner are in the centre
            post_id = a.get_sub_pop(volume_type=a.volume_type,
                                    volume_part="centre",
                                    side_len=self.side_len,
                                    neuron_id=post_id)

            print("Calculating max synapses")
            max_synapses = con_mat[pre_id, :][:, post_id].max()

            # The prune tuning func might set data to 0, we want to exclude those
            mean_synapses = float(np.sum(con_mat[pre_id, :][:, post_id].data)) \
                            / np.sum(con_mat[pre_id, :][:, post_id].data != 0)

            con = con_mat[pre_id, :][:, post_id]

            # con = con.toarray()
            # existingCon = con[con != 0]

            existing_con = con[np.nonzero(con)].transpose()

            # Any connections? Otherwise skip plot
            if ((type(existing_con) == np.matrixlib.defmatrix.matrix
                 and len(existing_con) == 0)
                    or (type(existing_con) != np.matrixlib.defmatrix.matrix
                        and existing_con.getnnz() == 0)):
                continue

            print(f"Plotting {existing_con.shape[0]} connections")

            matplotlib.rcParams.update({'font.size': 22})

            plt.hist(existing_con,
                     # range(0,1+maxSynapses),
                     range(0, 21),
                     density=False,
                     align="left",
                     color=self.plot_colour[idx],
                     linestyle=self.plot_style[idx],
                     linewidth=2,
                     histtype=u"step")

            plt.xlabel(f"Number of {connection_type}")
            plt.ylabel('Count')

            plt.title(f"{self.analysis[0].neuron_name(pre_type)} to {self.analysis[0].neuron_name(post_type)}")

        plt.tight_layout()

        fig_name = f"Summary-network-number-of-{connection_type}-from-{pre_type}-to-{post_type}-per-cell"
        self.analysis[0].save_figure(plt, fig_name, fig_type="png")

        plt.show()
        plt.pause(0.001)


    ############################################################################

    def summary_plot_cum_dist(self, pre_type, post_type):

        fig, ax = plt.subplots(1)
        matplotlib.rcParams.update({'font.size': 24})

        pair = ("iSPN", "dSPN")
        connection_type = "synapses"

        for idx, a in enumerate(self.analysis):
            pair_id = tuple([a.allTypes.index(x) for x in pair])

            cum_dist = np.cumsum(a.dend_position_bin[pair_id]) / np.sum(a.dend_position_bin[pair_id])

            # Dont plot the full range
            end_idx = np.where(a.dend_position_edges <= 400e-6)[0][-1]

            ax.plot(a.dend_position_edges[:end_idx] * 1e6, cum_dist[:end_idx],
                    color=self.plot_colour[idx], label=self.legend[idx],
                    linestyle=self.plot_style[idx], linewidth=3)

        ax.set_xlabel('Distance from soma ($\mu$m)')
        ax.set_ylabel('Cumulative distrib.')

        fig_name = f"Summary-cumDist-of-{connection_type}-from-{pre_type}-to-{post_type}-per-cell"
        self.analysis[0].save_figure(plt, fig_name, fig_type="png")
        plt.show()

    ############################################################################

    @staticmethod
    def setup_network(network_path, config_name, network_type, n_neurons=2000, random_seed=None):

        # These parameters should be cleared
        param_lookup = {'No pruning': ["f1", "softMax", "mu2", "a3", "distPruning"],
                        'DP': ["f1", "softMax", "mu2", "a3"],
                        'DP, f1': ["softMax", "mu2", "a3"],
                        'DP, f1, SM': ["mu2", "a3"],
                        'DP, f1, SM, mu2': ["a3"],
                        'DP, f1, SM, mu2, a3': []}

        # Increase these for figure later
        n_dspn = int(n_neurons / 2)
        n_ispn = int(n_neurons / 2)

        neurons_dir = "$DATA/neurons/"
        config_file = os.path.join(network_path, config_name)

        si = SnuddaInit(config_file=config_file, struct_def={}, random_seed=random_seed)
        si.define_striatum(num_dSPN=n_dspn, num_iSPN=n_ispn, num_FS=0, num_LTS=0, num_ChIN=0,
                           volume_type="cube", neurons_dir=neurons_dir)
        si.network_data = MethodsPaperFigure2.remove_pruning(network_data=si.network_data,
                                                             pre_neuron="iSPN", post_neuron="dSPN",
                                                             synapse_type="GABA",
                                                             parameter_list=param_lookup[network_type])
        si.write_json()

    @staticmethod
    def remove_pruning(network_data, pre_neuron, post_neuron, synapse_type, parameter_list):

        for p in parameter_list:
            network_data["Connectivity"][f"{pre_neuron},{post_neuron}"][synapse_type]["pruning"][p] = None

        return network_data


if __name__ == '__main__':
    files = ['Neuroinformatics2020/Net10062-var-1/network-pruned-synapses.hdf5',
             'Neuroinformatics2020/Net10062-var-2/network-pruned-synapses.hdf5',
             'Neuroinformatics2020/Net10062-var-3/network-pruned-synapses.hdf5',
             'Neuroinformatics2020/Net10062-var-4/network-pruned-synapses.hdf5',
             'Neuroinformatics2020/Net10062-var-5/network-pruned-synapses.hdf5']

    legends = ['No pruning',
               'DP',
               'DP, f1',
               'DP, f1, SM',
               'DP, f1, SM, mu2',
               'DP, f1, SM, mu2, a3']

    mpf = MethodsPaperFigure2(files,
                              legends)

    mpf.make_connection_probability_summary('iSPN', 'dSPN')

    mpf.make_num_synapses_summary_figure('iSPN', 'dSPN')

    mpf.summary_plot_cum_dist('iSPN', 'dSPN')

    import pdb

    pdb.set_trace()
