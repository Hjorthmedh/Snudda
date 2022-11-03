#!/usr/bin/env python3

import os
import numpy as np
import h5py
from snudda.utils import SnuddaLoad


class DegreeDistribution:

    def __init__(self, network_file):

        self.network_loader = SnuddaLoad(network_file=network_file, load_synapses=True)
        self.network_file = self.network_loader.network_file

        self.connection_matrix = self.network_loader.create_connection_matrix(sparse_matrix=False)
        self.degree_distribution = dict()

    def get_degree_correlation(self, pre_neuron, post_neuron):

        pre_neuron_id = self.network_loader.get_neuron_id_of_type(neuron_type=pre_neuron)
        post_neuron_id = self.network_loader.get_neuron_id_of_type(neuron_type=post_neuron)

        degree_correlation_histogram = np.zeros((len(post_neuron_id) + 1), dtype=np.uint)
        degree_correlation = np.zeros((len(pre_neuron_id), len(pre_neuron_id)), dtype=np.uint)

        map_idx = np.full((max(pre_neuron_id)+1), -1, dtype=int)
        for idx, id in enumerate(pre_neuron_id):
            map_idx[id] = idx

        sub_matrix = self.connection_matrix[pre_neuron_id, :][:, post_neuron_id]
        for pre_id_A in pre_neuron_id:
            for pre_id_B in pre_neuron_id:
                if pre_id_A >= pre_id_B:
                    continue

                # If connection matrix was sparse
                # con_vect_A = (self.connection_matrix[pre_id_A, :][:, post_neuron_id] > 0).toarray()
                # con_vect_B = (self.connection_matrix[pre_id_B, :][:, post_neuron_id] > 0).toarray()

                con_vect_A = (self.connection_matrix[pre_id_A, :][post_neuron_id] > 0)
                con_vect_B = (self.connection_matrix[pre_id_B, :][post_neuron_id] > 0)

                try:
                    degree = np.sum(np.logical_and(con_vect_A, con_vect_B))
                    degree_correlation_histogram[degree] += 1
                    degree_correlation[map_idx[pre_id_A], map_idx[pre_id_B]] = degree
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

        degree_correlation_histogram = np.trim_zeros(degree_correlation_histogram, trim="b")

        return degree_correlation, degree_correlation_histogram

    def write_all_degree_correlations(self, plot=True):

        file_name = os.path.join(os.path.dirname(self.network_file), f"degree-correlation.hdf5")

        print(f"Writing degree correlation to {file_name}")
        output = h5py.File(file_name, "w")

        neuron_types = self.network_loader.get_neuron_types(return_set=True)

        for pre_neuron in neuron_types:
            for post_neuron in neuron_types:
                print(f"Processing {pre_neuron} -> {post_neuron}")
                deg_corr, deg_corr_hist = self.get_degree_correlation(pre_neuron=pre_neuron, post_neuron=post_neuron)

                if np.count_nonzero(deg_corr) > 0:

                    try:
                        dc = output.create_dataset(f"{pre_neuron}_{post_neuron}", data=deg_corr)
                        dch = output.create_dataset(f"{pre_neuron}_{post_neuron}-hist", data=deg_corr_hist)
                        dc.attrs["pre"] = pre_neuron
                        dc.attrs["post"] = post_neuron
                        dc.attrs["type"] = "degree_correlation"
                        dch.attrs["pre"] = pre_neuron
                        dch.attrs["post"] = post_neuron
                        dch.attrs["type"] = "degree_correlation_histogram"
                    except:
                        import traceback
                        print(traceback.format_exc())
                        import pdb
                        pdb.set_trace()

        output.close()

    def write_in_out_degree(self):

        file_name = os.path.join(os.path.dirname(self.network_file), f"in-out-degree.hdf5")
        print(f"Write in/out degrees to {file_name}")
        output = h5py.File(file_name, "w")

        neuron_types = self.network_loader.get_neuron_types(return_set=True)

        for pre_neuron in neuron_types:
            for post_neuron in neuron_types:
                print(f"Counting connections {pre_neuron} -> {post_neuron}")

                in_deg = self.get_n_in(pre_type=pre_neuron, post_type=post_neuron)
                out_deg = self.get_n_out(pre_type=pre_neuron, post_type=post_neuron)

                output.create_dataset(f"in_degree_{pre_neuron}_to_{post_neuron}", data=in_deg)
                output.create_dataset(f"out_degree_{pre_neuron}_to_{post_neuron}", data=out_deg)

        output.close()

    def get_n_out(self, pre_type, post_type):

        pre_neuron_id = self.network_loader.get_neuron_id_of_type(neuron_type=pre_type)
        post_neuron_id = self.network_loader.get_neuron_id_of_type(neuron_type=post_type)

        n = np.sum(self.connection_matrix[pre_neuron_id, :][:, post_neuron_id] > 0, axis=1)
        return n

    def get_n_in(self, pre_type, post_type):

        pre_neuron_id = self.network_loader.get_neuron_id_of_type(neuron_type=pre_type)
        post_neuron_id = self.network_loader.get_neuron_id_of_type(neuron_type=post_type)

        n = np.sum(self.connection_matrix[pre_neuron_id, :][:, post_neuron_id] > 0, axis=0)
        return n

    def plot_all_degree_histogram(self):
        neuron_types = self.network_loader.get_neuron_types(return_set=True)

        for pre_neuron in neuron_types:
            self.plot_degree_histogram(pre_neuron)

    def plot_degree_histogram(self, pre_neuron):

        file_name = os.path.join(os.path.dirname(self.network_file), f"degree-correlation.hdf5")

        data = h5py.File(file_name, "r")

        import matplotlib.pyplot as plt
        fig = plt.figure()

        for d in (ds for name, ds in data.items()
                  if ds.attrs["type"] == "degree_correlation_histogram" and ds.attrs["pre"] == pre_neuron):

            degree = np.arange(0, d[()].size)
            plt.plot(degree, d[()], label=f"{pre_neuron}->{d.attrs['post']}")

        plt.legend(prop={'size': 6})
        plt.xlabel("Degree")
        plt.ylabel("Count")

        fig_name = os.path.join(os.path.dirname(self.network_file), f"degree_distribution_{pre_neuron}.png")
        plt.savefig(fig_name, dpi=300)
        print(f"Writing figure to {fig_name}")


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser("Extract degree distributions")
    parser.add_argument("network_file")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--process", action="store_true")

    args = parser.parse_args()

    dd = DegreeDistribution(args.network_file)

    if args.process:
        dd.write_all_degree_correlations()
        dd.write_in_out_degree()

    if args.plot:
        dd.plot_all_degree_histogram()

    if not args.plot and not args.process:
        print("Use --plot or --process to do something.")
