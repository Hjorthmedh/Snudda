#!/usr/bin/env python3

import os
import numpy as np

from snudda.utils.load import SnuddaLoad


class SnuddaExportConnectionMatrix(object):
    """ Exports a connection matrix from network. """

    def __init__(self, in_file, out_file, save_sparse=True, save_on_init=True):

        """ Constructor.

        Args:
            in_file : Network file
            out_file : Connection matrix file
            save_sparse : Should data be saved in sparse format?
        """

        self.sl = SnuddaLoad(in_file)

        self.out_file = out_file
        self.out_file_meta = f"{out_file}-meta"
        self.save_sparse = save_sparse

        data = self.sl.data

        self.con_mat = self.create_connection_matrix()
        self.neuron_type = [x["type"] for x in data["neurons"]]
        self.neuron_name = [x["name"] for x in data["neurons"]]
        self.neuron_morph = [x["morphology"] for x in data["neurons"]]
        self.population_unit = data["population_unit"]

        self.pos = data["neuron_positions"]

        if save_on_init:
            self.save()

    def save(self):

        print(f"Writing {self.out_file} (row = src, column=dest)")
        print(f"Saving {np.sum(np.sum(self.con_mat))} synapses, over {np.count_nonzero(self.con_mat)} coupled pairs.")

        if self.save_sparse:
            x_pos, y_pos = np.where(self.con_mat)
            sparse_data = np.zeros((len(x_pos), 3), dtype=int)
            for idx, (x, y) in enumerate(zip(x_pos, y_pos)):
                sparse_data[idx, :] = [x, y, self.con_mat[x, y]]

            # np.savetxt(self.out_file, sparse_data, delimiter=",", fmt="%d")
            np.save(self.out_file, sparse_data.astype(np.int32))

            # Test to verify
            for row in sparse_data:
                assert self.con_mat[row[0], row[1]] == row[2]
        else:
            np.save(self.out_file, self.con_mat.astype(np.int32))
            # np.savetxt(self.out_file, self.con_mat, delimiter=",", fmt="%d")

        print("Writing " + self.out_file_meta)
        with open(self.out_file_meta, "w") as f_out_meta:
            for i, (nt, nn, p, mf, pu) in enumerate(zip(self.neuron_type, self.neuron_name, self.pos, self.neuron_morph,
                                                        self.population_unit)):
                s = "%d,%s,%s,%f,%f,%f,%s, %d\n" % (i, nt, nn, p[0], p[1], p[2], mf, pu)
                f_out_meta.write(s)
            f_out_meta.close()

    ############################################################################

    def create_axon_dend_distance_matrix(self):
        # Axon speed 0.8m/s  # Tepper and Lee 2007, Wilson 1986, Wilson 1990

        num_neurons = self.sl.data["num_neurons"]

        axon_dend_dist_matrix = np.full((num_neurons, num_neurons), fill_value=np.nan, dtype=float)

        # This calculates the average delay (due to axon distance) between two neurons
        pre_id = None
        post_id = None
        synapse_distance = []

        for synapse_set in self.sl.synapse_iterator():
            for synapse in synapse_set:
                if pre_id != synapse[0] or post_id != synapse[1]:
                    if len(synapse_distance) > 0:
                        axon_dend_dist_matrix[pre_id, post_id] = np.mean(synapse_distance)
                    pre_id = synapse[0]
                    post_id = synapse[1]
                    synapse_distance = [(synapse[7] + synapse[8]) * 1e-6]
                else:
                    synapse_distance.append((synapse[7] + synapse[8]) * 1e-6)

        if len(synapse_distance) > 0:
            axon_dend_dist_matrix[pre_id, post_id] = np.mean(synapse_distance)

        return axon_dend_dist_matrix

    ############################################################################

    def plot_matrix(self, matrix, hide_nan=True):

        import matplotlib.pyplot as plt

        if hide_nan:
            plot_matrix = matrix.copy()
            plot_matrix[np.isnan(plot_matrix)] = 0
        else:
            plot_matrix = matrix

        plt.spy(plot_matrix)
        plt.show()

    ############################################################################

    def save_axon_dend_distance_matrix(self, plot=False):

        delay_file = f"{self.out_file}-path-distance"
        print(f"Writing path distance file: {delay_file}")

        axon_dend_distance = self.create_axon_dend_distance_matrix()

        if plot:
            self.plot_matrix(axon_dend_distance, hide_nan=True)

        save_matrix = axon_dend_distance.copy()
        save_matrix[np.isnan(save_matrix)] = np.inf

        # np.savetxt(delay_file, save_matrix, delimiter=",", fmt="%d")
        np.save(delay_file, save_matrix.astype(np.float32))

    ############################################################################

    def save_distance_matrix(self, plot=False):

        dist_file = f"{self.out_file}-dist"
        print(f"Writing distance file: {dist_file}")

        dist_matrix = self.sl.create_distance_matrix()

        if plot:
            self.plot_matrix(dist_matrix)

        np.save(dist_file, dist_matrix.astype(np.float32))

        # np.savetxt(dist_file, dist_matrix, delimiter=",", fmt="%d")

    ############################################################################

    def create_connection_matrix(self):

        """ Creates the connection matrix from the synapse matrix data. """

        num_neurons = self.sl.data["num_neurons"]

        con_mat = np.zeros((num_neurons, num_neurons), dtype=int)
        cnt = 0
        pre, post = 0, 0

        for syn_chunk in self.sl.synapse_iterator(data_type="synapses"):
            for syn in syn_chunk:
                p1 = syn[0]
                p2 = syn[1]

                if p1 == pre and p2 == post:
                    cnt += 1
                else:
                    con_mat[pre, post] += cnt
                    pre = p1
                    post = p2
                    cnt = 1

        con_mat[pre, post] += cnt
        cnt = 0

        assert np.sum(np.sum(con_mat)) == self.sl.data["num_synapses"], \
            "Synapse numbers in connection matrix does not match"

        return con_mat

    ############################################################################

    def plot_comparison(self, include_neuron_types=None):

        distance_matrix = self.sl.create_distance_matrix()
        connection_matrix = self.create_connection_matrix()
        axon_dend_dist = self.create_axon_dend_distance_matrix()

        if include_neuron_types:
            neuron_types = self.sl.get_neuron_types()

            include_neuron = np.array([n in include_neuron_types for n in neuron_types], dtype=bool)

            distance_matrix = distance_matrix[include_neuron, :][:, include_neuron]
            connection_matrix = connection_matrix[include_neuron, :][:, include_neuron]
            axon_dend_dist = axon_dend_dist[include_neuron, :][:, include_neuron]

        import matplotlib.pyplot as plt

        plt.figure()
        plt.scatter(distance_matrix[:]*1e6, connection_matrix[:])
        plt.xlabel("Soma-soma distance (mum)")
        plt.ylabel("# Synapses")
        plt.show()

        plt.figure()
        plt.scatter(distance_matrix[:]*1e6, axon_dend_dist[:]*1e6)
        plt.xlabel("Soma-soma distance (mum)")
        plt.ylabel("Axon-dend path distance (mum)")
        plt.show()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Export connection matrix to CSV file")
    parser.add_argument("inFile", help="Snudda HDF5 file with network")
    parser.add_argument("outFile", help="CSV output file")
    parser.add_argument("--full", action="store_false", dest="sparse")
    args = parser.parse_args()

    secm = SnuddaExportConnectionMatrix(args.inFile, args.outFile, save_sparse=args.sparse)
