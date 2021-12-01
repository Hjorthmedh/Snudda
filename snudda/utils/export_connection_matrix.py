#!/usr/bin/env python3

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

        self.outFile = out_file
        self.out_file_meta = f"{out_file}-meta"
        self.save_sparse = save_sparse

        data = self.sl.data

        self.con_mat = self.create_con_mat()
        self.neuron_type = [x["type"] for x in data["neurons"]]
        self.neuron_name = [x["name"] for x in data["neurons"]]
        self.neuron_morph = [x["morphology"] for x in data["neurons"]]

        self.pos = data["neuronPositions"]

        if save_on_init:
            self.save()

    def save(self):

        print(f"Writing {self.outFile} (row = src, column=dest)")
        print(f"Saving {np.sum(np.sum(self.con_mat))} synapses, over {np.count_nonzero(self.con_mat)} coupled pairs.")

        if self.save_sparse:
            x_pos, y_pos = np.where(self.con_mat)
            sparse_data = np.zeros((len(x_pos), 3), dtype=int)
            for idx, (x, y) in enumerate(zip(x_pos, y_pos)):
                sparse_data[idx, :] = [x, y, self.con_mat[x, y]]

            np.savetxt(self.outFile, sparse_data, delimiter=",", fmt="%d")

            # Test to verify
            for row in sparse_data:
                assert self.con_mat[row[0], row[1]] == row[2]
        else:
            np.savetxt(self.outFile, self.con_mat, delimiter=",", fmt="%d")

        print("Writing " + self.out_file_meta)
        with open(self.out_file_meta, "w") as f_out_meta:
            for i, (nt, nn, p, mf) in enumerate(zip(self.neuron_type, self.neuron_name, self.pos, self.neuron_morph)):
                s = "%d,%s,%s,%f,%f,%f,%s\n" % (i, nt, nn, p[0], p[1], p[2], mf)
                f_out_meta.write(s)
            f_out_meta.close()

    ############################################################################

    def create_con_mat(self):

        """ Creates the connection matrix from the synapse matrix data. """

        num_neurons = self.sl.data["nNeurons"]

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

        assert np.sum(np.sum(con_mat)) == self.sl.data["nSynapses"], \
            "Synapse numbers in connection matrix does not match"

        return con_mat


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Export connection matrix to CSV file")
    parser.add_argument("inFile", help="Snudda HDF5 file with network")
    parser.add_argument("outFile", help="CSV output file")
    parser.add_argument("--full", action="store_false", dest="sparse")
    args = parser.parse_args()

    secm = SnuddaExportConnectionMatrix(args.inFile, args.outFile, save_sparse=args.sparse)
