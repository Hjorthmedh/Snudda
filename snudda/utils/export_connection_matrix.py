import numpy as np
from snudda.utils.load import SnuddaLoad

class SnuddaExportConnectionMatrix(object):

    def __init__(self, in_file, out_file, save_sparse=True):

        self.sl = SnuddaLoad(in_file)

        self.outFile = out_file
        self.out_file_meta = out_file + "-meta"

        data = self.sl.data

        con_mat = self.create_con_mat()
        neuron_type = [x["type"] for x in data["neurons"]]
        neuron_name = [x["name"] for x in data["neurons"]]
        morph_file = [data["morph"][nn]["location"] for nn in neuron_name]
        pos = data["neuronPositions"]

        print("Writing " + self.outFile + " (row = src, column=dest)")
        if save_sparse:
            x_pos, y_pos = np.where(con_mat)
            sparse_data = np.zeros((len(x_pos), 3), dtype=int)
            for idx, (x, y) in enumerate(zip(x_pos, y_pos)):
                sparse_data[idx, :] = [x, y, con_mat[x, y]]

            np.savetxt(self.outFile, sparse_data, delimiter=",", fmt="%d")

            # Test to verify
            for row in sparse_data:
                assert con_mat[row[0], row[1]] == row[2]
        else:
            np.savetxt(self.outFile, con_mat, delimiter=",", fmt="%d")

        print("Writing " + self.out_file_meta)
        with open(self.out_file_meta, "w") as f_out_meta:
            for i, (nt, nn, p, mf) in enumerate(zip(neuron_type, neuron_name, pos, morph_file)):
                s = "%d,%s,%s,%f,%f,%f,%s\n" % (i, nt, nn, p[0], p[1], p[2], mf)
                f_out_meta.write(s)
            f_out_meta.close()

        # import pdb
        # pdb.set_trace()

    ############################################################################

    def create_con_mat(self):

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
    parser.add_argument("--sparse", action="store_true", default=False)
    args = parser.parse_args()

    secm = SnuddaExportConnectionMatrix(args.inFile, args.outFile, save_sparse=args.sparse)
