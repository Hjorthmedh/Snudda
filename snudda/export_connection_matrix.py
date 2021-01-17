import numpy as np
from snudda.load import SnuddaLoad

class SnuddaExportConnectionMatrix(object):

    def __init__(self, in_file, out_file):

        self.sl = SnuddaLoad(in_file)

        self.outFile = out_file
        self.out_file_meta = out_file + "-meta"

        data = self.sl.data

        con_mat = self.create_con_mat()
        neuron_type = [x["type"] for x in data["neurons"]]
        pos = data["neuronPositions"]

        print("Writing " + self.outFile + " (row = src, column=dest)")
        np.savetxt(self.outFile, con_mat, delimiter=",", fmt="%d")

        print("Writing " + self.out_file_meta)
        with open(self.out_file_meta, "w") as f_out_meta:
            for i, (nt, p) in enumerate(zip(neuron_type, pos)):
                s = "%d,%s,%f,%f,%f\n" % (i, nt, p[0], p[1], p[2])
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
    args = parser.parse_args()

    secm = SnuddaExportConnectionMatrix(args.inFile, args.outFile)
