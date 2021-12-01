#!/usr/bin/env python3

import numpy as np

from snudda.utils.export_connection_matrix import SnuddaExportConnectionMatrix


class SnuddaExportErodedConnectionMatrix(SnuddaExportConnectionMatrix):

    def __init__(self, in_file, out_file, fraction_kept=1.0, save_sparse=True, permute=False):

        super().__init__(in_file=in_file, out_file=out_file, save_sparse=save_sparse, save_on_init=False)

        self.con_mat = self.erode(mat=self.con_mat, fraction_kept=fraction_kept)

        if permute:
            self.con_mat = self.permute_all(self.con_mat)

        self.save()

    @staticmethod
    def erode(mat, fraction_kept=1.0):

        erode_mat = np.random.uniform(size=mat.shape) <= fraction_kept
        mat = np.multiply(erode_mat, mat)

        return mat

    @staticmethod
    def permute_all(mat):

        row_permute = np.random.permutation(mat.shape[0])
        col_permute = np.random.permutation(mat.shape[1])

        return mat[row_permute, :][:, col_permute]


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Export ERODED connection matrix to CSV file")
    parser.add_argument("inFile", help="Snudda HDF5 file with network")
    parser.add_argument("outFile", help="CSV output file")
    parser.add_argument("fraction", help="Fraction of connections kept", type=float)
    parser.add_argument("--full", action="store_false", dest="sparse")
    parser.add_argument("--permute", action="store_true")
    args = parser.parse_args()

    s = SnuddaExportErodedConnectionMatrix(args.inFile, args.outFile, save_sparse=args.sparse,
                                           fraction_kept=args.fraction, permute=args.permute)
