#!/usr/bin/env python3

import numpy as np

from snudda.utils.export_connection_matrix import SnuddaExportConnectionMatrix


class SnuddaExportErodedConnectionMatrix(SnuddaExportConnectionMatrix):

    def __init__(self, in_file, out_file, fraction_kept=1.0, save_sparse=True, permute=False, permute_type="types"):
        super().__init__(in_file=in_file, out_file=out_file, save_sparse=save_sparse, save_on_init=False)

        self.con_mat = self.erode(mat=self.con_mat, fraction_kept=fraction_kept)

        if permute:
            if permute_type == "all":
                self.con_mat = self.permute_all(self.con_mat)
            elif permute_type == "types":
                self.con_mat = self.permute_within_neuron_types(self.con_mat)
            elif permute_type == "noself":
                self.con_mat = self.permute_all_no_self_connections(self.con_mat)
            else:
                raise ValueError(f"Unknown permutation type {permute_type}, please use 'all', 'types', 'noself'")

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

    def permute_all_no_self_connections(self, mat):

        s = mat.shape
        assert s[0] == s[1] and len(s) == 2

        non_diagonal_idx = np.where(np.diag(np.ones((s[0],))).flatten() == 0)
        permuted_idx = np.random.permutation(non_diagonal_idx)

        new_mat = mat.flatten()

        new_mat[non_diagonal_idx] = mat.flatten()[permuted_idx]
        new_mat = new_mat.reshape(mat.shape)

        return new_mat

    def get_id_of_all_types(self):

        idx = []
        for nt in set(self.sl.get_neuron_types()):
            idx.append(self.sl.get_neuron_id_of_type(neuron_type=nt))

        return idx

    def permute_within_neuron_types(self, mat, list_of_type_idx=None):

        """ Args:
                mat : Connection matrix
                list_of_type_idx : [(1,2,3,4), (5,6,7,8), (9,10)] if neurons of 1,2,3,4 are same type, 5,6,7,8 same
                                   and 9, 10 are of same type.
        """

        if list_of_type_idx is None:
            list_of_type_idx = self.get_id_of_all_types()

        all_idx = np.concatenate(list_of_type_idx)
        all_idx_sorted = np.sort(all_idx)
        assert (np.diff(all_idx_sorted) == 1).all()
        assert all_idx_sorted[0] == 0 and len(all_idx_sorted) == mat.shape[0] == mat.shape[1]

        new_mat = np.full(mat.shape, np.nan).flatten()

        for type_ctr_pre, pre_idx in enumerate(list_of_type_idx):
            for type_ctr_post, post_idx in enumerate(list_of_type_idx):

                flat_idx = np.ravel_multi_index(np.ix_(pre_idx, post_idx), mat.shape).flatten()

                if type_ctr_pre == type_ctr_post:
                    new_mat[flat_idx] = self.permute_all_no_self_connections(mat[pre_idx, :][:, post_idx]).flatten()
                else:
                    new_mat[flat_idx] = self.permute_all(mat[pre_idx, :][:, post_idx]).flatten()

        new_mat = new_mat.reshape(mat.shape)

        assert np.sum(np.isnan(new_mat)) == 0, f"Internal error, not all indexes are given"

        return new_mat


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Export ERODED connection matrix to CSV file")
    parser.add_argument("in_file", help="Snudda HDF5 file with network")
    parser.add_argument("out_file", help="CSV output file")
    parser.add_argument("fraction", help="Fraction of connections kept", type=float)
    parser.add_argument("--full", action="store_false", dest="sparse")
    parser.add_argument("--permute", action="store_true")
    parser.add_argument("--permutation_type",
                        help="Permutation type: 'all', 'types' (preserve types), "
                             "'noself' (dont preserve types, but avoid self connections",
                        default="types")
    args = parser.parse_args()

    s = SnuddaExportErodedConnectionMatrix(args.in_file, args.out_file, save_sparse=args.sparse,
                                           fraction_kept=args.fraction, permute=args.permute,
                                           permute_type=args.permutation_type)
