# This code allows the user to remove a portion of the network by cutting a
# slice of the data.

# 0. Load the network
# 1. Define the cut plane
# 2. Find if any of the neurons are above the cut plane, remove those
# 3. Remove any synapses above the cut plane
# 4. Write a new hdf5
#

import os
import h5py
import numpy as np
import time
import numexpr

# TODO: Check that no population unit list needs to be remapped also, this is old code that was just refactored...


class SnuddaCut(object):

    def __init__(self, network_file, cut_equation="z>0",
                 out_file_name=None,
                 plot_only=False, show_plot=True):

        self.cut_equation = cut_equation

        self.h5libver = "latest"
        self.h5driver = "sec2"
        self.in_file = None
        self.out_file = None

        if out_file_name is None:
            base_path = os.path.dirname(network_file)
            self.out_file_name = os.path.join(base_path, "network-cut-slice.hdf5")
        else:
            self.out_file_name = out_file_name

        # We create a lambda expression, to avoid calling eval each time
        self.cut_equation_lambda = eval(f"lambda x,y,z : numexpr.evaluate({cut_equation})")

        self.open_input_file(network_file)

        if plot_only:
            self.plot_cut(include_synapses=False, include_gap_junctions=False)
            exit(0)

        self.setup_output_file(self.out_file_name)
        self.write_cut_slice(self.cut_equation_lambda)

        self.plot_cut(include_synapses=True, include_gap_junctions=True,
                      show_plot=show_plot)

        if False:
            self.in_file.close()
            self.out_file.close()

    ############################################################################

    def write_cut_slice(self, cut_equation_lambda):

        # Remove the neurons from the data
        soma_keep_flag = self.somas_inside(cut_equation_lambda)
        soma_keep_id = np.where(soma_keep_flag)[0]
        soma_remove_id = np.where(soma_keep_flag == False)[0]
        num_soma_keep = np.sum(soma_keep_flag)

        if num_soma_keep == 0:
            print("No somas left, aborting!")
            exit(-1)

        print(f"Keeping {num_soma_keep} out of {len(soma_keep_flag)}" 
              "neurons (the others have soma outside of cut plane)")

        # We need to remap neuronID in the synapses and gap junction matrix
        remap_id = dict([])
        for new_id, old_id in enumerate(soma_keep_id):
            remap_id[old_id] = new_id

        network_group = self.out_file.create_group("network")
        neuron_group = network_group.create_group("neurons")

        for var_name in self.in_file["network/neurons"]:

            data = self.in_file[f"network/neurons/{var_name}"]

            if len(data.shape) == 0:
                # Scalar data, just copy
                self.in_file.copy(f"network/neurons/{var_name}", neuron_group)
                continue

            elif len(data.shape) == 1:
                # 1D data, we only keep nSomaKeep of them
                data_shape = (num_soma_keep,)
            elif len(data.shape) == 2:
                # 2D data, need to make sure to maintain dimensions
                data_shape = (num_soma_keep, data.shape[1])
            else:
                print("writeCutSlice: Only handle 0D, 1D and 2D data, update code!")
                exit(-1)

            if var_name == "neuronID":
                # We need to remap
                neuron_group.create_dataset(var_name, data_shape, data.dtype,
                                            [remap_id[data[x]] for x in soma_keep_id],
                                            compression=data.compression)

                # Double check that it is OK, should be in order after
                assert (np.diff(neuron_group["neuronID"][()]) == 1).all(), \
                    "Problem with neuron remapping!"

            else:
                try:
                    neuron_group.create_dataset(var_name, data_shape, data.dtype,
                                                [data[x] for x in soma_keep_id],
                                                compression=data.compression)
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    exit(-1)

        # Next deal with synapses
        keep_syn_flag = self.synapses_inside(cut_equation_lambda, data_type="synapses")
        keep_syn_flag = self.filter_neurons_synapses(soma_remove_id,
                                                     keep_flag=keep_syn_flag,
                                                     data_type="synapses")

        # Lastly deal with gap junctions
        keep_gj_flag = self.synapses_inside(cut_equation_lambda, data_type="gapJunctions")
        keep_gj_flag = self.filter_neurons_synapses(soma_remove_id,
                                                    keep_flag=keep_gj_flag,
                                                    data_type="gapJunctions")

        num_syn = np.sum(keep_syn_flag)
        num_synapses = np.zeros((1,), dtype=np.uint64) + num_syn

        num_gj = np.sum(keep_gj_flag)
        num_gap_junctions = np.zeros((1,), dtype=np.uint64) + num_gj

        network_group.create_dataset("nSynapses", data=num_synapses, dtype=np.uint64)
        network_group.create_dataset("nGapJunctions", data=num_gap_junctions,
                                     dtype=np.uint64)

        # TODO: !!!! might need to handle chunk size differently based on size...

        syn_mat = self.in_file["network/synapses"]
        gj_mat = self.in_file["network/gapJunctions"]

        print("Copying synapses and gap junctions")

        network_group.create_dataset("synapses", dtype=np.int32, shape=(num_syn, 13),
                                     chunks=syn_mat.chunks, maxshape=(None, 13),
                                     compression=syn_mat.compression)

        network_group.create_dataset("gapJunctions", dtype=np.int32, shape=(num_gj, 11),
                                     chunks=gj_mat.chunks, maxshape=(None, 11),
                                     compression=gj_mat.compression)

        for idx, row_idx in enumerate(np.where(keep_syn_flag)[0]):
            # We need to remap the neuronID if some neurons have been removed!!
            row = syn_mat[row_idx, :]
            row[0] = remap_id[row[0]]
            row[1] = remap_id[row[1]]
            network_group["synapses"][idx, :] = row

        print(f"Keeping {num_syn} synapses (out of {syn_mat.shape[0]})")

        for idx, row_idx in enumerate(np.where(keep_gj_flag)[0]):
            # We need to remap the neuronID if some neurons have been removed!!
            row = gj_mat[row_idx, :]
            row[0] = remap_id[row[0]]
            row[1] = remap_id[row[1]]
            network_group["gapJunctions"][idx, :] = row

        print(f"Keeping {num_gj}  gap junctions (out of {gj_mat.shape[0]})")

    ############################################################################

    # Tells which somas are inside the cut equation

    def somas_inside(self, cut_equation_lambda):

        pos = self.in_file["network/neurons/position"][()]
        inside_flag = np.array([cut_equation_lambda(x, y, z) for x, y, z in pos], dtype=bool)

        return inside_flag

    ############################################################################

    def synapses_inside(self, cut_equation_lambda, data_type="synapses"):

        voxel_size = self.in_file["meta/voxelSize"][()]
        sim_origo = self.in_file["meta/simulationOrigo"][()]

        if data_type == "synapses":
            pos = self.in_file["network/synapses"][:, 2:5] * voxel_size + sim_origo
        elif data_type == "gapJunctions":
            pos = self.in_file["network/gapJunctions"][:, 6:9] * voxel_size + sim_origo
        else:
            print(f"filterNeuronsSynapses: Unknown data_type: {data_type} (valid are 'synapses' or 'gapJunctions'")
            exit(-1)

        inside_flag = np.array([cut_equation_lambda(x, y, z) for x, y, z in pos], dtype=bool)

        return inside_flag

    ############################################################################

    # Returns the row numbers that do not contain the neuronID, ie filters
    # the synapses belonging to neuronID out...

    # dataType = "synapses" or "gapJunctions"

    def filter_neurons_synapses(self, neuron_id, keep_flag=None, data_type="synapses"):

        if data_type == "synapses":
            data_str = "network/synapses"
        elif data_type == "gapJunctions":
            data_str = "network/gapJunctions"
        else:
            print(f"filter_neurons_synapses: Unknown data_type: {data_type} (valid are 'synapses', 'gapJunctions'")
            exit(-1)

        if keep_flag is None:
            keep_flag = np.ones((self.in_file[data_str].shape[0],), dtype=bool)

        src_id = self.in_file[data_str][:, 0]
        dest_id = self.in_file[data_str][:, 1]

        for n_id in neuron_id:
            keep_flag = np.logical_and(keep_flag,
                                       np.logical_and(src_id != n_id, dest_id != n_id))

        return keep_flag

    ############################################################################

    def open_input_file(self, network_file):

        self.in_file = h5py.File(network_file, "r", libver=self.h5libver, driver=self.h5driver)

    ############################################################################

    # This sets up the output file, copies the config and meta data,
    # but does not copy over the neurons, synapses or gap junctions

    def setup_output_file(self, out_file_name):

        print(f"Writing to {out_file_name}")

        self.out_file = h5py.File(out_file_name, "w", libver=self.h5libver, driver=self.h5driver)

        print("Copying 'config' and 'meta'")
        self.in_file.copy("config", self.out_file)
        self.in_file.copy("meta", self.out_file)

        if "morphologies" in self.in_file:
            print("Copying morphologies")
            self.in_file.copy("morphologies", self.out_file)

    ############################################################################

    # This is just used to verify

    def plot_cut(self, include_synapses=True, include_gap_junctions=True, show_plot=True):

        print("Plotting verification figure")

        voxel_size = self.in_file["meta/voxelSize"][()]
        sim_origo = self.in_file["meta/simulationOrigo"][()]

        in_pos = self.in_file["network/neurons/position"][()]
        in_syn = self.in_file["network/synapses"][:, 2:5] * voxel_size + sim_origo
        in_gj = self.in_file["network/gapJunctions"][:, 6:9] * voxel_size + sim_origo

        if self.out_file is not None:
            # Just double check that they match
            assert self.out_file["meta/voxelSize"][()] == voxel_size
            assert (self.out_file["meta/simulationOrigo"][()] == sim_origo).all()

            out_pos = self.out_file["network/neurons/position"][()]
            out_syn = self.out_file["network/synapses"][:, 2:5] * voxel_size + sim_origo
            out_gj = self.out_file["network/gapJunctions"][:, 6:9] * voxel_size + sim_origo

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        a_val = 0.5
        a_val2 = 0.2

        ax.scatter(in_pos[:, 0], in_pos[:, 1], in_pos[:, 2], c='black', s=25, alpha=a_val2)
        if self.out_file is not None:
            ax.scatter(out_pos[:, 0], out_pos[:, 1], out_pos[:, 2], c='red', s=15, alpha=a_val)

        if include_synapses:
            ax.scatter(in_syn[:, 0], in_syn[:, 1], in_syn[:, 2], c='black', s=9, alpha=a_val2)

            if self.out_file is not None:
                ax.scatter(out_syn[:, 0], out_syn[:, 1], out_syn[:, 2], c='red', s=4, alpha=a_val)

        if include_gap_junctions:
            ax.scatter(in_gj[:, 0], in_gj[:, 1], in_gj[:, 2], c='blue', s=6, alpha=a_val2)
            if self.out_file is not None:
                ax.scatter(out_gj[:, 0], out_gj[:, 1], out_gj[:, 2], c='green', s=3, alpha=a_val)

        ax.view_init(elev=-3, azim=-95)

        if show_plot:
            plt.ion()
            plt.show()

        fig_name = f"{self.out_file_name}.pdf"
        print(f"Writing to figure {fig_name}")
        plt.savefig(fig_name, dpi=300)

        if show_plot:
            print("Inspect plot, then quit debug")
            print("The viewing angle might not be good for your try, so leave it interactive")
            import pdb
            pdb.set_trace()
        else:
            print("Waiting 5 seconds")
            time.sleep(5)

    ############################################################################


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Cut a slice from a network")
    parser.add_argument("networkFile", help="Network file (hdf5)", type=str)
    parser.add_argument("cutEquation",
                        help="Equation defining parts left after cut, " \
                             + "e.g. 'z>0' or 'x+y>100e-6' (SI units)",
                        type=str)
    parser.add_argument("--plotOnly", \
                        help="Plots the network without cutting",
                        action="store_true")
    parser.add_argument("--hidePlot", help="Hide plot.", action="store_true")

    args = parser.parse_args()

    if args.hidePlot:
        showPlot = False
    else:
        showPlot = True

    sc = SnuddaCut(network_file=args.networkFile,
                   cut_equation=args.cutEquation,
                   plot_only=args.plotOnly,
                   show_plot=showPlot)
