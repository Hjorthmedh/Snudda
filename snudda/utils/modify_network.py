#!/usr/bin/env python3

from snudda.utils.load import SnuddaLoad
import h5py
import numpy as np
import sys


class SnuddaModifyNetwork:

    """ Modify a network in different ways. """

    def __init__(self, network_file, verbose=False):

        self.snudda_load = SnuddaLoad(network_file=network_file, load_synapses=False, verbose=verbose)
        self.in_file = self.snudda_load.hdf5_file
        self.out_file = None
        self.h5libver = "latest"
        self.h5driver = "sec2"

        self.keep_neuron_id = None
        self.removed_connection_type = None

        self.reset_network()

    def reset_network(self):

        """ Mars all neurons to be kept. """

        self.keep_neuron_id = set(self.in_file["network/neurons/neuronID"][:])
        self.removed_connection_type = []

    def remove_neuron_id(self, neuron_id):

        print(f"Removing neuron_id={set(neuron_id)}")
        self.keep_neuron_id = self.keep_neuron_id - set(neuron_id)

    def remove_neuron_type(self, neuron_type):

        remove_cell_id = self.snudda_load.get_cell_id_of_type(neuron_type=neuron_type)

        if len(remove_cell_id) > 0:
            print(f"Marking {neuron_type} ({len(remove_cell_id)}) for removal")
        else:
            available_neuron_types = sorted(list(set([x["type"] for x in self.snudda_load.data["neurons"]])))
            print(f"No {neuron_type} found in network. Available types are {', '.join(available_neuron_types)}")

        self.keep_neuron_id = self.keep_neuron_id - set(remove_cell_id)

    def remove_neuron_name(self, neuron_name):
        remove_cell_id = self.snudda_load.get_cell_id_with_name(neuron_name=neuron_name)

        if len(remove_cell_id) > 0:
            print(f"Marking {neuron_name} ({len(remove_cell_id)}) for removal")
        else:
            available_neuron_names = sorted(list(set([x["name"] for x in self.snudda_load.data["neurons"]])))
            print(f"No {neuron_name} found in network. Available types are {', '.join(available_neuron_names)}")

        self.keep_neuron_id = self.keep_neuron_id - set(remove_cell_id)

    def remove_connection(self, pre_neuron_type, post_neuron_type):

        available_neuron_types = sorted(list(set([x["type"] for x in self.snudda_load.data["neurons"]])))
        if pre_neuron_type not in available_neuron_types or post_neuron_type not in available_neuron_types:
            print(f"ERROR: Bad connection type {pre_neuron_type},{post_neuron_type}\n"
                  f"Available neuron types: {available_neuron_types}")
            return

        print(f"Marking {pre_neuron_type}, {post_neuron_type} synapses for removal.")
        self.removed_connection_type.append((pre_neuron_type, post_neuron_type))

    def filter_synapses(self, data_type):

        """ Filters synapses, data_type is either 'synapses' or 'gapJunctions' """

        synapse_data = self.in_file[f"network/{data_type}"]

        keep_flag = np.zeros((synapse_data.shape[0],), dtype=bool)

        neuron_types = [n["type"] for n in self.snudda_load.data["neurons"]]

        prev_source = None
        prev_dest = None
        prev_status = None

        for idx, row in enumerate(synapse_data):
            if row[0] == prev_source and row[1] == prev_dest:
                keep_flag[idx] = prev_status
            else:
                prev_source = row[0]
                prev_dest = row[1]

                if row[0] in self.keep_neuron_id and row[1] in self.keep_neuron_id:

                    row_status = 1

                    if self.removed_connection_type:
                        for con_type in self.removed_connection_type:
                            if neuron_types[row[0]] == con_type[0] and neuron_types[row[1]] == con_type[1]:
                                row_status = 0
                                break

                    keep_flag[idx] = row_status
                    prev_status = row_status
                else:
                    prev_status = 0

        return keep_flag

    def write_network(self, out_file_name=None):

        if not out_file_name:
            out_file_name = f"{self.in_file.filename}-modified.hdf5"

        assert out_file_name != self.in_file.filename, f"In and out file must be different."

        print(f"Writing to {out_file_name}")
        self.out_file = h5py.File(out_file_name, "w", libver=self.h5libver, driver=self.h5driver)

        if "config" in self.out_file:
            self.in_file.copy("config", self.out_file)

        self.in_file.copy("meta", self.out_file)

        if "morphologies" in self.in_file:
            print("Copying morphologies")
            self.in_file.copy("morphologies", self.out_file)

        soma_keep_id = list(self.keep_neuron_id)
        num_soma_keep = len(soma_keep_id)

        print(f"Keeping {num_soma_keep} neurons.")

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
                sys.exit(-1)

            if var_name == "neuronID":
                # We need to remap
                neuron_group.create_dataset(var_name, data_shape, data.dtype,
                                            [remap_id[data[x]] for x in soma_keep_id],
                                            compression=data.compression)

                # Double check that it is OK, should be in order after
                assert (np.diff(neuron_group["neuronID"][()]) == 1).all(), "Problem with neuron remapping!"

            else:
                try:
                    neuron_group.create_dataset(var_name, data_shape, data.dtype,
                                                [data[x] for x in soma_keep_id],
                                                compression=data.compression)
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    sys.exit(-1)

        if "synapses" in self.in_file["network"]:

            # Next deal with synapses
            keep_syn_flag = self.filter_synapses(data_type="synapses")

            # Lastly deal with gap junctions
            keep_gj_flag = self.filter_synapses(data_type="gapJunctions")

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

            network_group.create_dataset("synapses", dtype=np.int32, shape=(num_syn, syn_mat.shape[1]),
                                         chunks=syn_mat.chunks, maxshape=(None, syn_mat.shape[1]),
                                         compression=syn_mat.compression)

            network_group.create_dataset("gapJunctions", dtype=np.int32, shape=(num_gj, gj_mat.shape[1]),
                                         chunks=gj_mat.chunks, maxshape=(None, gj_mat.shape[1]),
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

        else:
            print("No synapses found (assuming this was a save file with only position information).")


def snudda_modify_network_cli():

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Modify connections in network.")
    parser.add_argument("original_network", type=str, help="Input network hdf5 file")
    parser.add_argument("output_network", type=str, help="Output network hdf5 file", default=None)
    parser.add_argument("--remove_neuron_type", type=str, help="Neuron type to remove", default=None)
    parser.add_argument("--remove_neuron_name", type=str, help="Neuron name to remove", default=None)
    parser.add_argument("--remove_neuron_id", type=str, help="Neuron ID to remove (e.g. 4,5,6)", default=None)
    parser.add_argument("--remove_connection", type=str, help="Connection to remove (e.g. 'dSPN','iSPN'", default=None)
    args = parser.parse_args()

    mod_network = SnuddaModifyNetwork(network_file=args.original_network)

    if args.remove_neuron_type:
        mod_network.remove_neuron_type(neuron_type=args.remove_neuron_type)

    if args.remove_neuron_name:
        mod_network.remove_neuron_name(neuron_name=args.remove_neuron_name)

    if args.remove_neuron_id:
        neuron_id = [int(x) for x in args.remove_neuron_id.split(",")]
        mod_network.remove_neuron_id(neuron_id=neuron_id)

    if args.remove_connection:
        assert args.remove_connection.count(",") == 1, "Format is --remove_connection pre_neuron_type,post_neuron_type"

        pre_type, post_type = args.remove_connection.split(",")
        mod_network.remove_connection(pre_neuron_type=pre_type, post_neuron_type=post_type)

    mod_network.write_network(out_file_name=args.output_network)


if __name__ == "__main__":
    snudda_modify_network_cli()
