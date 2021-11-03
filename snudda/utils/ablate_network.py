#!/usr/bin/env python3
import os.path
from collections import OrderedDict

from snudda.utils.load import SnuddaLoad
import h5py
import numpy as np
import sys
import json


class SnuddaAblateNetwork:

    """ Ablate neurons or synapses from a network.

        For an example config see Snudda/examples/config/example-ablation-config.json

          {
            "ablate_neurons": ["FS",
		       ["iSPN", 0.3],
		       [1,2,3,4]],
            "ablate_synapses" : [["dSPN", "iSPN"],
			 ["iSPN", "dSPN", 0.1]]
          }

        Here "ablate_neurons" is a list of neurons to remove. The list can contain neuron types (e.g. "FS"),
        or neuron type and removal probability (e.g. ["iSPN", 0.3], or neuron ID to remove (e.g. [1,2,3,4]).

        With "ablate_synapses" the user can specify what synapses to remove (e.g. ["dSPN", "iSPN"]] to remove
        all synapses between dSPN and iSPN). Alternatively the removal probability can also be specified
        (e.g. ["iSPN", "dSPN", 0.1] where 10% of connections between those neuron types are removed).


    """

    def __init__(self, network_file, verbose=False):

        self.snudda_load = SnuddaLoad(network_file=network_file, load_synapses=False, verbose=verbose)
        self.in_file = self.snudda_load.hdf5_file
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

    def remove_neuron_type(self, neuron_type, p_remove=1):

        """ Remove neuron of type neuron_type with probability p_remove (default 1)"""

        remove_cell_id = self.snudda_load.get_neuron_id_of_type(neuron_type=neuron_type)
        remove_flag = np.random.uniform(size=(len(remove_cell_id),)) <= p_remove
        remove_cell_id = remove_cell_id[remove_flag]

        if len(remove_flag) > 0:
            print(f"Marking {neuron_type} ({np.sum(remove_flag)} out of {len(remove_flag)}) for removal (P={p_remove})")
        else:
            available_neuron_types = sorted(list(set([x["type"] for x in self.snudda_load.data["neurons"]])))
            print(f"No {neuron_type} found in network. Available types are {', '.join(available_neuron_types)}")

        self.keep_neuron_id = self.keep_neuron_id - set(remove_cell_id)

    def remove_neuron_name(self, neuron_name):

        """ Remove neuron with name neuron_name """

        remove_cell_id = self.snudda_load.get_neuron_id_with_name(neuron_name=neuron_name)

        if len(remove_cell_id) > 0:
            print(f"Marking {neuron_name} ({len(remove_cell_id)}) for removal")
        else:
            available_neuron_names = sorted(list(set([x["name"] for x in self.snudda_load.data["neurons"]])))
            print(f"No {neuron_name} found in network. Available types are {', '.join(available_neuron_names)}")

        self.keep_neuron_id = self.keep_neuron_id - set(remove_cell_id)

    def remove_connection(self, pre_neuron_type, post_neuron_type, p_remove=1):

        """ Removes connections between specified neuron types. """

        available_neuron_types = sorted(list(set([x["type"] for x in self.snudda_load.data["neurons"]])))
        if pre_neuron_type not in available_neuron_types or post_neuron_type not in available_neuron_types:
            print(f"ERROR: Bad connection type {pre_neuron_type},{post_neuron_type}\n"
                  f"Available neuron types: {available_neuron_types}")
            return

        print(f"Marking {pre_neuron_type}, {post_neuron_type} synapses for removal (P={p_remove}).")
        self.removed_connection_type.append((pre_neuron_type, post_neuron_type, p_remove))

    def filter_synapses(self, data_type):

        """ Filters synapses, data_type is either 'synapses' or 'gapJunctions' """

        synapse_data = self.in_file[f"network/{data_type}"]

        keep_flag = np.zeros((synapse_data.shape[0],), dtype=bool)

        neuron_types = [n["type"] for n in self.snudda_load.data["neurons"]]

        prev_source = None
        prev_dest = None
        prev_status = None

        n_original_synapses = synapse_data.shape[0]

        for idx, (pre_id, post_id) in enumerate(zip(synapse_data[:,0], synapse_data[:,1])):

            if idx % 10000000 == 0:
                print(f"{idx}/{n_original_synapses} synapses processed")

            if pre_id == prev_source and post_id == prev_dest:
                keep_flag[idx] = prev_status
            else:
                prev_source = pre_id
                prev_dest = post_id

                if pre_id in self.keep_neuron_id and post_id in self.keep_neuron_id:

                    row_status = 1

                    if self.removed_connection_type:
                        for con_type in self.removed_connection_type:
                            if neuron_types[pre_id] == con_type[0] and neuron_types[post_id] == con_type[1]:
                                if np.random.uniform() <= con_type[2]:
                                    # All synapses between a given neuron pair is removed together
                                    row_status = 0
                                    break

                    keep_flag[idx] = row_status
                    prev_status = row_status
                else:
                    prev_status = 0

        print(f"{n_original_synapses}/{n_original_synapses} synapses processed")
        print("Filtering done.")

        return keep_flag

    def write_network(self, out_file_name=None):

        """ Write network to hdf5 file: output_file_name """

        if not out_file_name:
            out_file_name = f"{self.in_file.filename}-modified.hdf5"

        assert out_file_name != self.in_file.filename, f"In and out file must be different."

        print(f"Writing to {out_file_name}")
        out_file = h5py.File(out_file_name, "w", libver=self.h5libver, driver=self.h5driver)

        if "config" in out_file:
            self.in_file.copy("config", out_file)

        self.in_file.copy("meta", out_file)

        if "morphologies" in self.in_file:
            print("Copying morphologies")
            self.in_file.copy("morphologies", out_file)

        soma_keep_id = list(self.keep_neuron_id)
        num_soma_keep = len(soma_keep_id)

        print(f"Keeping {num_soma_keep} neurons.")

        # We need to remap neuronID in the synapses and gap junction matrix
        # remap_id = dict([])
        # Try using a np array for lookup instead of dict, faster?
        remap_id = np.full((len(self.snudda_load.data["neurons"]),), np.nan, dtype=int)
        for new_id, old_id in enumerate(soma_keep_id):
            remap_id[old_id] = new_id

        network_group = out_file.create_group("network")
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

            n_synapses = syn_mat.shape[0]

            temp_syn_mat = np.zeros((num_syn, syn_mat.shape[1]), dtype=np.int32)
            temp_gj_mat = np.zeros((num_gj, gj_mat.shape[1]), dtype=np.int32)

            syn_keep_idx = np.where(keep_syn_flag)[0]
            for idx, row_idx in enumerate(syn_keep_idx):

                if idx % 50000 == 0:
                    print(f"{idx} / {num_syn} synapse rows parsed")

                # We need to remap the neuronID if some neurons have been removed!!
                row = syn_mat[row_idx, :]
                row[0] = remap_id[row[0]]
                row[1] = remap_id[row[1]]
                temp_syn_mat[idx, :] = row

            network_group.create_dataset("synapses",
                                         data=temp_syn_mat,
                                         dtype=np.int32, shape=(num_syn, syn_mat.shape[1]),
                                         chunks=syn_mat.chunks, maxshape=(None, syn_mat.shape[1]),
                                         compression=syn_mat.compression)

            print(f"{n_synapses} / {num_syn} synapse rows parsed")
            print("Synapse matrix written.")

            print(f"Keeping {num_syn} synapses (out of {syn_mat.shape[0]})")

            n_gj = gj_mat.shape[0]

            for idx, row_idx in enumerate(np.where(keep_gj_flag)[0]):

                if idx % 50000 == 0:
                    print(f"{idx} / {num_gj} gap junction rows parsed")

                # We need to remap the neuronID if some neurons have been removed!!
                row = gj_mat[row_idx, :]
                row[0] = remap_id[row[0]]
                row[1] = remap_id[row[1]]
                temp_gj_mat[idx, :] = row

            network_group.create_dataset("gapJunctions",
                                         data=temp_gj_mat,
                                         dtype=np.int32, shape=(num_gj, gj_mat.shape[1]),
                                         chunks=gj_mat.chunks, maxshape=(None, gj_mat.shape[1]),
                                         compression=gj_mat.compression)

            print(f"{num_gj} / {num_gj} gap junction rows parsed")
            print("Gap junction matrix written.")
            print(f"Keeping {num_gj}  gap junctions (out of {gj_mat.shape[0]})")

        else:
            print("No synapses found (assuming this was a save file with only position information).")

        out_file.close()


def snudda_ablate_network_cli():

    from argparse import ArgumentParser

    # TODO: Fix so ablation can be specified using a json file for more complex ablations

    parser = ArgumentParser(description="Modify connections in network.")
    parser.add_argument("original_network", type=str, help="Input network hdf5 file")
    parser.add_argument("output_network", type=str, help="Output network hdf5 file", default=None)
    parser.add_argument("--config", type=str, help="Ablation config file", default=None)
    parser.add_argument("--remove_neuron_type", type=str, help="Neuron type to remove", default=None)
    parser.add_argument("--p_remove_neuron_type", type=float, help="Probability to remove neuron of type", default=1.0)
    parser.add_argument("--remove_neuron_name", type=str, help="Neuron name to remove", default=None)
    parser.add_argument("--remove_neuron_id", type=str, help="Neuron ID to remove (e.g. 4,5,6)", default=None)
    parser.add_argument("--remove_connection", type=str, help="Connection to remove (e.g. 'dSPN','iSPN'", default=None)
    parser.add_argument("--p_remove_connection", type=float, help="Probability to remove connection", default=1.0)
    args = parser.parse_args()

    mod_network = SnuddaAblateNetwork(network_file=args.original_network)

    if args.config:
        if not os.path.isfile(args.config):
            print(f"There is no config file {args.config}")
            exit(-1)

        with open(args.config) as f:
            config_data = json.load(f, object_pairs_hook=OrderedDict)

        if "ablate_neurons" in config_data:
            for ablate in config_data["ablate_neurons"]:
                if type(ablate[0]) == int:
                    # Ablate neurons with number
                    mod_network.remove_neuron_id(ablate)
                elif "_" in ablate[0]:
                    # Ablate neuron with name
                    mod_network.remove_neuron_name(neuron_name=ablate[0])
                else:
                    neuron_type = ablate[0]
                    if type(ablate) == list and len(ablate) > 1:
                        ablate_p = float(ablate[1])
                    else:
                        ablate_p = 1

                    mod_network.remove_neuron_type(neuron_type=neuron_type, p_remove=ablate_p)

            for con in config_data["ablate_synapses"]:
                pre_type = con[0]
                post_type = con[1]
                if len(con) > 2:
                    p_remove = con[2]
                else:
                    p_remove = 1
                mod_network.remove_connection(pre_neuron_type=pre_type, post_neuron_type=post_type, p_remove=p_remove)

    if args.remove_neuron_type:
        mod_network.remove_neuron_type(neuron_type=args.remove_neuron_type,
                                       p_remove=args.p_remove_neuron_type)

    if args.remove_neuron_name:
        mod_network.remove_neuron_name(neuron_name=args.remove_neuron_name)

    if args.remove_neuron_id:
        neuron_id = [int(x) for x in args.remove_neuron_id.split(",")]
        mod_network.remove_neuron_id(neuron_id=neuron_id)

    if args.remove_connection:

        if args.remove_connection.count(";") > 0:

            for c in args.remove_connection.split(";"):
                
                assert args.remove_connection.count(",") == 1, "Format is --remove_connection pre_neuron_type,post_neuron_type"

                pre_type, post_type = c.split(",")
                mod_network.remove_connection(pre_neuron_type=pre_type, post_neuron_type=post_type,
                                              p_remove=args.p_remove_connection)


        else:
            assert args.remove_connection.count(",") == 1, "Format is --remove_connection pre_neuron_type,post_neuron_type"

            pre_type, post_type = args.remove_connection.split(",")
            mod_network.remove_connection(pre_neuron_type=pre_type, post_neuron_type=post_type,
                                          p_remove=args.p_remove_connection)

    mod_network.write_network(out_file_name=args.output_network)


if __name__ == "__main__":
    snudda_ablate_network_cli()
