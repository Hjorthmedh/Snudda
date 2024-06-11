# Must have a __init__.py file in module directory

import csv
import json
import os
from collections import OrderedDict

import h5py
import numpy as np


# TODO: This needs to be updated to match the latest SONATA release.

class ConvHurt(object):

    def __init__(self, simulation_structures, base_dir="TEST/", target_simulator="NEST", has_input=False):

        self.base_dir = base_dir
        self.network_dir = os.path.join(base_dir, 'networks')

        self.target_simulator = target_simulator

        self.h5py_libver = "earliest"  # "latest"

        self.setup_directories(base_dir=base_dir)

        self.write_main_config(simulation_structures=simulation_structures,
                               base_dir=base_dir, target_simulator=target_simulator, has_input=has_input)

    ############################################################################

    @staticmethod
    def create_dir(directory):

        if not os.path.exists(directory):
            print(f"Creating directory : {directory}")
            os.makedirs(directory)

    ############################################################################

    def setup_directories(self,
                          base_dir=None):

        if base_dir is None:
            base_dir = "."

        self.create_dir(base_dir)

        dir_list = ['components',
                    os.path.join("components", "biophysical_neuron_dynamics"),
                    os.path.join("components", "hoc_templates"),
                    os.path.join("components", "mechanisms"),
                    os.path.join("components", "morphologies"),
                    os.path.join("components", "point_neuron_dynamics"),
                    os.path.join("components", "synapse_dynamics"),
                    'networks', 'inputs', 'output']

        for x in dir_list:
            self.create_dir(os.path.join(base_dir, x))

    ############################################################################
    #
    # Writing the main configuration file, linking everything together

    def write_main_config(self,
                          base_dir="TEST/",
                          out_file="circuit_config.json",
                          simulation_structures=[],
                          target_simulator="NEURON",
                          has_input=False):

        config = OrderedDict([])

        config["target_simulator"] = target_simulator

        manifest = OrderedDict([("$BASE_DIR", "."),
                                ("$NETWORK_DIR", "$BASE_DIR/networks"),
                                ("$COMPONENT_DIR", "$BASE_DIR/components")])

        config["manifest"] = manifest

        self.base_dir = base_dir
        self.network_dir = os.path.join(base_dir, 'networks')

        components = OrderedDict([("morphologies_dir", os.path.join("$COMPONENT_DIR", "morphologies")),
                                  ("synaptic_models_dir", os.path.join("$COMPONENT_DIR", "synapse_dynamics")),
                                  ("mechanisms_dir", os.path.join("$COMPONENT_DIR", "mechanisms")),
                                  ("biophysical_neuron_models_dir",
                                   os.path.join("$COMPONENT_DIR", "biophysical_neuron_dynamics")),
                                  ("point_neuron_models_dir", os.path.join("$COMPONENT_DIR", "point_neuron_dynamics")),
                                  ("templates_dir", os.path.join("$COMPONENT_DIR", "hoc_templates"))])

        config["components"] = components

        nodes = []
        edges = []
        for ss in simulation_structures:
            node_info = {"nodes_file": os.path.join("$NETWORK_DIR", f"{ss}_nodes.hdf5"),
                         "node_types_file": os.path.join("$NETWORK_DIR", f"{ss}_node_types.csv") }
            nodes.append(node_info)

            edge_info = {"edges_file": os.path.join("$NETWORK_DIR", f"{ss}_edges.hdf5"),
                         "edge_types_file": os.path.join("$NETWORK_DIR", f"{ss}_edge_types.csv") }
            edges.append(edge_info)

            if has_input:
                node_info = {"nodes_file": os.path.join("$NETWORK_DIR", f"{ss}-input_nodes.hdf5"),
                             "node_types_file": os.path.join("$NETWORK_DIR", f"{ss}-input_node_types.csv")}
                nodes.append(node_info)

                edge_info = {"edges_file": os.path.join("$NETWORK_DIR", f"{ss}-input_edges.hdf5"),
                             "edge_types_file": os.path.join("$NETWORK_DIR", f"{ss}-input_edge_types.csv")}
                edges.append(edge_info)

        config["networks"] = OrderedDict([("nodes", nodes), ("edges", edges)])

        with open(os.path.join(base_dir, out_file), 'wt') as f:
            json.dump(config, f, indent=4)

    ############################################################################

    # nodeID - vector with IDs for all neurons
    # nodeTypeID - vector with Type ID for each neuron
    # nodeGroupID - vector with group membership for each neuron
    # nodeGroupIndex - vector with index of each neuron within the given group

    # data = dictionary with "position", "rotation_angle_zaxis", ...
    #        etc that should be stored
    #        names should match what is in the group data
    def write_nodes(self,
                    node_file,
                    population_name,
                    node_id,
                    node_type_id,
                    node_group_id,
                    node_group_index,
                    data,
                    model_type=None,
                    model_template=None,
                    close_file=True):

        if isinstance(node_file, h5py._hl.files.File):
            f = node_file
            n_group = f["nodes"]
        else:
            f = h5py.File(os.path.join(self.network_dir, node_file), 'w', libver=self.h5py_libver)
            self.add_version(f)
            n_group = f.create_group("nodes")

        print(f"Creating nodes/{population_name}")
        nodes_group = n_group.create_group(population_name)

        nodes_group.create_dataset("node_id", data=node_id)
        nodes_group.create_dataset("node_type_id", data=node_type_id)
        nodes_group.create_dataset("node_group_id", data=node_group_id)
        nodes_group.create_dataset("node_group_index", data=node_group_index)

        groups = np.unique(node_group_id)

        for g in groups:
            idx = np.where(node_group_id == g)

            group_group = nodes_group.create_group(str(g))
            # dynGroup = groupGroup.create_group("dynamics_params")

            for data_type in data.keys():
                # This assumes all our data is matrices
                # if there are string data, we need to handle that separately
                if len(data[data_type].shape) == 1:
                    group_group.create_dataset(data_type, data=data[data_type][idx])
                elif len(data[data_type].shape) == 2:
                    group_group.create_dataset(data_type, data=data[data_type][idx, :])
                else:
                    print("Data has too many columns, require 1 or 2 columns max.")

            # These are optional, would overwrite the defaults in the CSV file
            if model_type is not None and model_type != [None]:
                # If this line fails, try to use .encode() on each element in the list
                str_type = f"S{max([len(x) if x is not None else 1 for x in model_type[idx]])}"

                nodes_group.create_dataset("model_type", (len(idx),), str_type,
                                           data=model_type[idx],
                                           compression="gzip")

            if model_template is not None and model_template != [None]:
                assert model_type is not None, "model_type must be set if model_template set"

                str_type2 = f"S{max([len(x) for x in model_template[idx]])}"
                nodes_group.create_dataset("model_template", (len(idx),), str_type2, data=model_template[idx])

            # model_type_id should be added here also

        if close_file:
            f.close()
            f = None

        return f

    ############################################################################

    # node_csv_file is the name of the file that will be stored in the networks dir
    # node_type_id is the node type ID for each row
    # data is a dictionary, where each key corresponds to a column, and it must
    # contain a list of len(nodeTypeID).

    def write_node_csv(self, node_csv_file, node_type_id, data):

        with open(os.path.join(self.network_dir, node_csv_file), 'wt') as csv_file:

            # csvFile.write("# First column is node_type_id, remaining columns specify data. First row is header.\n")

            csv_writer = csv.writer(csv_file, delimiter=" ")

            # Write header
            row = ['node_type_id']
            for k in data:
                if type(k) == bytes:
                    row.append(k.decode())
                else:
                    row.append(k)

            csv_writer.writerow(row)

            # Write data to csv
            for (i, typeID) in enumerate(node_type_id):
                row = [typeID]

                for v in data.values():
                    if type(v[i]) == bytes:
                        row.append(v[i].decode())
                    else:
                        row.append(v[i])

                csv_writer.writerow(row)

    ############################################################################

    # !!! WHAT HAPPENS IF WE CANT KEEP EVERYTHING IN MEMORY

    def write_edges(self,
                    edge_file,
                    population_rows,
                    edge_type_id,
                    source_id,
                    target_id,
                    data,
                    include_index=False):

        # We need to have the targetID vector sorted
        # How should we handle cells that do not have any edges at all?
        sort_idx = np.argsort(target_id)

        # EVERYTHING we write needs to use sort_idx (idx) so it is in the right order
        with h5py.File(os.path.join(self.network_dir, edge_file), 'w', libver=self.h5py_libver) as f:
            self.add_version(f)

            edg_group = f.create_group("edges")

            for pop_name, pop_rows in population_rows.items():
                source_population_name, target_population_name = pop_name.split("_")
                n_rows = len(pop_rows)

                e_group = edg_group.create_group(pop_name)

                sort_idx = np.argsort(target_id[pop_rows])
                idx = pop_rows[sort_idx]

                e_group.create_dataset("edge_group_id", data=np.zeros((n_rows,), dtype=int))
                e_group.create_dataset("edge_group_index", data=np.arange(n_rows), dtype=int)
                e_group.create_dataset("edge_type_id", data=edge_type_id[idx])
                e_group.create_dataset("source_node_id", data=source_id[idx])
                e_group.create_dataset("target_node_id", data=target_id[idx])

                e_group["source_node_id"].attrs["node_population"] = source_population_name
                e_group["target_node_id"].attrs["node_population"] = target_population_name

                group_group = e_group.create_group("0")

                for data_type in data.keys():
                    if len(data[data_type].shape) == 1:
                        group_group.create_dataset(data_type, data=data[data_type][idx])
                    elif len(data[data_type].shape) == 2:
                        group_group.create_dataset(data_type, data=data[data_type][idx, :])
                    else:
                        print("Unsupported width of data column.")

                if include_index:
                    # We need to create the indices needed by Allen Institute
                    self.create_index(e_group["source_node_id"], e_group, index_source=True)
                    self.create_index(e_group["target_node_id"], e_group, index_source=False)

    ############################################################################

    # createIndex provided by Kael Dai from Allen Institute, 2018-11-27

    def create_index(self, node_ids_ds, output_grp, index_source=0):

        # TODO: Verify that this is correct! -- I am not sure it is...

        if not index_source:
            edge_nodes = np.array(node_ids_ds, dtype=np.int64)
            output_grp = output_grp.create_group('indices/target_to_source')
        else:
            edge_nodes = np.array(node_ids_ds, dtype=np.int64)
            output_grp = output_grp.create_group('indices/source_to_target')

        edge_nodes = np.append(edge_nodes, [-1])
        n_targets = np.max(edge_nodes)
        ranges_list = [[] for _ in range(n_targets + 1)]

        n_ranges = 0
        begin_index = 0
        cur_trg = edge_nodes[begin_index]
        for end_index, trg_gid in enumerate(edge_nodes):
            if cur_trg != trg_gid:
                ranges_list[cur_trg].append((begin_index, end_index))
                cur_trg = int(trg_gid)
                begin_index = end_index
                n_ranges += 1

        node_id_to_range = np.zeros((n_targets + 1, 2))
        range_to_edge_id = np.zeros((n_ranges, 2))
        range_index = 0
        for node_index, trg_ranges in enumerate(ranges_list):
            if len(trg_ranges) > 0:
                node_id_to_range[node_index, 0] = range_index
                for r in trg_ranges:
                    range_to_edge_id[range_index, :] = r
                    range_index += 1
                node_id_to_range[node_index, 1] = range_index
            else:
                node_id_to_range[node_index, :] = -1

        output_grp.create_dataset('range_to_edge_id', data=range_to_edge_id, dtype='uint64')
        output_grp.create_dataset('node_id_to_range', data=node_id_to_range, dtype='uint64')

    ############################################################################

    def write_edges_csv(self,
                        edge_csv_file,
                        edge_type_id,
                        data):

        with open(os.path.join(self.network_dir, edge_csv_file), 'wt') as csv_file:

            # csvFile.write("# First column is edge_type_id, remaining columns specify data. First row is header.\n")

            csv_writer = csv.writer(csv_file, delimiter=" ")

            # Write header
            row = ['edge_type_id']
            for k in data:
                if type(k) == bytes:
                    row.append(k.decode())
                else:
                    row.append(k)

            csv_writer.writerow(row)

            # Write data to csv
            for (i, typeID) in enumerate(edge_type_id):
                row = [typeID]

                for v in data.values():
                    if type(v[i]) == bytes:
                        row.append(v[i].decode())
                    else:
                        row.append(v[i])

                csv_writer.writerow(row)

    ############################################################################

    def write_input(self, spike_file_name, spike_times, gids):

        f_name = os.path.join(self.base_dir, spike_file_name)

        print(f"Writing spikes to {f_name}")

        with h5py.File(f_name, 'w', libver=self.h5py_libver) as f:
            self.add_version(f)

            print(f"Writing file {f_name}")

            s_group = f.create_group("spikes")
            s_group.attrs["sorting"] = "gid"
            s_group.create_dataset("gids", data=gids)
            s_group.create_dataset("timestamps", data=spike_times*1e3)  # Convert to ms

        return f_name

    ############################################################################

    # OBS, check that the files go in right directory

    def write_input_neurodamus(self, spike_file_name, spikes):

        # File standard is NWB, info about the neurodata without borders api:
        # http://neurodatawithoutborders.github.io/api-python/build/html/index.html

        # However, we use the simpler data file format:
        # with h5py.File(self.networkDir + fileNWB, 'w', libver=self.h5pyLibver) as f:
        #   print("Work in progress !!!")

        # with open(self.networkDir + spikeFileName,'wt') as f:
        with open(spike_file_name, 'wt') as f:
            f.write("/scatter")

            for time, gid in spikes:
                f.write(f"{time} {gid}")

    ############################################################################

    def add_version(self, hdf5_file):

        hdf5_file.attrs["version"] = [0, 1]
        hdf5_file.attrs["magic"] = 0x0A7A
