# Must have a __init__.py file in module directory

import json
import h5py
import csv
from collections import OrderedDict
import os
import numpy as np

# TODO: This needs to be updated to match the latest SONATA release.

class ConvHurt(object):

    def __init__(self, simulation_structure, input_structures, base_dir="TEST/"):

        self.base_dir = base_dir
        self.network_dir = os.path.join(base_dir, 'networks')

        self.h5py_libver = "earliest"  # "latest"

        self.setup_directories(simulation_structure=simulation_structure,
                               input_structures=input_structures, base_dir=base_dir)

        self.write_main_config(simulation_structure=simulation_structure,
                               input_structures=input_structures,
                               base_dir=base_dir)

    ############################################################################

    @staticmethod
    def create_dir(directory):

        if not os.path.exists(directory):
            print(f"Creating directory : {directory}")
            os.makedirs(directory)

    ############################################################################

    def setup_directories(self,
                          simulation_structure,
                          input_structures,
                          base_dir=None):

        if base_dir is None:
            base_dir = "."

        self.create_dir(base_dir)

        dir_list = ['components', 'components/biophysical_neuron_dynamics', 'components/hoc_templates',
                    'components/mechanisms', 'components/morphologies', 'components/point_neuron_dynamics',
                    'components/synapse_dynamics', 'networks', 'inputs', 'output', f"networks/{simulation_structure}"]

        for x in input_structures:
            dir_list.append(f"networks/{x}")

        for x in dir_list:
            self.create_dir(f"{base_dir}{x}")

    ############################################################################
    #
    # Writing the main configuration file, linking everything together

    def write_main_config(self,
                          base_dir="TEST/",
                          out_file="circuit_config.json",
                          simulation_structure="striatum",
                          input_structures=None,
                          target_simulator="NEURON"):

        if input_structures is None:
            input_structures = ["thalamus", "cortex"]

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

        sim_node = os.path.join("$NETWORK_DIR", simulation_structure, f"{simulation_structure}_nodes.hdf5")
        sim_node_type = os.path.join("$NETWORK_DIR", simulation_structure, f"{simulation_structure}_node_types.csv")

        sim_edge = os.path.join("$NETWORK_DIR", simulation_structure, f"{simulation_structure}_edges.hdf5")
        sim_edge_type = os.path.join("$NETWORK_DIR", simulation_structure, f"{simulation_structure}_edge_types.csv")

        # We create a list where first element contains the simulated node,
        # subsequent elements are the inputs (can be more than one)

        node_files = [OrderedDict([("nodes_file", sim_node),
                                  ("node_types_file", sim_node_type)])]
        node_files += [OrderedDict([("nodes_file", os.path.join("$NETWORK_DIR", x, f"{x}_nodes.hdf5")),
                                   ("node_types_file", os.path.join("$NETWORK_DIR", x, f"{x}_node_types.csv"))])
                       for x in input_structures]

        edge_files = [OrderedDict([("edges_file", sim_edge), ("edge_types_file", sim_edge_type)])]
        edge_files += [OrderedDict([("edges_file", os.path.join("$NETWORK_DIR", x, f"{x}_edges.hdf5")),
                                    ("edge_types_file", os.path.join("$NETWORK_DIR", x, f"{x}_edge_types.csv"))])
                       for x in input_structures]

        config["networks"] = OrderedDict([("nodes", node_files),
                                          ("edges", edge_files)])

        with open(os.path.join(base_dir, out_file), 'wt') as f:
            json.dump(config, f, indent=4)

    ############################################################################

    # nodeID - vector with IDs for all neurons
    # nodeTypeID - vector with Type ID for each neuron
    # nodeGroupID - vector with group memebership for each neuron
    # nodeGroupIndex - vector with index of each neuron within the given group

    # data = dictionary with "position", "rotaton_angle_zaxis", ...
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
                    model_template=None):

        with h5py.File(os.path.join(self.network_dir, node_file), 'w', libver=self.h5py_libver) as f:

            n_group = f.create_group("nodes")

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

    def write_edges(self, edge_file,
                    edge_population_name,
                    edge_group,
                    edge_group_index,
                    edge_type_id,
                    source_id,
                    target_id,
                    data,
                    source_population_name=None,
                    target_population_name=None, ):

        if source_population_name is None:
            source_population_name = edge_population_name
            print(f"No source population name given, using edge_population_name: {edge_population_name}")

        if target_population_name is None:
            target_population_name = edge_population_name
            print(f"No source population name given, using edgePopulationName: {edge_population_name}")

            # We need to have the targetID vector sorted
        # How should we handle cells that do not have any edges at all?
        sort_idx = np.argsort(target_id)

        # This finds the positions where the sorter targetID increases
        all_diffs = np.diff(target_id[sort_idx])
        diff_idx = np.where(all_diffs > 0)[0]
        index_pointer = np.zeros(len(diff_idx) + 1)
        index_pointer[1:] = diff_idx + 1

        # import pdb
        # pdb.set_trace()

        # We assume all neurons has at least one synapse, the diff should
        # never jump by two
        if not (all_diffs < 2).all():
            print("!!! Not all neurons have synapses!")

        # EVERYTHING we write needs to use sortIdx so it is in the right order
        with h5py.File(os.path.join(self.network_dir, edge_file), 'w', libver=self.h5py_libver) as f:
            edg_group = f.create_group("edges")
            e_group = edg_group.create_group(edge_population_name)

            e_group.create_dataset("edge_group", data=edge_group[sort_idx])
            e_group.create_dataset("edge_group_index", data=edge_group_index[sort_idx])
            e_group.create_dataset("edge_type_id", data=edge_type_id[sort_idx])
            e_group.create_dataset("index_pointer", data=index_pointer)
            e_group.create_dataset("source_node_id", data=source_id[sort_idx])
            e_group.create_dataset("target_node_id", data=target_id[sort_idx])

            e_group["source_node_id"].attrs["node_population"] = source_population_name
            e_group["target_node_id"].attrs["node_population"] = target_population_name
            groups = np.unique(edge_group)

            for g in groups:
                idx = np.where(edge_group[sort_idx] == g)

                group_group = e_group.create_group(str(g))

                for data_type in data.keys():
                    if len(data[data_type].shape) == 1:
                        group_group.create_dataset(data_type,
                                                   data=data[data_type][sort_idx][idx])
                    elif len(data[data_type].shape) == 2:
                        group_group.create_dataset(data_type,
                                                   data=data[data_type][sort_idx, :][idx, :])
                    else:
                        print("Unsupported width of data column.")

            # We need to create the indices needed by Allen Institute
            self.create_index(e_group["source_node_id"], e_group, index_source=True)
            self.create_index(e_group["target_node_id"], e_group, index_source=False)
            # inGroup = eGroup.create_group("indices")
            # inGroup.create_group("source_to_target")
            # inGroup.create_group("target_to_source")

    ############################################################################

    # createIndex provided by Kael Dai from Allen Institute, 2018-11-27

    def create_index(self, node_ids_ds, output_grp, index_source=0):

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

        output_grp.create_dataset('range_to_edge_id', data=range_to_edge_id, dtype='uint64')
        output_grp.create_dataset('node_id_to_ranges', data=node_id_to_range, dtype='uint64')

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

    def write_input(self, spike_file_name, spikes):

        if spikes is None:
            print(f"No spikes specified, not writing {spike_file_name}")
            print("Use python3 Network_input.py yourinput.json yournetwork.hdf5 input-spikes.hdf5")
            return

        f_name = os.path.join(self.base_dir, spike_file_name)

        with h5py.File(f_name, 'w', libver=self.h5py_libver) as f:
            print(f"Writing file {f_name}")

            s_group = f.create_group("spikes")
            s_group.create_dataset("gids", data=spikes[:, 1])
            s_group.create_dataset("timestamps", data=spikes[:, 0])

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


if __name__ == "__main__":
    # ch = ConvHurt()
    # ch = ConvHurt(simulationStructure="cerebellum",
    #              inputStructures=["pons","cortex"])

    ch = ConvHurt(simulation_structure="striatum",
                  input_structures=["cortex", "thalamus"])

    # Test example, we have 5 neurons, big network
    # two groups

    node_data = {"positions": np.array([[1., 2., 3.], [4., 5., 6.], [7., 8., 9.],
                                        [1., 8., 9.], [2., 3., 2.]]),
                 "rotation_angle_zaxis": np.array([0.1, 0.2, 0.3, 0.4, 0.5])}

    ch.write_nodes(node_file='striatum_nodes.hdf5',
                   data=node_data,
                   node_id=np.array([0, 1, 2, 3, 4]),
                   population_name="striatum_nodes",
                   node_type_id=np.array([0, 1, 0, 1, 0]),
                   node_group_id=np.array([0, 0, 1, 1, 1]),
                   node_group_index=np.array([0, 1, 0, 1, 2]))

    node_type_id = np.array([0, 1])
    node_data_csv = OrderedDict([('name', ['A', 'B']),
                                 ('location', ['structA', 'structB'])])

    ch.write_node_csv(node_csv_file='striatum_node_types.csv',
                      node_type_id=node_type_id,
                      data=node_data_csv)

    edge_group = np.array([5, 5, 11, 11, 11])
    edge_group_index = np.array([0, 1, 0, 1, 2])
    edge_type_id = np.array([0, 1, 0, 1, 0])
    source_gid = np.array([1, 2, 3, 3, 4])
    target_gid = np.array([2, 3, 4, 0, 1])  # THESE ARE SORTED ... HAHAHA

    # Delay should be in ms (bad bad people, real scientists use SI units)
    edge_data = OrderedDict([("sec_id", np.array([10, 22, 33, 24, 15])),
                             ("sec_x", np.array([0.1, 0.3, 0.5, 0.2, 0])),
                             ("syn_weight", np.array([0.1e-9, 2e-9, 3e-9,
                                                     0.3e-9, 0.1e-9])),
                             ("delay", 1e3 * np.array([1e-3, 4e-3, 2e-3, 5e-3, 1e-3]))])

    ch.write_edges(edge_file="striatum_edges.hdf5",
                   edge_group=edge_group,
                   edge_group_index=edge_group_index,
                   edge_type_id=edge_type_id,
                   edge_population_name="striatum_edges",
                   source_id=source_gid,
                   target_id=target_gid,
                   data=edge_data)

    edge_type_id = np.array([0, 1])
    edge_csv_data = OrderedDict([('template', ['Exp2Syn', 'NULL']),
                                 ('dynamics_params', ["mysyn.json", 'yoursyn.json'])])

    ch.write_edges_csv(edge_csv_file="striatum_edge_types.csv",
                       edge_type_id=edge_type_id,
                       data=edge_csv_data)
