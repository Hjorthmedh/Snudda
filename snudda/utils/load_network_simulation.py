#!/usr/bin/env python3
import os

import h5py
import numpy as np
from collections import OrderedDict
from numba import jit

from snudda.utils.load import SnuddaLoad


class SnuddaLoadNetworkSimulation:

    def __init__(self, network_simulation_output_file=None,
                 network_path=None,
                 do_test=True,
                 verbose=False):

        self.verbose = verbose

        if network_simulation_output_file:
            self.network_simulation_output_file_name = network_simulation_output_file
        elif network_path:
            self.network_simulation_output_file_name = os.path.join(network_path, "simulation", "output.hdf5")
        else:
            self.network_simulation_output_file_name = None

        if network_path:
            self.network_path = network_path
        elif self.network_simulation_output_file_name:
            self.network_path = os.path.dirname(os.path.dirname(self.network_simulation_output_file_name))
        else:
            self.network_path = None

        self.network_simulation_file = None
        self.depolarisation_block = None
        self.test_data = do_test

        if self.network_simulation_output_file_name:
            self.load()

    def load(self, network_simulation_output_file=None):

        if not network_simulation_output_file:
            network_simulation_output_file = self.network_simulation_output_file_name

        print(f"Loading {network_simulation_output_file}")
        self.network_simulation_file = h5py.File(network_simulation_output_file, "r")

        if self.test_data:
            self.depolarisation_block = self.check_depolarisation_block()

    def close(self):
        if self.network_simulation_file:
            self.network_simulation_file.close()
            self.network_simulation_file = None

    def merge_spikes(self, spike_data=None):

        """ Merge spike_data dictionary into an array with spike time and neuron id as the columns,
            useful for plotting

            Args:
                spike_data : Dictionary with spike data (usually from get_spikes)
            """

        if spike_data is None:
            spike_data = self.get_spikes()

        n_spikes = 0
        for spikes in spike_data.values():
            n_spikes += spikes.size

        merged_spike_data = np.full((n_spikes, 2), np.nan)
        idx = 0

        for neuron_id, spikes in spike_data.items():
            if spikes.size > 0:
                merged_spike_data[idx:idx + spikes.size, 0] = spikes
                merged_spike_data[idx:idx + spikes.size, 1] = neuron_id
                idx += spikes.size

        # Sort the spikes after time
        idx = np.argsort(merged_spike_data[:, 0])
        merged_spike_data = merged_spike_data[idx, :]

        return merged_spike_data

    def get_spikes(self, neuron_id=None):

        """ Returns the spikes for neuron_id. If neuron_id is an integer, spike times are returned as an array.
            If neuron_id is a list or array, spike times are returned in a dictionary.

        Args:
            neuron_id : Neuron ID, either integer or list / array

        """

        if neuron_id is None:
            spike_data = OrderedDict()
            for nid in self.network_simulation_file["neurons"]:

                if "spikes" in self.network_simulation_file[f"neurons/{nid}"]:
                    spike_data[int(nid)] = self.network_simulation_file[f"neurons/{nid}/spikes"][()].copy()

            # If all neuronID not represented, add empty
            for nid in self.network_simulation_file["metaData/ID"]:
                if nid not in spike_data:
                    spike_data[nid] = np.array([])

        elif np.issubdtype(type(neuron_id), np.integer):
            if str(neuron_id) in self.network_simulation_file["neurons"] \
                    and "spikes" in self.network_simulation_file[f"neurons/{neuron_id}"]:
                spike_data = self.network_simulation_file[f"neurons/{neuron_id}/spikes"][()].copy()
            else:
                spike_data = np.array([])

        else:
            spike_data = OrderedDict()
            for nid in neuron_id:
                if str(nid) in self.network_simulation_file["neurons"] \
                        and "spikes" in self.network_simulation_file[f"neurons/{nid}"]:
                    spike_data[nid] = self.network_simulation_file[f"neurons/{nid}/spikes"][()].copy()
                else:
                    spike_data[nid] = np.array([])

        return spike_data

    def get_data(self, data_type, neuron_id=None):

        """ Returns data for neuron_id """

        data = OrderedDict()
        sec_id_x = OrderedDict()
        syn_info = OrderedDict()

        if neuron_id is None:
            neuron_id = self.network_simulation_file["neurons"].keys()

        for nid in neuron_id:
            snid = str(nid)
            inid = int(nid)

            try:
                if data_type in self.network_simulation_file["neurons"][snid].keys():
                    data[inid] = self.network_simulation_file["neurons"][snid][data_type][()].T.copy()
                    sec_id = self.network_simulation_file["neurons"][snid][data_type].attrs["sec_id"].copy()
                    sec_x = self.network_simulation_file["neurons"][snid][data_type].attrs["sec_x"].copy()

                    sec_id_x[inid] = (sec_id, sec_x)

                    if "synapse_type" in self.network_simulation_file["neurons"][snid][data_type].attrs:
                        synapse_type = self.network_simulation_file["neurons"][snid][data_type].attrs["synapse_type"].copy()
                        presynaptic_id = self.network_simulation_file["neurons"][snid][data_type].attrs["presynaptic_id"].copy()
                        cond = self.network_simulation_file["neurons"][snid][data_type].attrs["cond"].copy()
                        syn_info[inid] = (synapse_type, presynaptic_id, cond)
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

        return data, sec_id_x, syn_info

    def get_voltage(self, neuron_id=None):
        """ Return volt data for neuron_id. """

        orig_neuron_id = neuron_id

        if np.issubdtype(type(neuron_id), np.integer):
            neuron_id = [neuron_id]

        voltage, sec_id_x, _ = self.get_data("voltage", neuron_id=neuron_id)

        if np.issubdtype(type(orig_neuron_id), np.integer):
            return voltage[orig_neuron_id]
        else:
            return voltage

    @staticmethod
    @jit(nopython=True, fastmath=True, cache=True)
    def check_trace_depolarisation_block(neuron_id, time, voltage, threshold, max_duration):

        depolarisation_block = []

        ctr = 0
        t_start = 0
        dt = time[1] - time[0]
        n_limit = int(np.ceil(max_duration / dt))

        for idx in range(0, len(time)):
            if voltage[idx] > threshold:
                ctr += 1
            else:
                if ctr > n_limit:
                    depolarisation_block.append((neuron_id, t_start, time[idx]))

                t_start = time[idx]
                ctr = 0

        if ctr > n_limit:
            depolarisation_block.append((neuron_id, t_start, time[-1]))

        return depolarisation_block

    def check_depolarisation_block(self, threshold=-40e-3, max_duration=20e-3):

        if self.verbose:
            print("Checking neurons for depolarisation block")

        neuron_id_list = np.array(sorted([int(x) for x in self.network_simulation_file["neurons"].keys()]))

        if not ((np.diff(neuron_id_list) == 1).all() and neuron_id_list[0] == 0):
            print(f"Failed sanity check on neuron ID, not all neurons simulated? {neuron_id_list}")

        time = self.get_time()
        dt = time[1] - time[0]

        voltage, sec_id_x, _ = self.get_data("voltage", neuron_id=None)
        depolarisation_block = []

        n_limit = int(np.ceil(max_duration / dt))

        for neuron_id in neuron_id_list:

            depol_block = SnuddaLoadNetworkSimulation.check_trace_depolarisation_block(neuron_id=neuron_id,
                                                                                       time=time,
                                                                                       voltage=voltage[neuron_id],
                                                                                       threshold=threshold,
                                                                                       max_duration=max_duration)

            depolarisation_block = depolarisation_block + depol_block

        if self.verbose and len(depolarisation_block) > 0:
            for neuron_id, t_start, t_end in depolarisation_block:
                print(f"Neuron {neuron_id} has depolarisation block from {t_start:.3f} s to {t_end:.3f} s")

        bad_cells = sorted(list(set([x for x, ts, te in depolarisation_block])))
        bad_cell_str = [f"{x}: ({self.network_simulation_file['metaData/name'][x].decode()}, {self.network_simulation_file['metaData/parameterKey'][x].decode()}, {self.network_simulation_file['metaData/morphologyKey'][x].decode()})"
                        for x in bad_cells]
        bad_str = '\n'.join(bad_cell_str)

        if len(depolarisation_block) > 0:
            print(f"WARNING. Depolarisation block in neuron - neuron_id: (name, parameter_key, morphology_key):\n{bad_str}")

        return depolarisation_block

    def get_depolarisation_dictionary(self):

        depol_dict = dict()

        if self.depolarisation_block is None:
            self.depolarisation_block = self.check_depolarisation_block()

        for neuron_id, start_time, end_time in self.depolarisation_block:
            if neuron_id not in depol_dict:
                depol_dict[neuron_id] = []

            depol_dict[neuron_id].append((start_time, end_time))

        return depol_dict

    def get_time(self):

        t = self.network_simulation_file["time"][()].copy()
        return t

    def get_synaptic_current(self, pre_id=None, post_id=None):
        """ Return volt data for neuron_id. """

        current, sec_id_x, syn_info = self.get_data("synaptic_current", neuron_id=post_id)

        if pre_id is None:
            return current, sec_id_x, syn_info

        filtered_current = OrderedDict()
        filtered_sec_id_x = OrderedDict()
        filtered_syn_info = OrderedDict()

        for neuron_id, info in syn_info.items():
            idx = np.where(syn_info[1] == pre_id)[0]
            if len(idx) > 0:
                filtered_current[neuron_id] = current[neuron_id][:, idx]
                filtered_sec_id_x[neuron_id] = (sec_id_x[neuron_id][0][idx], sec_id_x[neuron_id][1][idx])
                filtered_syn_info[neuron_id] = (syn_info[neuron_id][0][idx], syn_info[neuron_id][1][idx],
                                                syn_info[neuron_id][2][idx])

        return filtered_current, filtered_sec_id_x, filtered_syn_info

    def get_neuron_positions(self, neuron_id=None):

        if neuron_id is None:
            pos_data = self.network_simulation_file["metaData/position"][()].copy()
        else:
            pos_data = self.network_simulation_file["metaData/position"][neuron_id, :].copy()

        return pos_data

    def get_id_of_neuron_type(self, neuron_type=None):

        if neuron_type:
            neuron_id = [x for x, y in zip(self.network_simulation_file["metaData/ID"][()],
                                           self.network_simulation_file["metaData/type"][()])
                         if SnuddaLoad.to_str(y).lower() == neuron_type.lower()]
        else:
            neuron_id = self.network_simulation_file["metaData/ID"][()].copy()

        return neuron_id

    def get_neuron_name(self, neuron_id=None):
        if neuron_id is not None:
            neuron_name = [SnuddaLoad.to_str(self.network_simulation_file["metaData/name"][x]) for x in neuron_id]
        else:
            neuron_name = [SnuddaLoad.to_str(x) for x in self.network_simulation_file["metaData/name"][()]]

        return neuron_name

    def get_neuron_type(self, neuron_id=None):
        if neuron_id is not None:
            neuron_type = [SnuddaLoad.to_str(self.network_simulation_file["metaData/type"][x]) for x in neuron_id]
        else:
            neuron_type = [SnuddaLoad.to_str(x) for x in self.network_simulation_file["metaData/type"][()]]

        return neuron_type

    def iter_neuron_type(self):
        neuron_types = list(set(self.get_neuron_type()))
        neuron_types.sort()
                
        for x in neuron_types:
            yield x

    def iter_neuron_id(self):
        for x in self.network_simulation_file["metaData/ID"]:
            yield x

    def get_neuron_keys(self, neuron_id):
        param_key = self.network_simulation_file["metaData/parameterKey"][neuron_id]
        morph_key = self.network_simulation_file["metaData/morphologyKey"][neuron_id]
        mod_key = self.network_simulation_file["metaData/modulationKey"][neuron_id]

        return param_key, morph_key, mod_key

    def export_to_txt(self, txt_file, neuron_id=None, time_scale=1.0):

        """
            Set time_scale to 1000 if you want ms

        """

        meta_file = f"{txt_file}-meta"

        print(f"Writing spikes to csv file {txt_file}")
        print(f"Writing metadata to {meta_file}")
        print("OBS, only neurons that have spikes are written to file.")

        with open(txt_file, "wt") as f, open(meta_file, "wt") as fm:

            spikes = self.get_spikes(neuron_id=neuron_id)

            for nid, spike_train in spikes.items():

                # Skip empty spike trains
                if spike_train.size > 0:
                    f.write(f"{' '.join([f'{time_scale*x:.5f}' for x in spike_train.flatten()])}\n")
                    fm.write(f"{nid}, {SnuddaLoad.to_str(self.network_simulation_file['metaData/name'][nid])}, "
                             f"{self.network_simulation_file['metaData/populationUnit'][nid]}, "
                             f"{','.join([str(x) for x in self.network_simulation_file['metaData/position'][nid,: ]])}\n")


def load_network_simulation_cli():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Load snudda activity data (spikes and or voltage)")
    parser.add_argument("dataFile", help="Data file")
    parser.add_argument("--export_spike_file", help="Name of csv file to export spikes to",
                        default=None)
    parser.add_argument("--time_scale", default=1.0, type=float)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--skip_test", help="Do tests on simulation data", action="store_true")
    args = parser.parse_args()

    slna = SnuddaLoadNetworkSimulation(network_simulation_output_file=args.dataFile, verbose=args.verbose,
                                       do_test=not args.skip_test)

    if args.export_spike_file is not None:
        slna.export_to_txt(txt_file=args.export_spike_file, time_scale=args.time_scale)


if __name__ == "__main__":
    load_network_simulation_cli()
