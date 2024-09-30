#!/usr/bin/env python3
import os

import h5py
import numpy as np
from numba import jit

from snudda.utils.load import SnuddaLoad


class SnuddaLoadSimulation:

    def __init__(self, network_simulation_output_file=None,
                 network_path=None,
                 do_test=True,
                 verbose=False,
                 quiet_load=False):

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
            self.load(quiet=quiet_load)

    def load(self, network_simulation_output_file=None, quiet=False):

        if not network_simulation_output_file:
            network_simulation_output_file = self.network_simulation_output_file_name

        print(f"Loading {network_simulation_output_file}")
        self.network_simulation_file = h5py.File(network_simulation_output_file, "r")

        if self.test_data:
            self.depolarisation_block = self.check_depolarisation_block(quiet=quiet)

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
            spike_data = dict()
            for nid in self.network_simulation_file["neurons"]:

                if "spikes" in self.network_simulation_file[f"neurons/{nid}"]:
                    spike_data[int(nid)] = self.network_simulation_file[f"neurons/{nid}/spikes"][()].copy()

            # If all neuronID not represented, add empty
            for nid in self.network_simulation_file["meta_data/id"]:
                if nid not in spike_data:
                    spike_data[nid] = np.array([])

        elif np.issubdtype(type(neuron_id), np.integer):
            if str(neuron_id) in self.network_simulation_file["neurons"] \
                    and "spikes" in self.network_simulation_file[f"neurons/{neuron_id}"]:
                spike_data = self.network_simulation_file[f"neurons/{neuron_id}/spikes"][()].copy()
            else:
                spike_data = np.array([])

        elif isinstance(neuron_id, np.ndarray) and neuron_id.size == 0:
            spike_data = dict()
        else:
            spike_data = dict()
            try:
                for nid in neuron_id:
                    if str(nid) in self.network_simulation_file["neurons"] \
                            and "spikes" in self.network_simulation_file[f"neurons/{nid}"]:
                        spike_data[nid] = self.network_simulation_file[f"neurons/{nid}/spikes"][()].copy()
                    else:
                        spike_data[nid] = np.array([])
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

        return spike_data

    def get_frequency(self, neuron_id, time_ranges=None):

        # neuron_id must be a list or an array, time_ranges is a list of tuples

        spike_data = self.get_spikes(neuron_id=neuron_id)

        if time_ranges is None:
            time_ranges = [0, np.max(self.get_time())]

        freq_table = np.zeros((len(spike_data), len(time_ranges)))

        for idx, (n_id, spikes) in enumerate(spike_data.items()):

            assert neuron_id is None or neuron_id[idx] == n_id, f"Order of neuron_id is not maintained"

            for t_idx, tr in enumerate(time_ranges):
                spike_count = np.sum(np.logical_and(tr[0] <= spikes, spikes <= tr[1]))
                freq_table[idx, t_idx] = spike_count / (tr[1] - tr[0])

        return freq_table

    def list_data_types(self, neuron_id):
        return list(self.network_simulation_file["neurons"][str(neuron_id)].keys())

    def get_all_data(self, neuron_id, exclude=None, include_time=False):

        data = dict()

        for data_type in self.list_data_types(neuron_id=neuron_id):

            if exclude is not None and data_type in exclude:
                continue

            if data_type not in data:
               data[data_type] = self.get_data(data_type=data_type, neuron_id=neuron_id)
            else:
                data[data_type].update(self.get_data(data_type=data_type, neuron_id=neuron_id))

        if include_time:
            data["time"] = self.get_time()

        return data

    def get_data(self, data_type, neuron_id=None):

        """ Returns data for neuron_id """

        data = dict()
        sec_id_x = dict()
        syn_info = dict()

        if neuron_id is None:
            neuron_id = self.network_simulation_file["neurons"].keys()
        elif np.issubdtype(type(neuron_id), np.integer):
            neuron_id = [neuron_id]

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

    def check_depolarisation_block(self, threshold=-40e-3, max_duration=200e-3, quiet=False):

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

            # We should only check for depolarisation block in the soma, sec_id = -1
            v_idx = np.where(sec_id_x[neuron_id][0] == -1)[0][0]
            depol_block = SnuddaLoadSimulation.check_trace_depolarisation_block(neuron_id=neuron_id,
                                                                                time=time,
                                                                                voltage=voltage[neuron_id][:, v_idx],
                                                                                threshold=threshold,
                                                                                max_duration=max_duration)

            depolarisation_block = depolarisation_block + depol_block

        if self.verbose and len(depolarisation_block) > 0:
            for neuron_id, t_start, t_end in depolarisation_block:
                print(f"Neuron {neuron_id} has depolarisation block from {t_start:.3f} s to {t_end:.3f} s")

        bad_cells = sorted(list(set([x for x, ts, te in depolarisation_block])))
        bad_cell_str = [f"{x}: ({self.network_simulation_file['meta_data/name'][x].decode()}, {self.network_simulation_file['meta_data/parameter_key'][x].decode()}, {self.network_simulation_file['meta_data/morphology_key'][x].decode()})"
                        for x in bad_cells]
        bad_str = '\n'.join(bad_cell_str)

        if len(depolarisation_block) > 0 and not quiet:
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

        filtered_current = dict()
        filtered_sec_id_x = dict()
        filtered_syn_info = dict()

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
            pos_data = self.network_simulation_file["meta_data/position"][()].copy()
        else:
            pos_data = self.network_simulation_file["meta_data/position"][neuron_id, :].copy()

        return pos_data

    def get_id_of_neuron_type(self, neuron_type=None):

        if neuron_type:
            neuron_id = [x for x, y in zip(self.network_simulation_file["meta_data/id"][()],
                                           self.network_simulation_file["meta_data/type"][()])
                         if SnuddaLoad.to_str(y).lower() == neuron_type.lower()]
        else:
            neuron_id = self.network_simulation_file["meta_data/id"][()].copy()

        return neuron_id

    def get_neuron_name(self, neuron_id=None):
        if neuron_id is not None:
            neuron_name = [SnuddaLoad.to_str(self.network_simulation_file["meta_data/name"][x]) for x in neuron_id]
        else:
            neuron_name = [SnuddaLoad.to_str(x) for x in self.network_simulation_file["meta_data/name"][()]]

        return neuron_name

    def get_neuron_type(self, neuron_id=None):
        if neuron_id is not None:
            neuron_type = [SnuddaLoad.to_str(self.network_simulation_file["meta_data/type"][x]) for x in neuron_id]
        else:
            neuron_type = [SnuddaLoad.to_str(x) for x in self.network_simulation_file["meta_data/type"][()]]

        return neuron_type

    def iter_neuron_type(self):
        neuron_types = list(set(self.get_neuron_type()))
        neuron_types.sort()
                
        for x in neuron_types:
            yield x

    def iter_neuron_id(self):
        for x in self.network_simulation_file["meta_data/id"]:
            yield x

    def get_neuron_keys(self, neuron_id):
        param_key = self.network_simulation_file["meta_data/parameter_key"][neuron_id]
        morph_key = self.network_simulation_file["meta_data/morphology_key"][neuron_id]
        mod_key = self.network_simulation_file["meta_data/modulation_key"][neuron_id]

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
                    fm.write(f"{nid}, {SnuddaLoad.to_str(self.network_simulation_file['meta_data/name'][nid])}, "
                             f"{self.network_simulation_file['meta_data/population_unit'][nid]}, "
                             f"{','.join([str(x) for x in self.network_simulation_file['meta_data/position'][nid,: ]])}\n")


class SnuddaLoadNetworkSimulation (SnuddaLoadSimulation):

    def __init__(self, *args, **kwargs):
        raise DeprecationWarning("Please use SnuddaLoadSimulation instead of SnuddaLoadNetworkSimulation")
        super(*args, **kwargs)


def load_network_simulation_cli():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Load snudda activity data (spikes and or voltage)")
    parser.add_argument("data_file", help="Data file")
    parser.add_argument("--export_spike_file", help="Name of csv file to export spikes to",
                        default=None)
    parser.add_argument("--time_scale", default=1.0, type=float)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--skip_test", help="Do tests on simulation data", action="store_true")
    args = parser.parse_args()

    slna = SnuddaLoadSimulation(network_simulation_output_file=args.data_file, verbose=args.verbose,
                                do_test=not args.skip_test)

    if args.export_spike_file is not None:
        slna.export_to_txt(txt_file=args.export_spike_file, time_scale=args.time_scale)


if __name__ == "__main__":
    load_network_simulation_cli()
