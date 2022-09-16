#!/usr/bin/env python3
import os

import h5py
import numpy as np
from collections import OrderedDict

from snudda.utils.load import SnuddaLoad


class SnuddaLoadNetworkSimulation:

    def __init__(self, network_simulation_output_file=None, network_path=None):

        if network_simulation_output_file:
            self.network_simulation_output_file_name = network_simulation_output_file
        elif network_path:
            self.network_simulation_output_file_name = os.path.join(network_path,
                                                                    "simulation",
                                                                    "output.hdf5")
        else:
            self.network_simulation_output_file_name = None

        if network_path:
            self.network_path = network_path
        elif self.network_simulation_output_file_name:
            self.network_path = os.path.dirname(os.path.dirname(self.network_simulation_output_file_name))
        else:
            self.network_path = None

        self.network_simulation_file = None

        if self.network_simulation_output_file_name:
            self.load()

    def load(self, network_simulation_output_file=None):

        if not network_simulation_output_file:
            network_simulation_output_file = self.network_simulation_output_file_name

        print(f"Loading {network_simulation_output_file}")
        self.network_simulation_file = h5py.File(network_simulation_output_file, "r")

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
                    spike_data[int(nid)] = self.network_simulation_file[f"neurons/{nid}/spikes/data"][()].copy()

            # If all neuronID not represented, add empty
            for nid in self.network_simulation_file["metaData/ID"]:
                if nid not in spike_data:
                    spike_data[nid] = np.array([])

        elif np.issubdtype(type(neuron_id), np.integer):
            if str(neuron_id) in self.network_simulation_file["neurons"] \
                    and "spikes" in self.network_simulation_file[f"neurons/{neuron_id}"]:
                spike_data = self.network_simulation_file[f"neurons/{neuron_id}/spikes/data"][()].copy()
            else:
                spike_data = np.array([])

        else:
            spike_data = OrderedDict()
            for nid in neuron_id:
                if str(nid) in self.network_simulation_file["neurons"] \
                        and "spikes" in self.network_simulation_file[f"neurons/{nid}"]:
                    spike_data[nid] = self.network_simulation_file[f"neurons/{nid}/spikes/data"][()].copy()
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
            if data_type in self.network_simulation_file["neurons"][snid]:
                data[inid] = self.network_simulation_file["neurons"][snid][data_type]["data"][()].T.copy()
                sec_id = self.network_simulation_file["neurons"][snid][data_type]["sec_id"][()].copy()
                sec_x = self.network_simulation_file["neurons"][snid][data_type]["sec_x"][()].copy()

                sec_id_x[inid] = (sec_id, sec_x)

                if "synapse_type" in self.network_simulation_file["neurons"][snid][data_type]:
                    synapse_type = self.network_simulation_file["neurons"][snid][data_type]["synapse_type"][()].copy()
                    presynaptic_id = self.network_simulation_file["neurons"][snid][data_type]["presynaptic_id"][()].copy()
                    cond = self.network_simulation_file["neurons"][snid][data_type]["cond"][()].copy()
                    syn_info[inid] = (synapse_type, presynaptic_id, cond)

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

    """
    def check_depolarisation_block(self, threshold=0, max_duration=100e-3):

        neuron_id = np.array(sorted(self.network_simulation_file["neurons"].keys()))
        assert (np.diff(neuron_id) == 1).all() and neuron_id[0] == 0, f"Failed sanity check on neuron ID"

        block_flag = np.zeros((len(neuron_id),), dtype=bool)

        time = self.get_time()
        dt = time[1] - time[0]
        csum_threshold = max_duration / dt

        voltage, sec_id_x, _ = self.get_data("voltage", neuron_id=None)

        for neuron_id, (volt, sidx) in enumerate(zip(voltage, sec_id_x)):
            v_thresh = volt > threshold

            v_idx = np.where(v_thresh)[0]
            np.diff(v_idx) 


        for idx, (volt, sidx) in enumerate(zip(voltage, sec_id_x)):
            v_idx = np.where(sidx[0] == 0)
            csum = 0
            max_csum = 0
            block_times = []
            t_start_block = 0

            for t_idx, vt in enumerate(voltage[:, v_idx]):

                if vt > threshold:
                    csum += 1
                    if csum > max_csum:
                        max_csum = csum
                else:
                    csum = 0
                    t_start_block = t_idx

            if max_csum > csum_threshold:
                block_flag[idx] = True

        blocked_neurons = np.where(block_flag)
        return blocked_neurons

    """

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
            neuron_name = [SnuddaLoad.to_str(x) for x, y in zip(self.network_simulation_file["metaData/name"][()],
                                                                self.network_simulation_file["metaData/ID"][()])
                           if y in neuron_id]
        else:
            neuron_name = [SnuddaLoad.to_str(x) for x in self.network_simulation_file["metaData/name"][()]]

        return neuron_name

    def get_neuron_type(self, neuron_id=None):
        if neuron_id is not None:
            neuron_type = [SnuddaLoad.to_str(x) for x, y in zip(self.network_simulation_file["metaData/type"][()],
                                                                self.network_simulation_file["metaData/ID"][()])
                           if y in neuron_id]
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

        network_loader = SnuddaLoad(network_file=self.network_path)
        network_data = network_loader.data

        print(f"Writing spikes to csv file {txt_file}")
        print(f"Writing metadata to {meta_file}")
        print("OBS, only neurons that have spikes are written to file.")

        with open(txt_file, "wt") as f, open(meta_file, "wt") as fm:

            spikes = self.get_spikes(neuron_id=neuron_id)

            for nid, spike_train in spikes.items():

                # Skip empty spike trains
                if spike_train.size > 0:
                    f.write(f"{' '.join([f'{time_scale*x:.5f}' for x in spike_train.flatten()])}\n")
                    fm.write(f"{nid}, {network_data['neurons'][nid]['name']}, "
                             f"{network_data['populationUnit'][nid]}, "
                             f"{','.join([str(x) for x in network_data['neurons'][nid]['position']])}\n")


def load_network_simulation_cli():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Load snudda activity data (spikes and or voltage)")
    parser.add_argument("dataFile", help="Data file")
    parser.add_argument("--export_spike_file", help="Name of csv file to export spikes to",
                        default=None)
    parser.add_argument("--time_scale", default=1.0, type=float)
    args = parser.parse_args()

    slna = SnuddaLoadNetworkSimulation(network_simulation_output_file=args.dataFile)
    slna.load()

    if args.export_spike_file is not None:
        slna.export_to_txt(txt_file=args.export_spike_file, time_scale=args.time_scale)


if __name__ == "__main__":
    load_network_simulation_cli()
