#!/usr/bin/env python3
import os

import numpy as np
import h5py

from snudda.utils.load import SnuddaLoad


class SnuddaLoadNetworkSimulation:

    def __init__(self, network_simulation_output_file=None, network_path=None):

        if network_simulation_output_file:
            self.network_simulation_output_file_name = network_simulation_output_file
        elif network_path:
            self.network_simulation_output_file_name = os.path.join(network_path, "simulation", "network-output.hdf5")
        else:
            self.network_simulation_output_file_name = None

        if network_path:
            self.network_path = network_path
        elif self.network_simulation_output_file_name:
            self.network_path = os.path.basename(os.path.basename(self.network_simulation_output_file_name))
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
            n_spikes += len(spikes)

        merged_spike_data = np.full((n_spikes, 2), np.nan)
        idx = 0

        for neuron_id, spikes in spike_data.items():
            merged_spike_data[idx:idx+len(spikes), 0] = spikes
            merged_spike_data[idx:idx+len(spikes), 1] = neuron_id
            idx += len(spikes)

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
            for nid in self.network_simulation_file["spikeData"]:
                spike_data[int(nid)] = self.network_simulation_file["spikeData"][nid][()]

        elif np.issubdtype(neuron_id, np.integer):
            if str(neuron_id) in self.network_simulation_file["spikeData"]:
                spike_data = self.network_simulation_file["spikeData"][str(neuron_id)][()]
            else:
                spike_data = np.array([])
        else:
            spike_data = dict()
            for nid in neuron_id:
                if str(nid) in self.network_simulation_file["spikeData"]:
                    spike_data[nid] = self.network_simulation_file["spikeData"][str(nid)][()].copy()
                else:
                    spike_data[nid] = np.array([])

        # If all neuronID not represented, add empty
        for nid in self.network_simulation_file["metaData/ID"]:
            if nid not in spike_data:
                spike_data[nid] = np.array([])

        return spike_data

    def get_voltage(self, neuron_id=None):
        """ Return volt data for neuron_id. """

        volt_data = dict()

        if neuron_id is None:
            for nid in self.network_simulation_file["voltData"]:
                if nid != "time":
                    volt_data[int(nid)] = self.network_simulation_file["voltData"][nid][()].copy()
        else:
            for nid in neuron_id:
                volt_data[int(nid)] = self.network_simulation_file["voltData"][str(nid)][()].copy()

        time_data = self.network_simulation_file["voltData"]["time"][()].copy()

        return volt_data, time_data

    def get_current(self, pre_id=None, post_id=None):
        assert "currentData" in self.network_simulation_file, "No synaptic current data in file"

        pre_id_all = self.network_simulation_file["currentData"]["preID"]
        post_id_all = self.network_simulation_file["currentData"]["postID"]

        if pre_id is not None and post_id is not None:
            idx = np.where(np.logical_and(pre_id_all == pre_id, post_id_all == post_id))[0]
        elif pre_id is not None:
            idx = np.where(pre_id_all == pre_id)
        elif post_id is not None:
            idx = np.where(post_id_all == post_id)
        else:
            idx = None

        if idx is not None:
            cur = self.network_simulation_file["currentData"]["current"][:, idx].copy()
        else:
            cur = self.network_simulation_file["currentData"]["current"].copy()

        time = self.network_simulation_file["currentData"]["time"].copy()

        return cur, time


    def get_neuron_positions(self, neuron_id=None):

        if neuron_id is None:
            pos_data = self.network_simulation_file["metaData/position"][()].copy()

        else:
            pos_data = self.network_simulation_file["metaData/position"][neuron_id, :].copy()

        return pos_data

    def get_id_of_neuron_type(self, neuron_type=None):

        if neuron_type:
            neuron_id = [x for x, y in zip(self.network_simulation_file["metaData/ID"],
                                           self.network_simulation_file["metaData/type"])
                         if y.lower() == neuron_type.lower()]
        else:
            neuron_id = self.network_simulation_file["metaData/ID"][()].copy()

        return neuron_id

    def get_neuron_name(self, neuron_id=None):
        if neuron_id:
            neuron_name = [SnuddaLoad.to_str(x) for x, y in zip(self.network_simulation_file["metaData/name"],
                                                                self.network_simulation_file["metaData/ID"])
                           if y in neuron_id]
        else:
            neuron_name = [SnuddaLoad.to_str(x) for x in self.network_simulation_file["metaData/name"]]

        return neuron_name

    def get_neuron_type(self, neuron_id=None):
        if neuron_id:
            neuron_type = [SnuddaLoad.to_str(x) for x, y in zip(self.network_simulation_file["metaData/type"],
                                                                self.network_simulation_file["metaData/ID"])
                           if y in neuron_id]
        else:
            neuron_type = [SnuddaLoad.to_str(x) for x in self.network_simulation_file["metaData/type"]]

        return neuron_type


def load_network_simulation_cli():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Load snudda activity data (spikes and or voltage)")
    parser.add_argument("dataFile", help="Data file")
    args = parser.parse_args()

    slna = SnuddaLoadNetworkSimulation(network_simulation_output_file=args.dataFile)
    slna.load()


if __name__ == "__main__":

    load_network_simulation_cli()