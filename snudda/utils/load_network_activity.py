import numpy as np
import h5py


class SnuddaLoadNetworkActivity:

    def __init__(self, network_activity_file):

        self.network_activity_file_name = network_activity_file
        self.activity_file = None

    def load(self, network_activitiy_file=None):

        if not network_activitiy_file:
            network_activitiy_file = self.network_activity_file_name

        self.activity_file = h5py.File(network_activitiy_file, "r")

    def close(self):
        if self.activity_file:
            self.activity_file.close()
            self.activity_file = None

    def merge_spikes(self, spike_data):

        """ Merge spike_data dictionary into an array with spike time and neuron id as the columns,
            useful for plotting

            Args:
                spike_data : Dictionary with spike data (usually from get_spikes)
            """

        n_spikes = 0
        for spikes in spike_data.values():
            n_spikes += len(spikes)

        merged_spike_data = np.full((n_spikes, 2), np.nan)
        idx = 0

        for neuron_id, spikes in spike_data.items():
            merged_spike_data[idx:idx+len(spikes), 0] = spikes
            merged_spike_data[idx:idx+len(spikes), 1] = neuron_id

        return merged_spike_data

    def get_spikes(self, neuron_id=None):

        """ Returns the spikes for neuron_id. If neuron_id is an integer, spike times are returned as an array.
            If neuron_id is a list or array, spike times are returned in a dictionary.

        Args:
            neuron_id : Neuron ID, either integer or list / array

        """

        if neuron_id is None:
            spike_data = dict()
            for nid in self.activity_file["spikeData"]:
                spike_data[int(nid)] = self.activity_file["spikeData"][nid]

        if np.issubdtype(neuron_id, np.integer):
            if str(neuron_id) in self.activity_file["spikeData"]:
                spike_data = self.activity_file["spikeData"][str(neuron_id)]
            else:
                spike_data = np.array([])
        else:
            spike_data = dict()
            for nid in neuron_id:
                if str(nid) in self.activity_file["spikeData"]:
                    spike_data[nid] = self.activity_file["spikeData"][str(nid)]
                else:
                    spike_data[nid] = np.array([])

        return spike_data

    def get_voltage(self, neuron_id=None):
        """ Return volt data for neuron_id. """

        if neuron_id is None:
            volt_data = self.activity_file["voltData"].copy()
        else:
            volt_data = self.activity_file["voltData"][:, neuron_id].copy()

        return volt_data

    def get_neuron_positions(self, neuron_id=None):

        if neuron_id is None:
            pos_data = self.activity_file["position"].copy()

        else:
            pos_data = self.activity_file["position"][neuron_id, :].copy()

        return pos_data

    def get_id_of_neuron_type(self, neuron_type):
        neuron_id = [x for x, y in zip(self.activity_file["neuronID"], self.activity_file["type"])
                     if y.lower() == neuron_type.lower()]

        return neuron_id

    def cli(self):
        from argparse import ArgumentParser

        parser = ArgumentParser(description="Load snudda activity data (spikes and or voltage)")
        parser.add_argument("dataFile", help="Data file", dest="data_file")
        args = parser.parse_args()

        SnuddaLoadNetworkActivity(network_activity_file=args.data_file)
