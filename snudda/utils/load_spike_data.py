#!/usr/bin/env python3

import numpy as np
from snudda.utils.load import SnuddaLoad


class SnuddaLoadSpikeData:

    def __init__(self, spike_data_file, network_file):

        assert False, "This is now obsolete. Please use load_network_simulation.py"

        self.spike_data_file = spike_data_file
        self.network_file = network_file

        self.snudda_load = SnuddaLoad(network_file=self.network_file)

        spike_data = np.genfromtxt(self.spike_data_file, delimiter='\t')
        self.time = spike_data[:, 0] * 1e-3
        self.spike_id = spike_data[:, 1].astype(int)

    def get_freq(self,
                 neuron_id=None,
                 neuron_name=None,
                 neuron_type=None,
                 time_range=None):

        spikes, nid = self.get_spikes_from(neuron_id=neuron_id,
                                           neuron_name=neuron_name,
                                           neuron_type=neuron_type,
                                           time_range=time_range)

        if time_range and time_range[0]:
            min_t = time_range[0]
        else:
            min_t = 0

        if time_range and time_range[1]:
            max_t = time_range[1]
        else:
            assert np.sum([len(x) for x in spikes]) > 0, "No spikes found for any neuron in simulation output file."
            max_t = np.max([np.max(x) for x in spikes if len(x) > 0])
            print(f"Assuming max time {max_t}")

        print(f"Time range {min_t} to {max_t}")

        freq_list = [len(x) / (max_t - min_t) for x in spikes]

        return freq_list, nid

    def get_spikes_from(self,
                        neuron_id=None,
                        neuron_name=None,
                        neuron_type=None,
                        time_range=None
                        ):

        if not neuron_id:
            if neuron_name:
                neuron_id = self.snudda_load.get_neuron_id_with_name(neuron_name=neuron_name)
            elif neuron_type:
                neuron_id = self.snudda_load.get_neuron_id_of_type(neuron_type=neuron_type)
            else:
                neuron_id = [x["neuronID"] for x in self.snudda_load.data["neurons"]]

        if type(neuron_id) == list:
            spikes = [self.get_spikes(neuron_id=nid, time_range=time_range) for nid in neuron_id]
        else:
            print("No neurons specified.")
            spikes = []

        return spikes, neuron_id

    def print_freq_info(self,
                      neuron_id=None,
                      neuron_name=None,
                      neuron_type=None,
                      time_range=None
                      ):

        freq_list, nid = self.get_freq(neuron_id=neuron_id,
                                       neuron_name=neuron_name,
                                       neuron_type=neuron_type,
                                       time_range=time_range)

        neuron_name = [self.snudda_load.data["neurons"][n]["name"] for n in nid]

        sort_idx = np.argsort(neuron_name)

        for si in sort_idx:
            freq = freq_list[si]
            nrn_id = nid[si]
            name = neuron_name[si]

            print(f"{name} ({nrn_id}) {freq} Hz in time range {time_range}")

    def get_spikes(self, neuron_id, time_range=None):

        idx = np.where(self.spike_id == neuron_id)[0]
        spikes = np.sort(self.time[idx])

        if time_range:
            t_min = time_range[0]
            t_max = time_range[1]

            if t_min and t_max:
                t_idx = np.where((t_min < spikes) & (spikes < t_max))[0]
            elif t_min:
                t_idx = np.where(t_min < spikes)[0]
            elif t_max:
                t_idx = np.where(spikes < t_max)[0]
            else:
                assert False, f"get_spikes: time_range={time_range}, expected (t_min,t_max) where both not None."

            spikes = spikes[t_idx]

        return spikes


def snudda_load_spike_data_cli():

    """ Command line parser for SnuddaLoadSpikeData script """

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Load snudda spike data, display firing frequency")
    parser.add_argument("networkFile", help="Network file")
    parser.add_argument("spikeFile", help="Spike data file")
    parser.add_argument("--neuronName", help="Name of neurons", type=str, default=None)
    parser.add_argument("--neuronType", help="Type of neurons", type=str, default=None)
    parser.add_argument("--neuronID", help="Neuron ID of neuron", type=int, default=None)
    parser.add_argument("--minT", help="Minimum time", type=float, default=None)
    parser.add_argument("--maxT", help="Maximum time", type=float, default=None)

    args = parser.parse_args()

    if args.minT or args.maxT:
        time_range = (args.minT, args.maxT)
    else:
        time_range = None

    sd = SnuddaLoadSpikeData(network_file=args.networkFile,
                             spike_data_file=args.spikeFile)

    sd.print_freq_info(neuron_id=args.neuronID,
                       neuron_name=args.neuronName,
                       neuron_type=args.neuronType,
                       time_range=time_range)


if __name__ == "__main__":

    snudda_load_spike_data_cli()