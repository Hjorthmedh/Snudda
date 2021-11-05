#!/usr/bin/env python3
import h5py


def upgrade_old_network_file(file_name):

    """ Adds morphologyID filed if missing to network_file """

    print(f"Opening {file_name}")
    with h5py.File(file_name, 'a') as f:

        neuron_group = f["network/neurons"]

        if "morphologyID" not in neuron_group:
            print("Adding morphologyID data to file, default value 0.")
            n_neurons = len(f["network/neurons/morphology"])

            neuron_morph_id_data = neuron_group.create_dataset("morphologyID", (n_neurons,), "int", compression="gzip")
            neuron_morph_id_data[:] = 0

        else:
            print("morphologyID already exists, ignoring")

        if "morphologyKey" not in neuron_group:
            print("Adding morphologyKey to file, default value empty")
            n_neurons = len(f["network/neurons/morphology"])

            neuron_morph_key_data = neuron_group.create_dataset("morphologyKey", n_neurons, "S5", compression="gzip")

        if "parameterKey" not in neuron_group:
            print("Adding parameterKey to file, default value empty")
            n_neurons = len(f["network/neurons/morphology"])

            neuron_morph_key_data = neuron_group.create_dataset("parameterKey", n_neurons, "S5", compression="gzip")

        if "modulationKey" not in neuron_group:
            print("Adding morphologyKey to file, default value empty")
            n_neurons = len(f["network/neurons/morphology"])

            neuron_morph_key_data = neuron_group.create_dataset("modulationKey", n_neurons, "S5", compression="gzip")


if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser(description="Upgrade snudda network file (hdf5), adding missing morphologyID field")
    parser.add_argument("network_file", help="Network file (hdf5)", type=str)
    args = parser.parse_args()

    upgrade_old_network_file(args.network_file)

