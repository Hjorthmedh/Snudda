import os

import h5py

from snudda.utils.load import SnuddaLoad


class InspectInput(object):

    def __init__(self, network_file, input_file):

        self.network_file = network_file
        self.input_file = input_file

        # We just need this to identify which neuron is which
        self.network = SnuddaLoad(self.network_file, load_synapses=False)

        self.input_data = h5py.File(input_file, 'r')

    def get_morphologies(self):

        return [os.path.basename(c["morphology"]) for c in self.network.data["neurons"]]

    def get_input_types(self, cell_id):

        s_cell_id = [str(c) for c in cell_id]
        input_types = set()

        for sc_id in s_cell_id:
            inp_t = set([inp for inp in self.input_data["input"][sc_id]])
            input_types = input_types.union(inp_t)

        return list(input_types)

    def check_input_ratio(self, neuron_type, verbose=True):

        print(f"Counting inputs for {neuron_type}")

        cell_id = self.network.get_neuron_id_of_type(neuron_type)
        cell_id_str = [str(c) for c in cell_id]

        input_count = dict()

        cell_morph = self.get_morphologies()

        unique_morph = set([cell_morph[c] for c in cell_id])

        input_type_list = self.get_input_types(cell_id)

        for inp in input_type_list:
            input_count[inp] = dict()
            for um in unique_morph:
                input_count[inp][um] = 0

        morph_counter = dict()

        for um in unique_morph:
            morph_counter[um] = 0

        # !!! TODO: We should split this by morphology also...

        for c_id in self.input_data['input']:
            if c_id in cell_id_str:
                morph = cell_morph[int(c_id)]
                morph_counter[morph] += 1

                for inp in input_type_list:
                    if inp in self.input_data['input'][c_id]:
                        n_input = len(self.input_data['input'][c_id][inp]['nSpikes'])
                        input_count[inp][morph] += n_input

        for inp in input_type_list:
            for um in unique_morph:
                avg_inp = round(input_count[inp][um] / morph_counter[um], 1)
                print(f"{inp} morphology {um}: {avg_inp} inputs")

        return input_count


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Inspect input")
    parser.add_argument("networkFile", help="Network file (hdf5)", type=str)
    parser.add_argument("inputFile", help="Input file (hdf5)", type=str)
    parser.add_argument("neuronType", help="Neuron type", type=str)

    args = parser.parse_args()

    inspector = InspectInput(network_file=args.networkFile,
                             input_file=args.inputFile)

    inspector.check_input_ratio(neuron_type=args.neuronType)
