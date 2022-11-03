import numpy as np
import h5py


class CloneNeurons:

    def __init__(self, source_network_file, dest_network_file, neuron_type):

        """ This code copies neurons of neuron_type from one network to another.
            The number of neurons must match in the two networks. """
        
        self.source_file = h5py.File(source_network_file,"r")
        self.dest_file = h5py.File(dest_network_file, "a")
        self.neuron_type = neuron_type

    def close(self):
        self.source_file.close()
        self.dest_file.close()
        
    def get_neuron_id(self, hdf5_file, neuron_type):

        neuron_id = [ x for x,y in zip(hdf5_file["network/neurons/neuronID"],
                                       hdf5_file["network/neurons/name"])
                      if y.decode().split("_")[0] == neuron_type]
        
        return np.array(neuron_id)

    def clone(self, source_id, dest_id):

        print(f"Cloning source_id={source_id} --> dest_id={dest_id}")
        
        for var_name in self.source_file["network/neurons"].keys():

            if var_name == "neuronID":
                # Keep original dest neuron ID
                continue

            if len(self.source_file["network/neurons"][var_name].shape) == 1:
                self.dest_file["network/neurons"][var_name][dest_id] = \
                    self.source_file["network/neurons"][var_name][source_id]
            else:
                self.dest_file["network/neurons"][var_name][dest_id, :] = \
                    self.source_file["network/neurons"][var_name][source_id, :]


    def clone_neurons(self, neuron_type=None):

        if neuron_type is None:
            neuron_type = self.neuron_type

        source_id = self.get_neuron_id(self.source_file, neuron_type)
        dest_id = self.get_neuron_id(self.dest_file, neuron_type)

        assert len(source_id) == len(dest_id), \
            f"There must be equal number of {neuron_type} in the two networks"
        
        for s_id, d_id in zip(source_id, dest_id):
            self.clone(source_id = s_id, dest_id = d_id)

            
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("source_file")
    parser.add_argument("dest_file")
    parser.add_argument("neuron_type")

    args = parser.parse_args()

    cn = CloneNeurons(source_network_file=args.source_file,
                      dest_network_file=args.dest_file,
                      neuron_type=args.neuron_type)

    cn.clone_neurons()
    cn.close()
