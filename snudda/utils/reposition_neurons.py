import h5py


class RepositionNeurons(object):

    def __init__(self, position_file):

        self.position_file = position_file
        self.hdf5_file = h5py.File(position_file, "r+")

        assert "synapses" not in self.hdf5_file["network"], "This code should not be used after detection."

    def __del__(self):
        self.close()

    def close(self):
        if self.hdf5_file:
            self.hdf5_file.close()
            self.hdf5_file = None

    def place(self, neuron_id, position, rotation=None, verbose=True):

        info = (f"Moving neuron {neuron_id} from {self.hdf5_file['network/neurons/position'][neuron_id]} " 
                f"to {position}")
        if rotation is not None:
            info += f" setting rotation to {rotation}"

        if verbose:
            print(info)

        self.hdf5_file["network/neurons/position"][neuron_id] = position

        if rotation is not None:
            self.hdf5_file["network/neurons/rotation"][neuron_id] = rotation.reshape((9, ))

