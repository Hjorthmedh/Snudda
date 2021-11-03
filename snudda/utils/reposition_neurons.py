import h5py


class RepositionNeurons(object):

    """ Reposition neurons in the network file. """

    def __init__(self, position_file):

        """ Constructor.

        Args:
            position_file : network position file
        """

        self.position_file = position_file
        self.hdf5_file = h5py.File(position_file, "r+")

        assert "synapses" not in self.hdf5_file["network"], "This code should not be used after detection."

    def __del__(self):
        self.close()

    def close(self):

        """ Close file."""

        if self.hdf5_file:
            self.hdf5_file.close()
            self.hdf5_file = None

    def place(self, neuron_id, position, rotation=None, verbose=True):

        """ Reposition neuron neuron_id to position with rotation.

        Args:
            neuron_id : Neuron ID affected
            position : New position of neuron
            rotation : New rotatoin of neuron
            verbose : Print info?
        """

        info = (f"Moving neuron {neuron_id} from {self.hdf5_file['network/neurons/position'][neuron_id]} " 
                f"to {position}")
        if rotation is not None:
            info += f" setting rotation to {rotation}"

        if verbose:
            print(info)

        self.hdf5_file["network/neurons/position"][neuron_id] = position

        if rotation is not None:
            self.hdf5_file["network/neurons/rotation"][neuron_id] = rotation.reshape((9, ))

    def set_morphology_id(self, neuron_id, morphology_id):

        """ Set morphologyID for neuron with neuron_id (neuron_id = None means all neurons) """
        if neuron_id:
            self.hdf5_file["network/neurons/morphologyID"][neuron_id] = morphology_id
        else:
            self.hdf5_file["network/neurons/morphologyID"][:] = morphology_id

    def set_parameter_id(self, neuron_id, parameter_id):

        """ Set parameterID for neuron with neuron_id (neuron_id = None means all neurons) """

        if neuron_id:
            self.hdf5_file["network/neurons/parameterID"][neuron_id] = parameter_id
        else:
            self.hdf5_file["network/neurons/parameterID"][:] = parameter_id

    def set_modulation_id(self, neuron_id, modulation_id):

        """ Set modulationID for neuron with neuron_id (neuron_id = None means all neurons) """

        if neuron_id:
            self.hdf5_file["network/neurons/modulationID"][neuron_id] = modulation_id
        else:
            self.hdf5_file["network/neurons/modulationID"][:] = modulation_id
