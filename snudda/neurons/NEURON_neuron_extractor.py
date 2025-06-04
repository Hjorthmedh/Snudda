import numpy as np
from snudda import SnuddaSimulate


class NEURONNeuronExtractor:

    def __init__(self, sim: SnuddaSimulate):

        self.sim = sim

    def extract_geometry_for_neuron(self, neuron_id):

        loc_info = []

        for sec_id, sec in self.sim.neurons[neuron_id].section_lookup.items():

            loc_info.append(np.array([[sec.x3d(i), sec.y3d(i), sec.z3d(i), sec.dim3d(i)]
                                      for i in range(sec.n3d)]))

        geometry = np.vstack(loc_info)

        return geometry

    def write_all_coordinates_to_hdf5(self, hdf5_file):

        print(f"Writing geometry to {hdf5_file.file.filename}")

        for neuron_id in self.sim.neurons.keys():

            geometry = self.extract_geometry_for_neuron(neuron_id=neuron_id)
            hdf5_file["neurons"][str(neuron_id)].create_dataset("geometry", geometry)
