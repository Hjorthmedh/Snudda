import numpy as np

class NEURONNeuronExtractor:

    def __init__(self, sim):

        self.sim = sim

    def extract_geometry_for_neuron(self, neuron_id):

        loc_info = []

        for sec_id, sec in self.sim.neurons[neuron_id].section_lookup.items():

            if sec.n3d() > 1:
                loc_info.append(np.array([[sec.x3d(i), sec.x3d(i+1),
                                           sec.y3d(i), sec.y3d(i+1),
                                           sec.z3d(i), sec.z3d(i+1),
                                           sec.diam3d(i), sec.diam3d(i+1)]
                                          for i in range(sec.n3d()-1)]))

        geometry = np.vstack(loc_info) * 1e-6  # We save in meter

        return geometry

    def extract_coordinates_for_neuron(self, neuron_id):

        loc_info = []

        for sec_id, sec in self.sim.neurons[neuron_id].section_lookup.items():

            if sec.n3d() > 1:
                loc_info.append(np.array([[sec.x3d(i), sec.y3d(i), sec.z3d(i), sec.diam3d(i)]
                                          for i in range(sec.n3d())]))

        coordinates = np.vstack(loc_info) * 1e-6  # We save in meter

        return coordinates

    def write_all_coordinates_to_hdf5(self, hdf5_file):

        print(f"Writing geometry to {hdf5_file.file.filename}")

        for neuron_id in self.sim.neurons.keys():

            geometry = self.extract_geometry_for_neuron(neuron_id=neuron_id)
            hdf5_file["neurons"][str(neuron_id)].create_dataset("geometry", data=geometry)

