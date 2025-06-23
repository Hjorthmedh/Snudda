import numpy as np

class NEURONNeuronExtractor:

    def __init__(self, sim):

        self.sim = sim

    def extract_geometry_for_neuron(self, neuron_id):

        loc_info = []

        for sec_id, sec in self.sim.neurons[neuron_id].section_lookup.items():

            if sec.n3d() <= 1:
                continue

            coords = np.array([[sec.x3d(i), sec.y3d(i), sec.z3d(i), sec.diam3d(i), sec.arc3d(i)]
                                for i in range(sec.n3d())])

            arc_points = np.linspace(0.0,1.0, num=sec.nseg+1, endpoint=True)

            xp = np.interp(arc_points, coords[:, 4], coords[:,0])
            yp = np.interp(arc_points, coords[:, 4], coords[:,1])
            zp = np.interp(arc_points, coords[:, 4], coords[:,2])
            dp = np.interp(arc_points, coords[:, 4], coords[:,3])

            for i in range(0, sec.nseg):
                loc_info.append(np.array([xp[i], xp[i+1],
                                          yp[i], yp[i+1],
                                          zp[i], zp[i+1],
                                          dp[i], dp[i+1]]))

        geometry = np.vstack(loc_info) * 1e-6  # We save in meter

        return geometry

    def extract_geometry_for_neuron_OLD(self, neuron_id):

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
                loc_info.append(np.array([[sec.x3d(i), sec.y3d(i), sec.z3d(i), sec.diam3d(i), sec.arc3d(i)]
                                          for i in range(sec.n3d())]))

        coordinates = np.vstack(loc_info) * 1e-6  # We save in meter

        return coordinates

    def write_all_coordinates_to_hdf5(self, hdf5_file):

        print(f"Writing geometry to {hdf5_file.file.filename}")

        for neuron_id in self.sim.neurons.keys():

            geometry = self.extract_geometry_for_neuron(neuron_id=neuron_id)
            hdf5_file["neurons"][str(neuron_id)].create_dataset("geometry", data=geometry)

