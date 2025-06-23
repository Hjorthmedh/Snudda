#!/usr/bin/env python3

import h5py
import numpy as np
from lfpykit import CellGeometry, PointSourcePotential
import matplotlib.pyplot as plt

# You need to install lfpykit for this part

class LFP:

    def __init__(self, output_file):

        self.file = h5py.File(output_file, "r")

        self.extracellular_electrode_parameters = None

    def get_data(self, neuron_id):

        geometry = self.file["neurons"][str(neuron_id)]["geometry"][()].copy()
        membrane_current = self.file["neurons"][str(neuron_id)]["membrane.i_membrane_"][()]

        return geometry, membrane_current

    def set_electrode(self, sigma, x, y, z):

        self.extracellular_electrode_parameters = {
            "sigma": sigma,
            "x": np.array(x)*1e6,
            "y": np.array(y)*1e6,
            "z": np.array(z)*1e6   # micrometers for LFPy
        }

    def calculate_potential(self, neuron_id):

        geometry, membrane_current = self.get_data(neuron_id=neuron_id)

        cell_geometry = CellGeometry(x=geometry[:, [0, 1]]*1e6,
                                     y=geometry[:, [2, 3]]*1e6,
                                     z=geometry[:, [4, 5]]*1e6,
                                     d=geometry[:, [6, 7]]*1e6)

        forward_model = PointSourcePotential(cell_geometry, **self.extracellular_electrode_parameters)

        M = forward_model.get_transformation_matrix()
        V_e = M @ (membrane_current * 1e9)

        return V_e * 1e-3   # Convert to SI units

    def calculate_potential_network(self):

        V_e_list = []
        for neuron_id in self.file["neurons"].keys():
            V_e = self.calculate_potential(neuron_id=neuron_id)
            V_e_list.append(V_e)

        if len(V_e_list) > 1:
            return sum(V_e_list)
        else:
            return V_e

    def get_time(self):

        t = self.file["time"][()].copy()
        return t

    def plot_lfp(self, v_e, t):

        plt.plot(t, v_e.T)
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
        plt.ion()
        plt.show()


def cli():

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Generate LFP")
    parser.add_argument("simulation_file")

    args = parser.parse_args()

    lfp = LFP(output_file=args.simulation_file)

    lfp.set_electrode(sigma=0.3, x=[0], y=[0], z=[0])
    v_e = lfp.calculate_potential_network()
    t = lfp.get_time()

    lfp.plot_lfp(v_e, t)

    input("Press a key to continue")

    import pdb
    pdb.set_trace()

if __name__ == "__main__":
    cli()