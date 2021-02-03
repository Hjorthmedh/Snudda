import os

import numpy as np
from snudda.create_cube_mesh import create_cube_mesh
from snudda.detect import SnuddaDetect
from snudda.place import SnuddaPlace


class PruningIllustration(object):

    def __init__(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        self.network_path = os.path.join(os.path.dirname(__file__), "pruning_illustration_network")
        self.config_file = os.path.join(self.network_path, "network-config.json")
        self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        create_cube_mesh(file_name=os.path.join(self.network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        sp = SnuddaPlace(config_file=self.config_file, d_view=None)

        print("Calling read_config")
        sp.read_config()
        print("Read done")
        sp.write_data(self.position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 200 micrometer distance from centre of soma.

        self.sd = SnuddaDetect(config_file=self.config_file, position_file=self.position_file,
                               save_file=self.save_file, rc=None,
                               hyper_voxel_size=150)

        # Reposition the neurons so we know how many synapses and where they will be located before pruning
        neuron_positions = np.array([[0, 60-1, 0],  # Postsynaptiska
                                     [0, 84-1, 0],
                                     [0, 108-1, 0],
                                     [0, 132-1, 0],
                                     [0, 156-1, 0],
                                     [0, 180-1, 0],
                                     [0, 204-1, 0],
                                     [0, 228-1, 0],
                                     [0, 252-1, 0],
                                     [0, 276-1, 0],
                                     [60-1, 0, 0],  # Presynaptiska
                                     [84-1, 0, 0],
                                     [108-1, 0, 0],
                                     [132-1, 0, 0],
                                     [156-1, 0, 0],
                                     [180-1, 0, 0],
                                     [204-1, 0, 0],
                                     [228-1, 0, 0],
                                     [252-1, 0, 0],
                                     [276-1, 0, 0],
                                     ]) * 1e-6

        # TODO: Add potential for gap junctions also by having 5 + 5 neurons in other grid

        for idx, pos in enumerate(neuron_positions):
            self.sd.neurons[idx]["position"] = pos

        ang = -np.pi / 2
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])

        ang = np.pi / 2
        R_y = np.array([[np.cos(ang), 0, np.sin(ang)],
                        [0, 1, 0],
                        [-np.sin(ang), 0, np.cos(ang)]])

        for idx in range(0, 10):  # Post synaptic neurons
            self.sd.neurons[idx]["rotation"] = R_x

        for idx in range(10, 20):  # Presynaptic neurons
            self.sd.neurons[idx]["rotation"] = R_y

        self.sd.detect(restart_detection_flag=True)

        if True:
            self.sd.process_hyper_voxel(1)
            plt, ax = self.sd.plot_hyper_voxel(plot_neurons=True, elev_azim=(90, 0),
                                               draw_axon_voxels=False, draw_dendrite_voxels=False,
                                               draw_axons=True, draw_dendrites=True,
                                               show_axis=True, title="", fig_file_name="Pruning-fig-1")
            # TODO: Check why soma is plotted in wrong place? Mistake with origo plotoffset?
            import pdb
            pdb.set_trace()


# Create a network with double axon and double dendrites (tuning fork style)

    def setup_network(self):

        pass
# Position the neurons


if __name__ == "__main__":

    pil = PruningIllustration()
    import pdb
    pdb.set_trace()