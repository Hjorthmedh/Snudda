import os
import unittest

from snudda.create_cube_mesh import create_cube_mesh
from snudda.detect import SnuddaDetect
from snudda.load import SnuddaLoad
from snudda.place import SnuddaPlace
import numpy as np

from snudda.prune import SnuddaPrune


class TestPrune(unittest.TestCase):
    
    def setUp(self):

        os.chdir(os.path.dirname(__file__))
        self.network_path = os.path.join(os.path.dirname(__file__), "tests", "network_testing_prune")

        create_cube_mesh(file_name=os.path.join(self.network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        config_file = os.path.join(self.network_path, "network-config.json")
        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        sp = SnuddaPlace(config_file=config_file, d_view=None)

        sp.read_config()
        sp.write_data(position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file,
                               save_file=save_file, rc=None,
                               hyper_voxel_size=120)

        # Reposition the neurons so we know how many synapses and where they will be located before pruning
        neuron_positions = np.array([[0, 20, 0],  # Postsynaptiska
                                     [0, 40, 0],
                                     [0, 60, 0],
                                     [0, 80, 0],
                                     [0, 100, 0],
                                     [0, 120, 0],
                                     [0, 140, 0],
                                     [0, 160, 0],
                                     [0, 180, 0],
                                     [0, 200, 0],
                                     [20, 0, 0],  # Presynaptiska
                                     [40, 0, 0],
                                     [60, 0, 0],
                                     [80, 0, 0],
                                     [100, 0, 0],
                                     [120, 0, 0],
                                     [140, 0, 0],
                                     [160, 0, 0],
                                     [180, 0, 0],
                                     [200, 0, 0],
                                     ]) * 1e-6
        
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
        
    def test_detect(self): 

        work_log = os.path.join(self.network_path, "log", "network-detect-worklog.hdf5")

        sp = SnuddaPrune(work_history_file=work_log)
        sp.prune(pre_merge_only=False)

        # Load the pruned data and check it

        pruned_output = os.path.join(self.network_path, "network-pruned-synapses.hdf5")
        sl = SnuddaLoad(pruned_output)

        # TODO: Add pruning tests
        # TODO: Add option to set pruning rules after detection, so we can rerun different kinds quickly
        pass
