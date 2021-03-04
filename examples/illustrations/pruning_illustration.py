import os

import numpy as np
from snudda.create_cube_mesh import create_cube_mesh
from snudda.detect import SnuddaDetect
from snudda.place import SnuddaPlace
from snudda.plotting.plot_network import PlotNetwork
from snudda.prune import SnuddaPrune
from snudda.utils.reposition_neurons import RepositionNeurons


class PruningIllustration(object):

    def __init__(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        self.network_path = "pruning_illustration_network"
        self.config_file = os.path.join(self.network_path, "network-config.json")
        self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        create_cube_mesh(file_name=os.path.join(self.network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        sp = SnuddaPlace(config_file=self.config_file, d_view=None)

        print("Calling read_config")
        sp.parse_config()
        print("Read done")
        sp.write_data(self.position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 200 micrometer distance from centre of soma.

        self.sd = SnuddaDetect(config_file=self.config_file, position_file=self.position_file,
                               save_file=self.save_file, rc=None,
                               hyper_voxel_size=150)

        # Reposition the neurons so we know how many synapses and where they will be located before pruning
        neuron_positions = np.array([[0, 59, 0],  # Postsynaptiska
                                     [0, 89, 0],
                                     [0, 119, 0],
                                     [0, 149, 0],
                                     [0, 179, 0],
                                     [0, 209, 0],
                                     [0, 239, 0],
                                     [0, 269, 0],
                                     [0, 299, 0],
                                     [0, 329, 0],
                                     [59, 0, 0],  # Presynaptiska
                                     [89, 0, 0],
                                     [119, 0, 0],
                                     [149, 0, 0],
                                     [179, 0, 0],
                                     [209, 0, 0],
                                     [239, 0, 0],
                                     [269, 0, 0],
                                     [299, 0, 0],
                                     [329, 0, 0],
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

        # Also update so that the new positions are saved in the place file
        rn = RepositionNeurons(self.position_file)
        for neuron_info in self.sd.neurons:
            rn.place(neuron_info["neuronID"], position=neuron_info["position"], rotation=neuron_info["rotation"],
                     verbose=False)
        rn.close()

        if False:
            self.sd.process_hyper_voxel(1)
            plt, ax = self.sd.plot_hyper_voxel(plot_neurons=True, elev_azim=(90, 0),
                                               draw_axon_voxels=False, draw_dendrite_voxels=False,
                                               draw_axons=True, draw_dendrites=True,
                                               show_axis=False, title="No pruning",
                                               fig_file_name="Pruning-fig-1-no-pruning")
            import pdb
            pdb.set_trace()

    def prune_network(self, pruning_config=None, fig_name=None, title=None):

        work_log = os.path.join(self.network_path, "log", "network-detect-worklog.hdf5")
        pruned_output = os.path.join(self.network_path, "network-pruned-synapses.hdf5")

        if pruning_config is not None and not os.path.exists(pruning_config):
            pruning_config = os.path.join(self.network_path, pruning_config)

        sp = SnuddaPrune(network_path=self.network_path, config_file=pruning_config)  # Use default config file
        sp.prune(pre_merge_only=False)
        sp = []

        plot_axon = True
        plot_dendrite = True
        #plot_axon = np.ones((20,), dtype=bool)
        #plot_dendrite = np.ones((20,), dtype=bool)
        #plot_axon[:10] = False
        #plot_dendrite[10:] = False

        pn = PlotNetwork(pruned_output)
        plt, ax = pn.plot(fig_name=fig_name, show_axis=False,
                          plot_axon=plot_axon, plot_dendrite=plot_dendrite,
                          title=title, title_pad=-14,
                          elev_azim=(90, 0))

        # Load the pruned data and check it
        # sl = SnuddaLoad(pruned_output)

if __name__ == "__main__":

    pil = PruningIllustration()
    pil.prune_network(fig_name="0-No-Pruning.pdf", title="No pruning")

    # Showing f1 pruning
    pil.prune_network(pruning_config="network-config-f1-1.json",
                      title="f1 = 1",
                      fig_name="1-Pruning-f1-1-No-pruning.pdf")

    pil.prune_network(pruning_config="network-config-f1-0.75.json",     # Not used in figure anymore
                      title="f1 = 0.75",
                      fig_name="1-Pruning-f1-0.75.pdf")

    pil.prune_network(pruning_config="network-config-f1-0.5.json", 
                      title="f1 = 0.5",
                      fig_name="1-Pruning-f1-0.5.pdf")

    pil.prune_network(pruning_config="network-config-f1-0.25.json",
                      title="f1 = 0.25",
                      fig_name="1-Pruning-f1-0.25.pdf")

    # Adding mu2 pruning
    pil.prune_network(pruning_config="network-config-f1-1-mu2-3.json",
                      title="f1 = 1, mu2 = 3",
                      fig_name="2-Pruning-f1-1-mu2-3.pdf")

    pil.prune_network(pruning_config="network-config-f1-0.75-mu2-3.json",   # Not used in figure anymore
                      title="f1 = 0.75, mu2 = 3",
                      fig_name="2-Pruning-f1-0.75-mu2-3.pdf")

    pil.prune_network(pruning_config="network-config-f1-0.5-mu2-3.json",
                      title="f1 = 0.5, mu2 = 3",
                      fig_name="2-Pruning-f1-0.5-mu2-3.pdf")

    pil.prune_network(pruning_config="network-config-f1-0.25-mu2-3.json",
                      title="f1 = 0.25, mu2 = 3",
                      fig_name="2-Pruning-f1-0.25-mu2-3.pdf")


    # Softmax
    pil.prune_network(pruning_config="network-config-softmax-4.json",   # Not used in figure anymore
                      title="softmax=4",
                      fig_name="3-Pruning-softmax-4.pdf")

    pil.prune_network(pruning_config="network-config-softmax-3.json",
                      title="softmax=3",
                      fig_name="3-Pruning-softmax-3.pdf")

    pil.prune_network(pruning_config="network-config-softmax-2.json",
                      title="softmax=2",
                      fig_name="3-Pruning-softmax-2.pdf")

    pil.prune_network(pruning_config="network-config-softmax-1.json",
                      title="softmax=1",
                      fig_name="3-Pruning-softmax-1.pdf")

    # Adding a3 pruning
    pil.prune_network(pruning_config="network-config-a3-1.json",    # Not used in figure anymore
                      title="a3 = 1",
                      fig_name="4-Pruning-a3-1.pdf")

    pil.prune_network(pruning_config="network-config-a3-0.75.json",
                      title="a3 = 0.75",
                      fig_name="4-Pruning-a3-0.75.pdf")

    pil.prune_network(pruning_config="network-config-a3-0.5.json",
                      title="a3 = 0.5",
                      fig_name="4-Pruning-a3-0.5.pdf")

    pil.prune_network(pruning_config="network-config-a3-0.25.json",
                      title="a3 = 0.25",
                      fig_name="4-Pruning-a3-0.25.pdf")


    # Distance dependent examples
    pil.prune_network(pruning_config="network-config-distance-dep-v3-1.json",
                      title='"exp(-(0.5*d/70e-6)^2)"',
                      fig_name="5-Pruning-distance-dep-1.pdf")

    pil.prune_network(pruning_config="network-config-distance-dep-v3-2.json",
                      title='"exp(-(d-190e-6)**2/(2*40e-6**2))"',
                      fig_name="5-Pruning-distance-dep-2.pdf")

    pil.prune_network(pruning_config="network-config-distance-dep-v3-3.json",
                      title='"1/(1+exp(-(d-280e-6)/30e-6))"',
                      fig_name="5-Pruning-distance-dep-3.pdf")


    print(f"\n --> Figures written to {os.path.join(pil.network_path,'figures')}\n")

    #import pdb
    #pdb.set_trace()

