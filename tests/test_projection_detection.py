import os
import time
import unittest


class TestProjectionDetection(unittest.TestCase):

    def setUp(self):
        from snudda.place.create_cube_mesh import create_cube_mesh

        # Create cube meshes
        self.network_path = os.path.join("networks", "network_testing_projection_detection")
        mesh_file_a = os.path.join(self.network_path, "mesh", "volume_A.obj")
        mesh_file_b = os.path.join(self.network_path, "mesh", "volume_B.obj")

        create_cube_mesh(mesh_file_a, [5e-3, 0, 0], 300e-6, "Volume A - connect structures example")
        create_cube_mesh(mesh_file_b, [-5e-3, 0, 0], 300e-6, "Volume B - connect structures example")

        # Define network
        from snudda.init.init import SnuddaInit

        cnc = SnuddaInit(network_path=self.network_path, random_seed=123)
        cnc.define_structure(struct_name="VolumeA", struct_mesh=mesh_file_a, d_min=15e-6, mesh_bin_width=50e-6)
        cnc.define_structure(struct_name="VolumeB", struct_mesh=mesh_file_b, d_min=15e-6, mesh_bin_width=50e-6)

        cnc.add_neurons(name="dSPN", num_neurons=20, volume_id="VolumeA",
                        neuron_dir=os.path.join("$DATA", "neurons", "striatum", "dspn"))
        cnc.add_neurons(name="iSPN", num_neurons=20, volume_id="VolumeB",
                        neuron_dir=os.path.join("$DATA", "neurons", "striatum", "ispn"))

        # Add the projection we want to test dSPN->iSPN
        proj_file = os.path.join("data", "ExampleProjectionDetection.json")

        cnc.add_neuron_target(neuron_name="dSPN",
                              target_name="iSPN",
                              connection_type="GABA",
                              dist_pruning=None,
                              f1=None,
                              soft_max=None,
                              mu2=None,
                              a3=None,
                              conductance=5e-9,
                              cluster_synapses=False,
                              mod_file="tmGlutA",
                              parameter_file=None,
                              channel_param_dictionary=None,
                              projection_file=proj_file)

        # Also add dSPN-dSPN and iSPN-iSPN synapses
        # Note we do NOT add dSPN-iSPN again this way, as that would overwrite the above connections
        # (The above neuron_projection will also do normal touch detection)

        SPN2SPNdistDepPruning = "1-exp(-(0.4*d/60e-6)**2)"

        MSD1gGABA = [0.24e-9, 0.1e-9]
        MSD2gGABA = [0.24e-9, 0.1e-9]

        MSD1GABAfailRate = 0.7  # Taverna 2008, figure 2
        MSD2GABAfailRate = 0.4  # Taverna 2008, 2mM

        pfdSPNdSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-DD-tmgaba-fit.json")
        pfdSPNiSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-DI-tmgaba-fit.json")
        pfiSPNdSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-ID-tmgaba-fit.json")
        pfiSPNiSPN = os.path.join("$DATA", "synapses", "striatum", "PlanertFitting-II-tmgaba-fit.json")

        cnc.add_neuron_target(neuron_name="dSPN",
                              target_name="dSPN",
                              connection_type="GABA",
                              dist_pruning=SPN2SPNdistDepPruning,
                              f1=0.38, soft_max=3, mu2=2.4,
                              a3=1.0,
                              conductance=MSD1gGABA,
                              parameter_file=pfdSPNdSPN,
                              mod_file="tmGabaA",
                              channel_param_dictionary={"tau1": (1.3e-3, 1e3),
                                                        "tau2": (12.4e-3, 1e3),
                                                        "failRate": MSD1GABAfailRate})

        cnc.add_neuron_target(neuron_name="iSPN",
                              target_name="iSPN",
                              connection_type="GABA",
                              dist_pruning=SPN2SPNdistDepPruning,
                              f1=0.55, soft_max=4, mu2=2.4,
                              a3=1.0,
                              conductance=MSD2gGABA,
                              parameter_file=pfiSPNiSPN,
                              mod_file="tmGabaA",
                              channel_param_dictionary={"tau1": (1.3e-3, 1e3),
                                                        "tau2": (12.4e-3, 1e3),
                                                        "failRate": MSD2GABAfailRate})

        cnc.write_json()

        # Place neurons, then detect, project and prune

        from snudda.place.place import SnuddaPlace
        sp = SnuddaPlace(network_path=self.network_path, verbose=True)
        sp.parse_config()
        sp.write_data()

        from snudda.detect.detect import SnuddaDetect
        sd = SnuddaDetect(network_path=self.network_path, hyper_voxel_size=100, verbose=True)
        sd.detect()

        from snudda.detect.project import SnuddaProject
        sp = SnuddaProject(network_path=self.network_path)
        sp.project()

        from snudda.detect.prune import SnuddaPrune
        sp = SnuddaPrune(network_path=self.network_path, verbose=True)
        sp.prune()

    def test_project(self):

        # Are there connections dSPN->iSPN
        from snudda.utils.load import SnuddaLoad
        network_file = os.path.join(self.network_path, "network-synapses.hdf5")
        sl = SnuddaLoad(network_file)

        dspn_id_list = sl.get_neuron_id_of_type("dSPN")
        ispn_id_list = sl.get_neuron_id_of_type("iSPN")

        tot_proj_ctr = 0

        for dspn_id in dspn_id_list:
            for ispn_id in ispn_id_list:

                synapses, synapse_coords = sl.find_synapses(pre_id=dspn_id, post_id=ispn_id)
                if synapses is not None:
                    tot_proj_ctr += synapses.shape[0]

        with self.subTest(stage="projection_exists"):
            # There should be projection synapses between dSPN and iSPN in this toy example
            self.assertTrue(tot_proj_ctr > 0)


if __name__ == '__main__':
    unittest.main()
