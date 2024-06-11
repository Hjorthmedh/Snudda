import unittest
import os
import time


class TestProject(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        from snudda.place.create_cube_mesh import create_cube_mesh

        # Create cube meshes
        self.network_path = os.path.join("networks", "network_testing_project")
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
        proj_file = os.path.join("data", "ExampleProjection.json")

        cnc.neuron_projection(neuron_name="dSPN",
                              target_name="iSPN",
                              projection_name="ExampleProjection",
                              projection_file=proj_file,
                              source_volume="VolumeA",
                              dest_volume="VolumeB",
                              projection_radius=100e-6,
                              number_of_targets=[10, 5],
                              number_of_synapses=[10, 5],
                              dendrite_synapse_density="1",
                              connection_type="GABA",
                              dist_pruning=None,
                              f1=0.9,
                              soft_max=None,
                              mu2=None,
                              a3=None)

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
                              region_name="VolumeA",
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
                              region_name="VolumeB",
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

        tot_dd_syn_ctr = 0
        for dspn_id in dspn_id_list:
            for dspn_id2 in dspn_id_list:

                synapses, synapse_coords = sl.find_synapses(pre_id=dspn_id, post_id=dspn_id2)
                if synapses is not None:
                    tot_dd_syn_ctr += synapses.shape[0]

        tot_ii_syn_ctr = 0
        for ispn_id in ispn_id_list:
            for ispn_id2 in ispn_id_list:

                synapses, synapse_coords = sl.find_synapses(pre_id=ispn_id, post_id=ispn_id2)
                if synapses is not None:
                    tot_ii_syn_ctr += synapses.shape[0]

        with self.subTest(stage="normal_synapses_exist"):
            # In this toy example neurons are quite sparsely placed, but we should have at least some
            # synapses
            self.assertTrue(tot_dd_syn_ctr > 0)
            self.assertTrue(tot_ii_syn_ctr > 0)

        # We need to run in parallel also to verify we get same result (same random seed)

        serial_synapses = sl.data["synapses"].copy()
        del sl  # Close old file so we can overwrite it

        os.environ["IPYTHONDIR"] = os.path.join(os.path.abspath(os.getcwd()), ".ipython")
        os.environ["IPYTHON_PROFILE"] = "default"
        os.system("ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&")
        time.sleep(15)

        # Run place, detect and prune in parallel by passing rc
        from ipyparallel import Client
        u_file = os.path.join(".ipython", "profile_default", "security", "ipcontroller-client.json")
        rc = Client(url_file=u_file, timeout=120, debug=False)
        d_view = rc.direct_view(targets='all')  # rc[:] # Direct view into clients

        from snudda.detect.detect import SnuddaDetect
        sd = SnuddaDetect(network_path=self.network_path, hyper_voxel_size=100, rc=rc, verbose=True)
        sd.detect()

        from snudda.detect.project import SnuddaProject
        # TODO: Currently SnuddaProject only runs in serial
        sp = SnuddaProject(network_path=self.network_path)
        sp.project()

        from snudda.detect.prune import SnuddaPrune
        # Prune has different methods for serial and parallel execution, important to test it!
        sp = SnuddaPrune(network_path=self.network_path, rc=rc, verbose=True)
        sp.prune()

        with self.subTest(stage="check-parallel-identical"):
            sl2 = SnuddaLoad(network_file)
            parallel_synapses = sl2.data["synapses"].copy()

            # ParameterID, sec_X etc are randomised in hyper voxel, so you need to use same
            # hypervoxel size for reproducability between serial and parallel execution

            # All synapses should be identical regardless of serial or parallel execution path
            self.assertTrue((serial_synapses == parallel_synapses).all())

        os.system("ipcluster stop")


if __name__ == '__main__':
    unittest.main()
