import unittest
import os

import neuron
import numpy as np

from snudda.utils import snudda_parse_path
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel

import bluepyopt.ephys as ephys
from snudda.neurons.neuron_model_extended import NeuronModel
from snudda.neurons.neuron_morphology import NeuronMorphology


class SegmentIdTestCase(unittest.TestCase):

    """NEURON renumberes all the branches into segments with separate id. This code verifies that Snudda's algorithm
       replicates what NEURON does internally. It will also warn if NEURON's algorithm at some point changes.

        If you want to test all morphologies in BasalGangliaData, first set SNUDDA_DATA:

        export SNUDDA_DATA=/home/hjorth/HBP/BasalGangliaData/data/
        python3 -m unittest test_segmentid.py

    """

    def setUp(self) -> None:
        self.sim = NrnSimulatorParallel(cvode_active=False)

    def test_segment_id_numbering(self, morph_file=None):

        # Load morphology into Snudda
        if not morph_file:
            print("No morphology file specified, skipping.")
            return

        print(f"Loading neuron {morph_file}")
        snudda_neuron = NeuronMorphology(name="fs", swc_filename=morph_file, use_cache=False)

        self.assertTrue((snudda_neuron.soma[0, :3] == 0).all(), f"Soma should be centered for {morph_file}.")

        # Load morphology into NEURON
        neuron_model = NeuronModel(param_file=os.path.join("data", "fake-parameters.json"),
                                   morph_path=morph_file,
                                   mech_file=os.path.join("data", "fake-mechanisms.json"),
                                   cell_name="fs",
                                   modulation_file=None,
                                   parameter_id=0,
                                   modulation_id=0)

        neuron_model.instantiate(sim=self.sim)

        ax = None

        for link, sec_id, sec_x, dend_sec \
            in zip(snudda_neuron.dend_links,
                   snudda_neuron.dend_sec_id,
                   snudda_neuron.dend_sec_x,
                   neuron_model.map_id_to_compartment(section_id=snudda_neuron.dend_sec_id)):

            assert sec_id > 0, f"sec id {sec_id}, for dendrites should be > 0"

            # Coordinates of segment in snudda NeuronMorphology -- convert to natural units micrometers for NEURON
            x0, y0, z0 = snudda_neuron.dend[link[0], :3] * 1e6
            x1, y1, z1 = snudda_neuron.dend[link[1], :3] * 1e6

            # Coordinates of segment in NEURON

            n_points = int(dend_sec.n3d())
            arc_dist = np.array([dend_sec.arc3d(x) for x in range(0, n_points)])
            norm_arc_dist = arc_dist / arc_dist[-1]

            # Find closest point
            closest_idx0 = np.argmin(np.abs(norm_arc_dist - sec_x[0]))
            closest_idx1 = np.argmin(np.abs(norm_arc_dist - sec_x[1]))

            x0_ref = dend_sec.x3d(closest_idx0)
            y0_ref = dend_sec.y3d(closest_idx0)
            z0_ref = dend_sec.z3d(closest_idx0)

            x1_ref = dend_sec.x3d(closest_idx1)
            y1_ref = dend_sec.y3d(closest_idx1)
            z1_ref = dend_sec.z3d(closest_idx1)

            error_cutoff = 10

            self.assertTrue(np.linalg.norm([x0-x0_ref, y0-y0_ref, z0-z0_ref]) < error_cutoff
                             and np.linalg.norm([x1-x1_ref, y1-y1_ref, z1-z1_ref]) < error_cutoff,
                            (f"Error when parsing {morph_file}\n"
                             f"Snudda morphology sec_id {sec_id}, sec_x {sec_x[0]} to {sec_x[1]} "
                             f"xyz = {x0}, {y0}, {z0} to {x1}, {y1}, {z1}\n"
                             f"NEURON coords {x0_ref}, {y0_ref}, {z0_ref} to {x1_ref}, {y1_ref}, {z1_ref}\n"
                             f"Distance: {np.linalg.norm([x0-x0_ref, y0-y0_ref, z0-z0_ref])} "
                             f"and {np.linalg.norm([x1-x1_ref, y1-y1_ref, z1-z1_ref])}"))

    def test_neurons_in_folder(self, neuron_dir=None):

        import glob

        if not neuron_dir:
            print("No neuron dir given, skipping.")
            return

        n_dirs = glob.glob(os.path.join(neuron_dir, '*'))

        # In case the user gave the neuron directory with SWC files
        swc_files = glob.glob(os.path.join(neuron_dir, '*swc'))

        for n_dir in n_dirs:

            swc_files1 = glob.glob(os.path.join(neuron_dir, n_dir, '*swc'))
            for swc_f in swc_files1:
                swc_files.append(os.path.join(neuron_dir, n_dir, swc_f))

            swc_files2 = glob.glob(os.path.join(neuron_dir, n_dir, "morphology", '*swc'))
            for swc_f in swc_files2:
                swc_files.append(os.path.join(neuron_dir, n_dir, "morphology", swc_f))

        for swc_file in swc_files:
            with self.subTest(msg=f"Testing {swc_file}"):
                self.test_segment_id_numbering(morph_file=swc_file)

    def test_all_dir(self):

        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "..", "snudda", "data")

        neuron_dirs = [snudda_parse_path(os.path.join("$SNUDDA_DATA", "neurons", "striatum", "dspn")),
                       snudda_parse_path(os.path.join("$SNUDDA_DATA", "neurons", "striatum", "ispn")),
                       snudda_parse_path(os.path.join("$SNUDDA_DATA", "neurons", "striatum", "fs")),
                       snudda_parse_path(os.path.join("$SNUDDA_DATA", "neurons", "striatum", "lts")),
                       ]

        for neuron_dir in neuron_dirs:
            self.test_neurons_in_folder(neuron_dir=neuron_dir)


if __name__ == '__main__':
    unittest.main()
