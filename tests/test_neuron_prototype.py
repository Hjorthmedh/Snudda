import json
import os
import unittest
import numpy as np

from snudda.neurons.neuron_prototype import NeuronPrototype


class NeuronPrototypeTestCase(unittest.TestCase):

    def setUp(self) -> None:

        pass

    def test_setup(self):

        neuron_path = os.path.join("validation", "striatum-var", "dspn",
                                   "str-dspn-e150917_c9_d1-mWT-1215MSN03-v20210212")

        with self.subTest(msg="neuron_path given"):
            np1 = NeuronPrototype(neuron_path=neuron_path,
                                  neuron_name="np1")

            nm1 = np1.clone(parameter_id=123, morphology_id=456, modulation_id=789)

            # Check if we load the right neuron
            self.assertEqual(os.path.basename(nm1.swc_filename), "WT-1215MSN03-cor-rep-ax-res3-var0.swc")

            par_file = os.path.join(neuron_path, "parameters.json")
            with open(par_file, "r") as f:
                par_file_data = json.load(f)

            self.assertEqual(np1.get_parameters(3), par_file_data["3"])
            self.assertEqual(np1.get_parameters(7), par_file_data["3"])  # Parameter ID is modulo

            mod_file = os.path.join(neuron_path, "modulation.json")
            with open(mod_file, "r") as f:
                mod_file_data = json.load(f)

            self.assertEqual(np1.get_modulation_parameters(2), mod_file_data[2])

        morphology_path = os.path.join(neuron_path, "morphology")
        parameter_path = os.path.join(neuron_path, "parameters.json")
        modulation_path = os.path.join(neuron_path, "modulation.json")
        mechanism_path = os.path.join(neuron_path, "mechanisms.json")
        meta_path = os.path.join(neuron_path, "meta.json")

        with self.subTest(msg="other paths given"):
            np2 = NeuronPrototype(morphology_path=morphology_path,
                                  parameter_path=parameter_path,
                                  mechanism_path=mechanism_path,
                                  modulation_path=modulation_path,
                                  meta_path=meta_path,
                                  neuron_name="np2", neuron_path=None)

            nm2 = np1.clone(parameter_id=124, morphology_id=10, modulation_id=11)

            # Check correct morphology was picked
            self.assertEqual(os.path.basename(nm2.swc_filename), "WT-1215MSN03-cor-rep-ax-res3-var5.swc")

        with self.subTest(msg="instantiate all test"):
            np3 = NeuronPrototype(neuron_path=neuron_path,
                                  neuron_name="np3")

            np3.instantiate()
            res = np3.apply("compartment_length", ["dend"])

            # Verify that correct
            res2 = []
            for nm in np3.morphology_cache.values():
                res2.append(nm.compartment_length("dend"))

            res_sum = np.sort([sum(x) for x in res])
            res2_sum = np.sort([sum(x) for x in res2])
            self.assertTrue((res_sum == res2_sum).all())

            print(f"Compartment length sums: {res_sum} vs {res2_sum}")


if __name__ == '__main__':
    unittest.main()
