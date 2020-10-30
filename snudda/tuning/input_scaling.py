import os
import glob
import collections
import json
from snudda.CreateCubeMesh import CreateCubeMesh
from snudda.Neuron_morphology import NeuronMorphology

class InputScaling(object):

    def __init__(self, network_dir, cellspec_dir):

        self.network_dir = network_dir
        self.cellspec_dir = cellspec_dir

        if not os.path.isdir(self.network_dir):
            os.mkdir(self.network_dir)

        assert os.path.isdir(self.cellspec_dir), f"Cellspec directory {self.cellspec_dir} does not exist."

        self.network_file_name = os.path.join(self.network_dir,"network-config.json")

    # Writes config files
    def setup_network(self):

        config_def = self.create_config(num_replicas=10)

        print(f"Writing network config file to {self.network_file_name}")
        with open(self.network_file_name, "w") as f:
            json.dump(config_def, f, indent=2)

        ccm = CreateCubeMesh("data/mesh/InputTestMesh.obj",[ 0, 0, 0], 1e-3,
                             description="Mesh file used for Input Scaling")

    # This loops through all neuron directories in cellspec in preparation of writing a network config file
    def gather_all_neurons(self):
        all_neurons = collections.OrderedDict()

        neuron_type_dir = [d for d in glob.glob(os.path.join(self.cellspec_dir, '*')) if os.path.isdir(d)]

        for ntd in neuron_type_dir:

            neuron_type = os.path.basename(os.path.normpath(ntd))

            neuron_dir = [d for d in glob.glob(os.path.join(ntd, '*')) if os.path.isdir(d)]
            neuron_ctr = 0

            for nd in neuron_dir:
                neuron_info = collections.OrderedDict()

                # Find neuron morphology swc file, obs currently assume lowercase(!)
                neuron_morph_list = glob.glob(os.path.join(nd, '*swc'))

                parameter_file = os.path.join(nd, "parameters.json")
                mechanism_file = os.path.join(nd, "mechanisms.json")
                modulation_file = os.path.join(nd, "modulation.json")  # Optional

                if len(neuron_morph_list) == 0:
                    assert (not os.path.isfile(parameter_file) and not os.path.isfile(mechanism_file)), \
                        f"Directory {nd} has parameter.json or mechanism.json but no swc file."

                    # No swc file, skipping directory
                    continue

                # Check if empty neuron_morph_list, or if more than one morphology
                assert len(neuron_morph_list) == 1, f"Should only be one swc file in {nd}"
                assert os.path.isfile(parameter_file), f"Missing parameter file {parameter_file}"
                assert os.path.isfile(mechanism_file), f"Missing mechanism file {mechanism_file}"

                neuron_info["morphology"] = neuron_morph_list[0]
                neuron_info["parameters"] = parameter_file
                neuron_info["mechanisms"] = mechanism_file

                # Modulation file is optional
                if os.path.isfile(modulation_file):
                    neuron_info["modulation"] = modulation_file

                neuron_ctr += 1
                neuron_name = f"{neuron_type}_{neuron_ctr}"

                all_neurons[neuron_name] = neuron_info

            if neuron_ctr > 0:
                print(f"Found {neuron_ctr} neurons in {ntd}")

        return all_neurons


    def has_axon(self, swc_file):
        nm = NeuronMorphology(swc_filename=swc_file)

        return len(nm.axon) > 0


    def create_config(self, num_replicas=10):

        config_def = collections.OrderedDict()

        volume_def = collections.OrderedDict()

        vol_name = "InputTest"

        volume_def[vol_name] = collections.OrderedDict()
        volume_def[vol_name]["type"] = "mesh"
        volume_def[vol_name]["dMin"] = 15e-6
        volume_def[vol_name]["meshFile"] = "data/mesh/InputTestMesh.obj"
        volume_def[vol_name]["meshBinWidth"] = 100e-6

        config_def["Volume"] = volume_def
        config_def["Connectivity"] = dict()  # Unconnected

        neuron_def = self.gather_all_neurons()

        fake_axon_density = ["r", "1", 10e-6]

        for n in neuron_def.keys():
            neuron_def[n]["num"] = num_replicas
            neuron_def[n]["volumeID"] = vol_name
            neuron_def[n]["rotationMode"] = "random"
            neuron_def[n]["hoc"] = None

            if not self.has_axon(neuron_def[n]["morphology"]):
                print(f"Morphology {neuron_def[n]['morphology']} has no axon, faking it.")
                # We will have no connections in this test network, so add empty density
                neuron_def[n]["axonDensity"] = fake_axon_density

        config_def["Neurons"] = neuron_def

        return config_def

    def setup_input(self):

        # TODO: Each replicate of neuron should have its own density
        # This script is a helper script, then we write a master script to do entire range.
        # 1. Generate network, then use same network and neurons for all simulations
        # 2. Separate input generation step, where network is not regenerated.

        # Let the user specify output simulation name
        # Let the user specify range of frequencies, and density


        pass
        # We can vary the input, but not the density during a simulation. To handle varying densities we have
        # the same morphology/model repeated in the simulation


if __name__ == "__main__":
    input_scaling = InputScaling("networks/input_scaling_v1", "data/cellspecs-v2/")

    input_scaling.setup_network()
