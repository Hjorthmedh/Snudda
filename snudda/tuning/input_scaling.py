import os
import glob
import collections


class InputScaling(object):

    def __init__(self, network_dir, cellspec_dir):

        self.network_dir = network_dir
        self.cellspec_dir = cellspec_dir

        if not os.path.isdir(self.network_dir):
            os.mkdir(self.network_dir)

        assert os.path.isdir(self.cellspec_dir), f"Cellspec directory {self.cellspec_dir} does not exist."

        all_neurons = self.gather_all_neurons()

        print(all_neurons)


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
                modulation_file = os.path.join(nd, "modulation.json")   # Optional

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

        return all_neurons





if __name__ == "__main__":
    input_scaling = InputScaling("networks/inp_scaling", "data/cellspecs-v2/")

