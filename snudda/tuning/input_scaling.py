import os
import glob
import collections
import json
from snudda.CreateCubeMesh import CreateCubeMesh
from snudda.Neuron_morphology import NeuronMorphology
from snudda.load import SnuddaLoad
from snudda.core import Snudda
import numpy as np


class InputScaling(object):

    def __init__(self, network_dir, cellspec_dir):

        self.network_dir = network_dir
        self.cellspec_dir = cellspec_dir

        self.neuron_types = None
        self.neuron_id = None
        self.input_info = None
        self.neuron_info = None

        if not os.path.isdir(self.network_dir):
            os.mkdir(self.network_dir)

        assert os.path.isdir(self.cellspec_dir), f"Cellspec directory {self.cellspec_dir} does not exist."

        self.network_config_file_name = os.path.join(self.network_dir, "network-config.json")
        self.network_file = os.path.join(self.network_dir, "network-pruned-synapses.hdf5")
        self.input_config_file = os.path.join(self.network_dir, "input-config.json")

        self.core = Snudda(self.network_dir)


    # Writes config files
    def setup_network(self):

        config_def = self.create_network_config(num_replicas=10)

        print(f"Writing network config file to {self.network_config_file_name}")
        with open(self.network_config_file_name, "w") as f:
            json.dump(config_def, f, indent=2)

        CreateCubeMesh("data/mesh/InputTestMesh.obj", [0, 0, 0], 1e-3,
                       description="Mesh file used for Input Scaling")

        from snudda.place import SnuddaPlace
        from snudda.detect import SnuddaDetect
        from snudda.prune import SnuddaPrune

        SnuddaPlace(config_file=self.network_config_file_name)
        SnuddaDetect(config_file=self.network_config_file_name,
                     position_file=os.path.join(self.network_dir,"network-neuron-positions.hdf5"),
                     save_file=os.path.join(self.network_dir, "voxels", "network-putative-synapses.hdf5"))
        SnuddaPrune(work_history_file=os.path.join(self.network_dir,"log","network-detect-worklog.hdf5"))

        # TODO: Also run snudda place, detect, and prune
        # TODO: Make it so that snudda detect and prune are fast if there are no allowed connections
        # TODO: Skip placing neurons that will not receive any inputs or distribute any inputs

    def setup_input(self, input_type=None):

        # TODO: use create_input_config to create config file for input
        # TODO: generate inputs

        synapse_density_cortical_input = "1.15*0.05/(1+np.exp(-(d-30e-6)/5e-6))"
        synapse_density_thalamic_input = "0.05*np.exp(-d/200e-6)"

        if input_type == 'cortical':
            synapse_density = synapse_density_cortical_input
            print("Using cortical synapse density for input.")
        elif input_type == 'thalamic':
            synapse_density = synapse_density_thalamic_input
            print("Using thalamic synapse density for input")
        else:
            synapse_density = "1"
            print("No density profile used for input.")

        self.create_input_config(input_config_file=self.input_config_file,
                                 input_frequency=1.0,
                                 n_input_min=0,
                                 n_input_max=1000,
                                 synapse_conductance=0.5e-9,
                                 synapse_density=synapse_density)


    def run_simulations(self):
        pass
        # TODO: Run simulation

    def analyse_results(self):
        pass
    # TODO: Extract spiking frequency (skip first second for each interval to let network settle)
    # TODO: Create summary graphs

    # This loops through all neuron directories in cellspec in preparation of writing a network config file
    def gather_all_neurons(self):
        all_neurons = collections.OrderedDict()

        neuron_type_dir = [d for d in glob.glob(os.path.join(self.cellspec_dir, '*')) if os.path.isdir(d)]

        self.neuron_types = []

        for ntd in neuron_type_dir:

            neuron_type = os.path.basename(os.path.normpath(ntd))
            self.neuron_types.append(neuron_type)

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

    @staticmethod
    def has_axon(swc_file):
        nm = NeuronMorphology(swc_filename=swc_file)

        return len(nm.axon) > 0

    def create_network_config(self, num_replicas=10):

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

    def load_network_info(self, network_file):

        self.neuron_info = SnuddaLoad(network_file).data["neurons"]

    def collect_neurons(self):

        if not self.neuron_info:
            self.load_network_info(network_file=self.network_file)

        neuron_id = self.neuron_info["neuronID"]
        neuron_type = self.neuron_info["type"]

        neuron_sets = dict()

        for nt in set(neuron_type):
            idx = np.where(neuron_type == nt)
            neuron_id = np.sort(neuron_id[idx])

            neuron_sets[nt] = neuron_id

        return neuron_sets

    # synapse_density and synapse_conductance can be either values (then same for all neurons)
    # or dictionaries, with neuron_type as key.

    def create_input_config(self,
                            input_config_file,
                            n_input_min, n_input_max, input_frequency,
                            synapse_density, synapse_conductance,
                            input_duration=10.0):

        neuron_sets = self.collect_neurons()
        n_inputs = dict()

        for neuron_type in neuron_sets:
            neuron_id_list = neuron_sets[neuron_type]
            num_range = np.linspace(n_input_min, n_input_max, num=len(neuron_id_list)).astype(int)

            for neuron_id, num_input in zip(neuron_id_list, num_range):

                if type(synapse_density) == dict:
                    sd = synapse_density[neuron_type]
                else:
                    sd = synapse_density

                if type(synapse_conductance) == dict:
                    sc = synapse_conductance[neuron_type]
                else:
                    sc = synapse_conductance

                self.add_input(input_target=neuron_id,
                               input_frequency=input_frequency,
                               input_duration=input_duration,
                               input_density=sd,
                               input_conductance=sc)

        with open(input_config_file, "w") as f:
            json.dump(self.input_info,f,indent=4)

    def add_input(self, input_target, input_frequency, input_duration,
                  input_density, input_conductance,
                  synapse_parameter_file):

        if type(input_target) != str:
            input_target = str(input_target)

        if not self.input_info:
            self.input_info = collections.OrderedDict()

        assert type(input_frequency) == list, f"add_input: input_freq should be a list: {input_frequency}"
        n_steps = len(input_frequency)
        input_start = input_duration * np.arange(0, n_steps)
        input_end = input_duration * np.arange(1, n_steps+1)

        self.input_info[input_target] = collections.OrderedDict()

        self.input_info[input_target]["generator"] = "poisson"
        self.input_info[input_target]["type"] = "AMPA_NMDA"
        self.input_info[input_target]["synapseDensity"] = input_density
        self.input_info[input_target]["frequency"] = input_frequency
        self.input_info[input_target]["start"] = input_start
        self.input_info[input_target]["end"] = input_end
        self.input_info[input_target]["populationUnitCorrelation"] = 0.0
        self.input_info[input_target]["jitter"] = 0.0
        self.input_info[input_target]["conductance"] = input_conductance
        self.input_info[input_target]["modFile"] = "tmGlut"
        self.input_info[input_target]["parameterFile"] = synapse_parameter_file





if __name__ == "__main__":
    input_scaling = InputScaling("networks/input_scaling_v1", "data/cellspecs-v2/")

    input_scaling.setup_network()
    input_scaling.setup_input(input_type="thalamic")