import os
from snudda.input.input_tuning import InputTuning

# Create a separate dir for Neuromodulation Basal Ganglia Data while tuning
snudda_data = os.getenv("SNUDDA_DATA")

network_path = os.path.join("..", "networks", "dspn_DA_bath")

print(f"Creating network in {network_path}")

input_tuning = InputTuning(network_path, snudda_data=snudda_data)

neurons_path = os.path.join("$DATA", "neurons", "striatum")

input_tuning.setup_network(neurons_path=neurons_path, 
                           num_replicas=1,
                           neuron_types="dspn",
                           reaction_diffusion_file="reaction_diffusion_D1_with_DA_decay.json",
                           network_random_seed=1234)
input_tuning = None


from snudda import Snudda
snd = Snudda(network_path=network_path)
snd.setup_input(input_config="input.json")
