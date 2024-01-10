import os

print("Starting setup_input_tuning_cortical_signal.py")

from snudda.input.input_tuning import InputTuning

print("Import done")

# Should be set by script calling setup_input_tuning_dspn
# os.environ["SNUDDA_DATA"] = "../../../../../BasalGangliaData/data/"


assert os.path.isdir(os.getenv("SNUDDA_DATA")), f"You need to have BasalGangliaData installed for this example."

if os.getenv("SNUDDA_TUNE_NEURON"):
    neuron_type=os.getenv("SNUDDA_TUNE_NEURON")
else:
    neuron_type="dspn"

input_type = "cortical"
    
print(f"Optimising input for neuron type {neuron_type}")

network_path = os.path.join("networks", f"input_tuning_{neuron_type}_{input_type}_signal")
input_tuning = InputTuning(network_path)

print("Constructor done, calling setup_network.")

neurons_path = os.path.join("$DATA", "neurons", "striatum")
input_tuning.setup_network(neurons_path=neurons_path, 
                           num_replicas=41,
                           neuron_types=neuron_type)

print("Calling setup_input")

input_tuning.setup_input(input_type=input_type,  # eg. "cortical" or "thalamic"
                         num_input_min=10,
                         num_input_max=400,
                         input_duration=10.0,
                         input_frequency_range=[10.0],
                         use_meta_input=True)


print("All done with setup_input_tuning.py")

