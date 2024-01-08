import os

print("Starting setup_input_tuning.py")

from snudda.input.input_tuning import InputTuning

print("Import done")

# Should be set by script calling setup_input_tuning_dspn
# os.environ["SNUDDA_DATA"] = "../../../../../BasalGangliaData/data/"


assert os.path.isdir(os.getenv("SNUDDA_DATA")), f"You need to have BasalGangliaData installed for this example."

if os.getenv("SNUDDA_TUNE_NEURON"):
    neuron_type=os.getenv("SNUDDA_TUNE_NEURON")
else:
    neuron_type="dspn"

print(f"Optimising input for neuron type {neuron_type}")

network_path = os.path.join("networks", f"input_tuning_{neuron_type}")
input_tuning = InputTuning(network_path)

print("Constructor done, calling setup_network.")

neurons_path = os.path.join("$DATA", "neurons", "striatum")
input_tuning.setup_network(neurons_path=neurons_path, 
                           num_replicas=21,
                           neuron_types=neuron_type)

print("Calling setup_input")


input_tuning.setup_background_input(input_types=["cortical_background", "thalamic_background"],
                                    input_density=["1.15*0.05/(1+exp(-(d-30e-6)/5e-6))", "0.05*exp(-d/200e-6)"],
                                    input_fraction=[0.5, 0.5],
                                    num_input_min=10, num_input_max=500,
                                    input_frequency=[1, 1], input_duration=10)

print("All done with setup_input_tuning.py")

