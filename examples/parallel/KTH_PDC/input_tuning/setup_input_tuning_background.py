import os
import ast

from snudda.input.input_tuning import InputTuning

print("Starting setup_input_tuning_background.py")

# Should be set by script calling setup_input_tuning_dspn
# os.environ["SNUDDA_DATA"] = "../../../../../BasalGangliaData/data/"


assert os.path.isdir(os.getenv("SNUDDA_DATA")), f"You need to have BasalGangliaData installed for this example."

if os.getenv("SNUDDA_TUNE_NEURON"):
    neuron_type = os.getenv("SNUDDA_TUNE_NEURON")
else:
    neuron_type = "dspn"

if os.getenv("SEED_LIST"):
    seed_list = ast.literal_eval(os.getenv("SEED_LIST"))
else:
    seed_list = None
    
print(f"Optimising input for neuron type {neuron_type}")

network_path = os.path.join("networks", f"input_tuning_{neuron_type}_background")

from ipyparallel import Client

try:
    ipython_profile = "default"
    ipython_dir = os.getenv('IPYTHONDIR')
    if not ipython_dir:
        ipython_dir = os.path.join(os.path.abspath(os.getcwd()), ".ipython")

    u_file = os.path.join(ipython_dir, f"profile_{ipython_profile}", "security", "ipcontroller-client.json")
    print(f"Reading IPYPARALLEL connection info from {u_file}\n")
    rc = Client(profile=ipython_profile, connection_info=u_file, timeout=120, debug=False)
except:
    print("\n!!! --> Failed to start ipyparallel\n")
    rc = None
    
input_tuning = InputTuning(network_path, rc=rc, input_seed_list=seed_list)

print("Constructor done, calling setup_network.")

neurons_path = os.path.join("$DATA", "neurons", "striatum")
input_tuning.setup_network(neurons_path=neurons_path, 
                           num_replicas=31,
                           neuron_types=neuron_type)

print("Calling setup_input")

if neuron_type in ["lts"]:
    input_types = ["cortical_background"]
    input_density = ["1.15*0.05/(1+exp(-(d-30e-6)/5e-6))"]
    input_fraction = [1]
    input_frequency = [2]
else:
    input_types = ["cortical_background", "thalamic_background"]
    input_density = ["1.15*0.05/(1+exp(-(d-30e-6)/5e-6))", "0.05*exp(-d/200e-6)"]
    input_fraction = [0.5, 0.5]
    input_frequency = [2, 2]


input_tuning.setup_background_input(input_types=input_types,
                                    input_density=input_density,
                                    input_fraction=input_fraction,
                                    num_input_min=10, num_input_max=1000,
                                    input_frequency=input_frequency,
                                    input_duration=10)

print("All done with setup_input_tuning_background.py")

