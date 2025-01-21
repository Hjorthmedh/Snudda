import os
import ast
from snudda.input.input_tuning import InputTuning

print("Starting setup_verify_input_tuning.py")

# Should be set by script calling setup_input_tuning_dspn
# os.environ["SNUDDA_DATA"] = "../../../../../BasalGangliaData/data/"


assert os.path.isdir(os.getenv("SNUDDA_DATA")), f"You need to have BasalGangliaData installed for this example."

if os.getenv("SNUDDA_TUNE_NEURON"):
    neuron_type = os.getenv("SNUDDA_TUNE_NEURON")
else:
    neuron_type = "dspn"

if os.getenv("INPUT_TYPE"):
    input_type = os.getenv("INPUT_TYPE")
else:
    input_type = "cortical"
        
print(f"Verify input for neuron type {neuron_type}")

network_path = os.path.join("networks", f"verify_input_{neuron_type}_{input_type}")

from ipyparallel import Client

ipython_profile = "default"
ipython_dir = os.getenv('IPYTHONDIR')
if not ipython_dir:
    ipython_dir = os.path.join(os.path.abspath(os.getcwd()), ".ipython")

u_file = os.path.join(ipython_dir, f"profile_{ipython_profile}", "security", "ipcontroller-client.json")
print(f"Reading IPYPARALLEL connection info from {u_file}\n")
# rc = Client(profile=ipython_profile, connection_info=u_file, timeout=120, debug=False)
rc = None

input_tuning = InputTuning(network_path, rc=rc)

print("Constructor done, calling setup_network.")

neurons_path = os.path.join("$DATA", "neurons", "striatum")
input_tuning.setup_network(neurons_path=neurons_path, 
                           num_replicas=25,  # Increased it 5 -> 25 for dspn
                           neuron_types=neuron_type)

print("Calling setup_input")

input_duration = 10
input_frequency_range = [0, 2, 4, 6, 8, 10]

input_tuning.setup_input_verification(input_type=input_type,
                                      neuron_type=neuron_type,
                                      input_frequency_range=input_frequency_range,
                                      input_duration=input_duration,
                                      generate=True, seed_list=None)
                                   

print("All done with setup_input_tuning_background.py")

