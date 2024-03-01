import os
import ast
from snudda.input.input_tuning import InputTuning

print("Starting setup_input_tuning_thalamic_signal.py")

# Should be set by script calling setup_input_tuning_dspn
# os.environ["SNUDDA_DATA"] = "../../../../../BasalGangliaData/data/"


assert os.path.isdir(os.getenv("SNUDDA_DATA")), f"You need to have BasalGangliaData installed for this example."

if os.getenv("SNUDDA_TUNE_NEURON"):
    neuron_type=os.getenv("SNUDDA_TUNE_NEURON")
else:
    neuron_type="dspn"

input_type = "thalamic"

if os.getenv("SEED_LIST"):
    seed_list = ast.literal_eval(os.getenv("SEED_LIST"))
else:
    seed_list = None
    
print(f"Optimising input for neuron type {neuron_type}")

network_path = os.path.join("networks", f"input_tuning_{neuron_type}_{input_type}_signal")


from ipyparallel import Client

ipython_profile = "default"
ipython_dir = os.getenv('IPYTHONDIR')
if not ipython_dir:
    ipython_dir = os.path.join(os.path.abspath(os.getcwd()), ".ipython")

u_file = os.path.join(ipython_dir, f"profile_{ipython_profile}", "security", "ipcontroller-client.json")
print(f"Reading IPYPARALLEL connection info from {u_file}\n")
rc = Client(profile=ipython_profile, connection_info=u_file, timeout=120, debug=False)

input_tuning = InputTuning(network_path, rc=rc, input_seed_list=seed_list)

print("Constructor done, calling setup_network.")

neurons_path = os.path.join("$DATA", "neurons", "striatum")
input_tuning.setup_network(neurons_path=neurons_path, 
                           num_replicas=20,
                           neuron_types=neuron_type)

print("Calling setup_input")

input_tuning.setup_input(input_type=input_type,  # eg. "cortical" or "thalamic"
                         num_input_min=1,
                         num_input_max=200,
                         input_duration=10.0,
                         input_frequency_range=[10.0],
                         use_meta_input=True)


print("All done with setup_input_tuning_cortical_signal.py")

