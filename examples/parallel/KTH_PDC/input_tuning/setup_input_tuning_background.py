import os
import ast
import json
from snudda.input.input_tuning import InputTuning
from snudda.utils.snudda_path import snudda_isdir, snudda_parse_path, snudda_simplify_path, get_snudda_data

print("Starting setup_input_tuning_background.py")

# Should be set by script calling setup_input_tuning_dspn
# os.environ["SNUDDA_DATA"] = "../../../../../BasalGangliaData/data/"


assert os.path.isdir(os.getenv("SNUDDA_DATA")), f"You need to have BasalGangliaData installed for this example."

if os.getenv("SNUDDA_TUNE_NEURON"):
    neuron_type = os.getenv("SNUDDA_TUNE_NEURON")
else:
    neuron_type = "dspn"

# https://www.jneurosci.org/content/27/16/4374 -- Lacey, Bolam, Magill 2007
input_fraction_lookup = {"dspn": [0.62, 0.38],
                         "ispn": [0.74, 0.26]}
    
input_fraction = input_fraction_lookup.get(neuron_type, [0.5, 0.5])
    
if os.getenv("SEED_LIST"):
    seed_list = ast.literal_eval(os.getenv("SEED_LIST"))
else:
    seed_list = None
    
print(f"Optimising input for neuron type {neuron_type}")

network_path = os.path.join("networks", f"input_tuning_{neuron_type}_background")

from ipyparallel import Client

if False:
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
else:
    print(f"Disabled PARALLEL")
    rc = None
        
input_tuning = InputTuning(network_path, rc=rc, input_seed_list=seed_list)

print("Constructor done, calling setup_network.")

neurons_path = os.path.join("$DATA", "neurons", "striatum")
input_tuning.setup_network(neurons_path=neurons_path, 
                           num_replicas=31,
                           neuron_types=neuron_type)

print("Calling setup_input")

input_info_file = snudda_parse_path(os.path.join("$SNUDDA_DATA", "input_config", "input_info.json"),
                                    os.getenv("SNUDDA_DATA"))

with open(input_info_file, "r") as f:
    extra_info = json.load(f)

if neuron_type.lower() == "dspn":
    extra_input = {"GPe_FoxP2": extra_info["dspn"]["GPe_FoxP2"] }
elif neuron_type.lower() == "ispn":
    extra_input = {"GPe_FoxP2": extra_info["ispn"]["GPe_FoxP2"] }
elif neuron_type.lower() == "fs":
    extra_input = {"GPe_PV": extra_info["fs"]["GPe_PV"] }
else:
    extra_input = None

input_tuning.setup_background_input(input_types=["cortical_background", "thalamic_background"],
                                    input_density=["1.15*0.05/(1+exp(-(d-30e-6)/5e-6))", "1"],
                                    input_fraction=input_fraction,
                                    num_input_min=10, num_input_max=1000,
                                    input_frequency=[2, 2],
                                    input_duration=10,
                                    extra_fixed_input=extra_input)

print("All done with setup_input_tuning_background.py")

