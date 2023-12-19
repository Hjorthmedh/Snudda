import os
from snudda.input.input_tuning import InputTuning

os.environ["SNUDDA_DATA"] = "../../../../../BasalGangliaData/data/"
assert os.path.isdir(os.getenv("SNUDDA_DATA")), f"You need to have BasalGangliaData installed for this example."

network_path = os.path.join("networks", "input_tuning")
input_tuning = InputTuning(network_path)

neurons_path = os.path.join("$DATA", "neurons", "striatum")
input_tuning.setup_network(neurons_path=neurons_path, 
                           num_replicas=7,
                           neuron_types="dspn")
input_tuning.setup_input(input_type="cortical",  # eg. "cortical" or "thalamic"
                         num_input_min=50,
                         num_input_max=200,
                         input_duration=3.0,
                         input_frequency_range=[2.0, 4.0, 6.0, 8.0, 10.0])

