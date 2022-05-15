import os
import matplotlib.pyplot as plt
import json
from check_features_traces import statistics_spiking

def number_inputs(model, network_path, input_type):
    volt_pass = os.path.join(network_path, "voltages_passed.json")
    with open(volt_pass, "r") as f:
        d = json.load(f)
        
    histogram_data = dict({int(1.0): list()})
    
    for k, r in d.items():
        histogram_data[int(r["input_config"][input_type]["frequency"][0])].append(int(r["input_config"][input_type]['nInputs']))
        
    plt.hist(histogram_data[1], bins=50)
    plt.savefig(os.path.join(network_path,"number_input_hist.png"), dpi=300)
    with open(f"inputs_population_{model}.json", "w") as f:
        json.dump(histogram_data, f)
        
def voltage_values(model, network_path, input_type):
    
    import numpy as np
    histogram_data = dict({int(1.0): list()})
    
    volt_pass = os.path.join(network_path, "voltages_passed.json")
    with open(volt_pass, "r") as f:
        d = json.load(f)

    for k, r in d.items():
        histogram_data[int(r["input_config"][input_type]["frequency"][0])].append(np.mean(r["voltages"]))
        
    plt.hist(histogram_data[1], bins=50)
    plt.savefig(os.path.join(network_path,"voltage_values_hist.png"), dpi=300)


def spike_values(model, network_path, input_type):
    import numpy as np
    histogram_data = dict({int(1.0): list()})

    with open(network_path, "r") as f:
        d = json.load(f)

    for k, r in d.items():
        freq = statistics_spiking(np.array(r["voltages"][0]))
        histogram_data[int(r["input_config"][input_type]["frequency"][0])].append(freq)

    plt.hist(histogram_data[1], bins=50)
    plt.show()