import json
import os


def count_passing(model, network_path, input_type):
    neurons = json.load(open(f"neuron_config_{model}.json"))
    volt_pass = os.path.join(network_path, "voltages_passed.json")
    with open(volt_pass, "r") as f:
        d = json.load(f)

    key_volt = dict()
    which_passed = dict()

    for r, data in d.items():

        key_n = "_".join([*neurons[r].values()])

        key_volt.update({key_n: data["voltages"]})
        if key_n in which_passed.keys():
            which_passed[key_n].append(data)
        else:
            which_passed.update({key_n: [data]})

    input_numbers_per_model = dict()
    for k, data in which_passed.items():
        inputs = list()
        for d in data:
            inputs.append(d["input_config"][input_type]["nInputs"])

        input_numbers_per_model.update({k: inputs})

    with open(f"input_map_{model}.json", "w") as f:
        json.dump(input_numbers_per_model, f)
