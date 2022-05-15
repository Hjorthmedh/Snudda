from .check_features_traces import *
import json
import os


def plot_all_frequency(network_path, input_type, verbose=False):
    volt_file_name = os.path.join(network_path, "voltage.npy")
    voltage = load_voltage_npy(file_name=volt_file_name)

    cell_ids_name = os.path.join(network_path, "cell_id_simulation.npy")
    cell_ids = load_voltage_npy(file_name=cell_ids_name)[0]

    ls = len(json.load(open(os.path.join(network_path, "tuning-info.json")))["FrequencyRange"])

    print(ls)
    exp_data = json.load(open("exp_data.json"))
    mean_frq = exp_data["frequency"]["mean"]
    std_frq = exp_data["frequency"]["std"]

    infos = json.load(open(os.path.join(network_path, "input_config.json")))

    if isinstance(input_type, str):
        input_type = [input_type]

    for i, vs in enumerate(voltage):
        info = infos[str(cell_ids[i])]

        for k, v in enumerate(split_voltage_trace(vs, ls)):
            print(" The input config \n")
            for input_t in input_type:
                if verbose:
                    print(f" Frequency {info[input_t]['frequency'][k]} \n")
                    print(f" Number input {info[input_t]['nInputs']} \n")
                frq = get_frequency(v)
                compare_normal_distribution(frq, mean_frq, std_frq,
                                            "_".join([f"cell_id_{i}", f"frequency_{info[input_t]['frequency'][k]},"
                                                     f"input_{info[input_t]['nInputs']}"]))


def plot_all_voltages(network_path, input_type):

    ls = len(json.load(open(os.path.join(network_path, "tuning-info.json")))["FrequencyRange"])

    exp_data = json.load(open("exp_data.json"))
    mean_resting_membrane = exp_data["background_activity_resting_membrane"]["mean"]
    std_resting_membrane = exp_data["background_activity_resting_membrane"]["std"]

    infos = json.load(open(os.path.join(network_path, "input_config.json")))

    volt_file_name = os.path.join(network_path, "voltage.npy")
    voltage = load_voltage_npy(file_name=volt_file_name)

    cell_ids_name = os.path.join(network_path, "cell_id_simulation.npy")
    cell_ids = load_voltage_npy(file_name=cell_ids_name)[0]

    if isinstance(input_type, str):
        input_type = [input_type]

    for i, vs in enumerate(voltage):
        info = infos[str(cell_ids[i])]

        for k, v in enumerate(split_voltage_trace(vs, ls)):
            print(" The input config \n")
            for input_t in input_type:
                print(f" Frequency {info[input_t]['frequency'][k]} \n")
                print(f" Number input {info[input_t]['nInputs']} \n")
            compare_normal_distribution(v, mean_resting_membrane, std_resting_membrane,
                                        "_".join([f"all_cell_id_{i}", f"frequency_{info[input_t]['frequency'][k]},"
                                                 f"input_{info[input_t]['nInputs']}"]))


def plot_passing_voltage_distributions(network_path, input_type):
    exp_data = json.load(open("exp_data.json"))
    mean_resting_membrane = exp_data["background_activity_resting_membrane"]["mean"]
    std_resting_membrane = exp_data["background_activity_resting_membrane"]["std"]

    voltages_passed = json.load(open(os.path.join(network_path, "voltages_passed.json")))

    if isinstance(input_type, str):
        input_type = [input_type]

    for i, info in voltages_passed.items():

        for k, v in enumerate(info["voltages"]):

            print(" The input config \n")
            for input_t in input_type:
                print(f" Frequency {info['input_config'][input_t]['frequency'][k]} \n")
                print(f" Number input {info['input_config'][input_t]['nInputs']} \n")

            compare_normal_distribution(v, mean_resting_membrane, std_resting_membrane, title_name=f" cell_{i}_frequency_{info['input_config'][input_t]['frequency'][k]}_"
                                                                                                   f"number_input_{info['input_config'][input_t]['nInputs']} \n")


def plot_passing_frequency_distributions(network_path, input_type):
    exp_data = json.load(open("exp_data.json"))
    mean_resting_membrane = exp_data["frequency"]["mean"]
    std_resting_membrane = exp_data["frequency"]["std"]

    voltages_passed = json.load(open(os.path.join(network_path, "voltages_passed.json")))

    if isinstance(input_type, str):
        input_type = [input_type]

    for i, info in voltages_passed.items():

        for k, v in enumerate(info["voltages"]):
            print(" The input config \n")
            for input_t in input_type:
                print(f" Frequency {info['input_config'][input_t]['frequency'][k]} \n")
                print(f" Number input {info['input_config'][input_t]['nInputs']} \n")

            compare_normal_distribution(v, mean_resting_membrane, std_resting_membrane)
