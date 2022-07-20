from check_features_traces import *
import json
import copy
import os
import numpy as np
from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation

def post_analysis(network_path,snudda_data_tools=os.path.join('BasalGangliaData','tools')):
    #Load voltage trace data and validate (pass/fail) against target stated in exp_data.json 
    exp_data = json.load(open(os.path.join(snudda_data_tools,"exp_data.json")))
    mean_resting_membrane = exp_data["background_activity_resting_membrane"]["mean"]
    std_resting_membrane = exp_data["background_activity_resting_membrane"]["std"]
    
    print("Loading input information")
    input_config = json.load(open(os.path.join(network_path, "input_config.json")))
    ls = len(json.load(open(os.path.join(network_path, "tuning-info.json")))["FrequencyRange"])
    
    print("Loading network information")
    network_file = os.path.join(network_path, "network-synapses.hdf5")
    network_info = SnuddaLoad(network_file=network_file)
    trace_id = [x["neuronID"] for x in network_info.data["neurons"]]
    
    print("Loading simulation data")
    simulation_file = os.path.join(network_path,'simulation', "output.hdf5")
    snudda_simulation_load = SnuddaLoadNetworkSimulation(network_simulation_output_file=simulation_file)
    traces = snudda_simulation_load.get_voltage()
    time = snudda_simulation_load.get_time()
    print("Simulation data loaded")
    
    results = dict()
    voltage_pass = dict()
    for r in trace_id:
        res = analyse_tuning(volt=traces[r], num_trials=ls, mean_resting_membrane=mean_resting_membrane,std_resting_membrane=std_resting_membrane)
        idx = np.where(res)[0]
        if len(idx) > 0:
            cell_input_definition = copy.deepcopy(input_config[str(r)])
            for k, inputs in cell_input_definition.items():
                cell_input_definition[k]["frequency"] = np.take(cell_input_definition[k]["frequency"], indices=idx, axis=0)
            list_voltages = split_voltage_trace(traces[r], ls)
            rs = np.take(a=list_voltages, indices=idx, axis=0)
            voltage_pass.update({int(r): {"voltages": rs, "input_config": cell_input_definition}})
        results.update({int(r): res})
    return results, voltage_pass

def post_analysis_spiking(network_path):

    exp_data = json.load(open("exp_data.json"))
    mean_frq = exp_data["frequency"]["mean"]
    std_frq = exp_data["frequency"]["std"]

    input_config = json.load(open(os.path.join(network_path, "input_config.json")))

    input_config = json.load(open(os.path.join(network_path, "input_config.json")))

    print(" Loading ")
    volt_file_name = os.path.join(network_path, "voltage.npy")
    voltage = load_voltage_npy(file_name=volt_file_name)

    cell_ids_name = os.path.join(network_path, "cell_id_simulation.npy")
    cell_ids = load_voltage_npy(file_name=cell_ids_name)[0]
    ls = len(json.load(open(os.path.join(network_path, "tuning-info.json")))["FrequencyRange"])
    print(ls)

    results = dict()
    voltage_pass = dict()
    for i, v in enumerate(voltage):

        if i % 1000 == 0:
            print(f" Running analysis on experiments {i}")

        res = analyse_tuning_spiking(volt=v, num_trials=ls, mean_frq=mean_frq,
                                    std_frq=std_frq)

        idx = np.where(res)[0]
        print(res)
        print(idx)
        if len(idx) > 0:
            c_id = cell_ids[i]
            cell_input_definition = copy.deepcopy(input_config[str(c_id)])
            for k, inputs in cell_input_definition.items():
                cell_input_definition[k]["frequency"] = np.take(cell_input_definition[k]["frequency"], indices=idx, axis=0)
            list_voltages = split_voltage_trace(v, ls)
            rs = np.take(a=list_voltages, indices=idx, axis=0)
            print("adding")
            voltage_pass.update({i: {"voltages": rs, "input_config": cell_input_definition}})

        results.update({i: res})

    return results, voltage_pass
