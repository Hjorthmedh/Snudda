import os,json
from collections import OrderedDict
import numpy as np
import copy
from snudda.input.input_tuning import InputTuning

network_name="curinj_network04"
network_name_pd2_norm=network_name+"_pd2_norm"
output_name='pair_recording_output'#network-output-6'
network_path = os.path.join("networks",network_name)
#input_config_path=os.path.join("..","..","..","snudda","data","input_config","input-v10-scaled.json")
neurons_dir=os.path.join("..","..","20220930_3","PD0","neurons")
experiment_config_file_template = os.path.join("..", "..", "..", "snudda", "data", "experiment_config", "experiment-config-7-bobek-templ.json")
experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "experiment_config", output_name+'_'+network_name+"_experiment-config-7-bobek.json")
experiment_config_file_norm = os.path.join("..", "..", "..", "snudda", "data", "experiment_config", output_name+'_'+network_name_pd2_norm+"_experiment-config-7-bobek.json")
#nrnivmodl C:\Users\bo.bekkouche\PycharmProjects\currentinjection\Snudda\examples\bgd01\parkinson\20211105\PD0\mechanisms
num_reps=1
num_para=3
###############################################################
with open(experiment_config_file_template, "r") as f:
    exp_proto = json.load(f, object_pairs_hook=OrderedDict)
stim_period = 1
no_stim_period = 0.1  # 0.5
start_time = 0.1  # 0.2
n_steps = 2
min_cur = 0.1 * 1e-9
max_cur = 1e-9
end_time = start_time + (stim_period + no_stim_period) * n_steps
time_step = (end_time - start_time) / n_steps
num_subtypes = 1 * num_para * num_reps
default_inj_proto = exp_proto["currentInjection"][0]
exp_proto["currentInjection"] = []
# for i in np.arange(0,num_subtypes):
i = 0
for i in np.arange(0, num_subtypes):
    exp_proto["currentInjection"].append(copy.deepcopy(default_inj_proto))
    exp_proto["currentInjection"][i]["neuronID"] = int(i)  # 3
    exp_proto["currentInjection"][i]["start"] = np.arange(start_time, end_time, time_step).tolist()
    exp_proto["currentInjection"][i]["end"] = (np.arange(start_time, end_time, time_step) + stim_period).tolist()
    exp_proto["currentInjection"][i]["amplitude"] = np.linspace(min_cur, max_cur, n_steps).tolist()
exp_proto["meta"]["simulationDuration"] = end_time
exp_proto["meta"]["pairRecordingOutputFile"] = output_name + '.hdf5'
with open(experiment_config_file, "w") as f:
    print("Writing to file: " + str(experiment_config_file))
    json.dump(exp_proto, f, indent=2)

###############################################################
#Generate current injection experiment files for PD2 with current amplitude normalized to lost dendritic length
###############################################################
with open(experiment_config_file_template, "r") as f:
    exp_proto = json.load(f, object_pairs_hook=OrderedDict)
# stim_period = 0.1#0.5  # 1
# no_stim_period = 0.1#0.5  # 0.5
# start_time = 0.1#0.2  # 0.2
# n_steps = 2
# min_cur = 0.1 * 1e-9
# max_cur = 1e-9
normfac=2085.93/3434.14# around 0.6074
end_time = start_time + (stim_period + no_stim_period) * n_steps
time_step = (end_time - start_time) / n_steps
num_subtypes = 1 * num_para * num_reps
default_inj_proto = exp_proto["currentInjection"][0]
exp_proto["currentInjection"] = []
# for i in np.arange(0,num_subtypes):
i = 0
for i in np.arange(0, num_subtypes):
    exp_proto["currentInjection"].append(copy.deepcopy(default_inj_proto))
    exp_proto["currentInjection"][i]["neuronID"] = int(i)  # 3
    exp_proto["currentInjection"][i]["start"] = np.arange(start_time, end_time, time_step).tolist()
    exp_proto["currentInjection"][i]["end"] = (np.arange(start_time, end_time, time_step) + stim_period).tolist()
    exp_proto["currentInjection"][i]["amplitude"] = (np.linspace(min_cur, max_cur, n_steps)*normfac).tolist()
exp_proto["meta"]["simulationDuration"] = end_time
exp_proto["meta"]["pairRecordingOutputFile"] = output_name + '.hdf5'
with open(experiment_config_file_norm, "w") as f:
    print("Writing to file: " + str(experiment_config_file_norm))
    json.dump(exp_proto, f, indent=2)

###############################################################
###################################################################
#Create PD0 "network"
###################################################################
network_name_pd0=network_name + "_pd0"
SNUDDA_DATA_pd0 = os.path.join("..","..","20220930_3","PD0")
network_path_pd0 = os.path.join("networks",network_name_pd0)
condition='PD0'
subtype='ispn'
subsubtype = 'str-ispn-e151123_c1_D2-mWT-P270-09'
single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
neurons_path_pd0 = os.path.join(SNUDDA_DATA_pd0, "neurons", "striatum")
input_tuning_pd0 = InputTuning(network_path_pd0,snudda_data=SNUDDA_DATA_pd0)
input_tuning_pd0.setup_network(neurons_path=neurons_path_pd0,
                           num_replicas=3,
                           neuron_types=subtype,
                           single_neuron_path=single_neuron_path,
                           parameter_key="p3ce323bd",
                           morphology_key="m688a16a9"
                           )

###################################################################
#Create PD2 "network"
###################################################################
SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
network_name_pd2=network_name + "_pd2"
network_path_pd2 = os.path.join("networks",network_name_pd2)
condition='PD2'
subtype='ispn'
subsubtype = 'str-ispn-e151123_c1_D2-mWT-P270-09'
single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
neurons_path_pd2 = os.path.join(SNUDDA_DATA_pd2, "neurons", "striatum")
input_tuning_pd2 = InputTuning(network_path_pd2,snudda_data=SNUDDA_DATA_pd2)
input_tuning_pd2.setup_network(neurons_path=neurons_path_pd2,
                           num_replicas=3,
                           neuron_types=subtype,
                           single_neuron_path=single_neuron_path,
                           parameter_key="p3ce323bd",
                           morphology_key="me4daa6bd"
                           )

###################################################################
#Swap PD0-PD2
###################################################################
from snudda.utils.swap_to_degenerated_morphologies import SwapToDegeneratedMorphologies
new_network_file = os.path.join(network_path_pd2, "network-synapses.hdf5")
original_network_file = os.path.join(network_path_pd0, "network-synapses.hdf5")
new_input_file = os.path.join(network_path_pd2, "input-spikes.hdf5")
original_input_file = os.path.join(network_path_pd0, "input-spikes.hdf5")
swap = SwapToDegeneratedMorphologies(original_network_file=original_network_file,
                                     new_network_file=new_network_file,
                                     original_snudda_data_dir=SNUDDA_DATA_pd0,
                                     new_snudda_data_dir=SNUDDA_DATA_pd2,
                                     )

from snudda.simulate.pair_recording import PairRecording

# pr = PairRecording(network_path=network_path_pd0, experiment_config_file=experiment_config_file)
# pr.run()

# pr_pd2 = PairRecording(network_path=network_path_pd2, experiment_config_file=experiment_config_file)
# pr_pd2.run()




print("end")
