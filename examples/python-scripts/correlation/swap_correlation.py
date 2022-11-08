import os, json
from collections import OrderedDict
import copy
from snudda.input.input_tuning import InputTuning
import numpy as np

experiment_config_file_template=os.path.join("..","..","..","snudda","data","input_config","testing_input.json")
experiment_config_file=os.path.join("..","..","..","snudda","data","input_config","corr_exp_input.json")
experiment_config_file_high=os.path.join("..","..","..","snudda","data","input_config","corr_exp_input_high.json")
#input_config_path=os.path.join("..","..","..","snudda","data","input_config","input-v10-scaled.json")

##############################################
#Network (folder) names
network_name_pd0="testnetwork6_pd0"
network_name_pd2="testnetwork6_pd2"
##############################################

###################################################################
#Make input file
###################################################################
with open(experiment_config_file_template, "r") as f:
    exp_proto = json.load(f, object_pairs_hook=OrderedDict)
stim_period = 0.1#2.0  # 1
no_stim_period = 0.1#0.2  # 0.5
start_time = 0.1#0.6  # 0.2
n_steps = 3#20
min_corr = 0.5
max_corr = 0.999
end_time = start_time + (stim_period + no_stim_period) * n_steps
time_step = (end_time - start_time) / n_steps
exp_proto['iSPN']["CorticalSignal"]["start"] = np.arange(start_time, end_time, time_step).tolist()
exp_proto['iSPN']["CorticalSignal"]["end"] = (np.arange(start_time, end_time, time_step) + stim_period).tolist()
exp_proto['iSPN']["CorticalSignal"]["populationUnitCorrelation"] = np.linspace(min_corr, max_corr, n_steps).tolist()
exp_proto['iSPN']["CorticalSignal"]["nInputs"]["iSPN"] = 100
with open(experiment_config_file, "w") as f:
    print("Writing to file: " + str(experiment_config_file))
    json.dump(exp_proto, f, indent=2)

exp_proto['iSPN']["CorticalSignal"]["nInputs"]["iSPN"] = 100
with open(experiment_config_file_high, "w") as f:
    print("Writing to file: " + str(experiment_config_file_high))
    json.dump(exp_proto, f, indent=2)

###################################################################
#Create PD0 "network"
###################################################################
SNUDDA_DATA_pd0 = os.path.join("..","..","20220930_3","PD0")
network_path_pd0 = os.path.join("networks",network_name_pd0)
condition='PD0'
subtype='iSPN'
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
from snudda.input import SnuddaInput
sii = SnuddaInput(network_path=network_path_pd0,input_config_file=experiment_config_file,verbose=True,random_seed=1234)
sii.generate()

from snudda.plotting import PlotInput
input_file = os.path.join(network_path_pd0, "input-spikes.hdf5")
spi = PlotInput(input_file)
spi.plot_input("iSPN", 10, fig_size=(15,5))
print("stop")


###################################################################
#Create PD2 "network"
###################################################################
SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
network_path_pd2 = os.path.join("networks",network_name_pd2)
condition='PD2'
subtype='iSPN'
subsubtype = 'str-ispn-e151123_c1_D2-mWT-P270-09'
single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
neurons_path_pd2 = os.path.join(SNUDDA_DATA_pd2, "neurons", "striatum")
input_tuning_pd2 = InputTuning(network_path_pd2,snudda_data=SNUDDA_DATA_pd2)
input_tuning_pd2.setup_network(neurons_path=neurons_path_pd2,
                           num_replicas=10,
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
                                     original_input_file=original_input_file,
                                     new_input_file=new_input_file
                                     )
swap.write_new_network_file()
swap.write_new_input_file()
swap.close()

if 0:#Pausing creation of strenghted network for now
    ###################################################################
    #Create PD2 "network" with strenghtened synapses (increased nInputs)
    ###################################################################
    SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
    network_name_pd2_str="testnetwork3_pd2_strength"
    network_path_pd2_str = os.path.join("networks",network_name_pd2_str)
    condition='PD2'
    subtype='iSPN'
    subsubtype = 'str-ispn-e151123_c1_D2-mWT-P270-09'
    single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
    neurons_path_pd2 = os.path.join(SNUDDA_DATA_pd2, "neurons", "striatum")
    input_tuning_pd2 = InputTuning(network_path_pd2_str,snudda_data=SNUDDA_DATA_pd2)
    input_tuning_pd2.setup_network(neurons_path=neurons_path_pd2,
                               num_replicas=10,
                               neuron_types=subtype,
                               single_neuron_path=single_neuron_path,
                               parameter_key="p3ce323bd",
                               morphology_key="me4daa6bd"
                               )

    ###################################################################
    #Swap PD2_strength with PD0 so that it gets same input synapses as PD2
    ###################################################################
    from snudda.utils.swap_to_degenerated_morphologies import SwapToDegeneratedMorphologies
    new_network_file = os.path.join(network_path_pd2_str, "network-synapses.hdf5")
    original_network_file = os.path.join(network_path_pd0, "network-synapses.hdf5")
    new_input_file = os.path.join(network_path_pd2_str, "input-spikes.hdf5")
    original_input_file = os.path.join(network_path_pd0, "input-spikes.hdf5")
    swap = SwapToDegeneratedMorphologies(original_network_file=original_network_file,
                                         new_network_file=new_network_file,
                                         original_snudda_data_dir=SNUDDA_DATA_pd0,
                                         new_snudda_data_dir=SNUDDA_DATA_pd2,
                                         original_input_file=original_input_file,
                                         new_input_file=new_input_file
                                         )

    swap.write_new_network_file()
    swap.write_new_input_file()
    swap.close()

    ###################################################################
    #Modify synapses of PD2 strengthed
    ###################################################################
    import h5py
    factor=2
    input_spikes = h5py.File(os.path.join(network_path_pd2_str,'input-spikes.hdf5'),'r+')
    nrnids = input_spikes['input'].keys()
    condbefore = input_spikes['input']['0']['CorticalSignal']['conductance'][()]
    for nrnid in nrnids:
        cond = input_spikes['input'][nrnid]['CorticalSignal']['conductance'][()]
        input_spikes['input'][nrnid]['CorticalSignal']['conductance'][()] = cond*factor
    input_spikes.close()

    ###################################################################
    #Create PD2 "network" with rewired synapses (increased density)
    ###################################################################
    SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
    network_name_pd2_rew="testnetwork3_pd2_rewired"
    network_path_pd2_rew = os.path.join("networks",network_name_pd2)
    condition='PD2'
    subtype='iSPN'
    subsubtype = 'str-ispn-e151123_c1_D2-mWT-P270-09'
    single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
    neurons_path_pd2 = os.path.join(SNUDDA_DATA_pd2, "neurons", "striatum")
    input_tuning_pd2 = InputTuning(network_path_pd2_rew,snudda_data=SNUDDA_DATA_pd2)
    input_tuning_pd2.setup_network(neurons_path=neurons_path_pd2,
                               num_replicas=10,
                               neuron_types=subtype,
                               single_neuron_path=single_neuron_path,
                               parameter_key="p3ce323bd",
                               morphology_key="me4daa6bd")

print("end")