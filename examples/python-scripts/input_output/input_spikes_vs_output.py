import os, json
from collections import OrderedDict
import copy
from snudda.input.input_tuning import InputTuning
import numpy as np
import plot_inpspk_outp


def create_networks(ninputs,network_name,numrep,stim_period,no_stim_period,start_time,n_steps,nrn_types):
    subtype_cap=nrn_types['subtype_cap']# = "iSPN"
    para_key=nrn_types['para_key']# = "p3ce323bd"  # "p3ce323bd"
    morph_key_pd0=nrn_types['morph_key_pd0']# = "m688a16a9"
    morph_key_pd2=nrn_types['morph_key_pd2']# = "me4daa6bd"
    subsubtype=nrn_types['subsubtype']# = 'str-ispn-e151123_c1_D2-mWT-P270-09'

    experiment_config_file_template = os.path.join("..", "..", "..", "snudda", "data", "input_config",
                                                   "input_output_" + subtype_cap.lower() + "_template.json")
    experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "input_config", "spk_input_vs_out.json")
    experiment_config_file_high = os.path.join("..", "..", "..", "snudda", "data", "input_config",
                                               "spk_input_vs_out.json")
    # input_config_path=os.path.join("..","..","..","snudda","data","input_config","input-v10-scaled.json")

    ###################################################################
    #Make input file
    ###################################################################
    with open(experiment_config_file_template, "r") as f:
        exp_proto = json.load(f, object_pairs_hook=OrderedDict)

    min_freq = 0
    max_freq = 4
    end_time = start_time + (stim_period + no_stim_period) * n_steps
    time_step = (end_time - start_time) / n_steps
    exp_proto[subtype_cap]["CorticalSignal"]["start"] = np.arange(start_time, end_time, time_step).tolist()
    exp_proto[subtype_cap]["CorticalSignal"]["end"] = (np.arange(start_time, end_time, time_step) + stim_period).tolist()
    #exp_proto[subtype_cap]["CorticalSignal"]["frequency"] = (np.arange(1, 4, 1)).tolist()
    exp_proto[subtype_cap]["CorticalSignal"]["frequency"] = np.linspace(min_freq, max_freq, n_steps).tolist()
    #exp_proto[subtype_cap]["CorticalSignal"]["populationUnitCorrelation"] = np.linspace(min_corr, max_corr, n_steps).tolist()
    exp_proto[subtype_cap]["CorticalSignal"]["nInputs"][subtype_cap] = ninputs
    with open(experiment_config_file, "w") as f:
        print("Writing to file: " + str(experiment_config_file))
        json.dump(exp_proto, f, indent=2)

    exp_proto[subtype_cap]["CorticalSignal"]["conductance"] = 1.4*exp_proto[subtype_cap]["CorticalSignal"]["conductance"]
    exp_proto[subtype_cap]["CorticalSignal"]["nInputs"][subtype_cap] = ninputs
    with open(experiment_config_file_high, "w") as f:
        print("Writing to file: " + str(experiment_config_file_high))
        json.dump(exp_proto, f, indent=2)

    ###################################################################
    #Create PD0 "network"
    ###################################################################
    network_name_pd0=network_name+"_pd0"
    SNUDDA_DATA_pd0 = os.path.join("..","..","20220930_3","PD0")
    network_path_pd0 = os.path.join("networks",network_name_pd0)
    condition='PD0'
    subtype=subtype_cap#.lower()
    #subsubtype = 'str-ispn-e151123_c1_D2-mWT-P270-09'#'str-ispn-e151123_c1_D2-mWT-P270-09'
    single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
    neurons_path_pd0 = os.path.join(SNUDDA_DATA_pd0, "neurons", "striatum")
    input_tuning_pd0 = InputTuning(network_path_pd0,snudda_data=SNUDDA_DATA_pd0)
    input_tuning_pd0.setup_network(neurons_path=neurons_path_pd0,
                               num_replicas=numrep,
                               neuron_types=subtype,
                               single_neuron_path=single_neuron_path,
                               parameter_key=para_key,
                               morphology_key=morph_key_pd0
                               )
    from snudda.input import SnuddaInput
    sii = SnuddaInput(network_path=network_path_pd0,input_config_file=experiment_config_file,verbose=True,random_seed=1234)
    sii.generate()

    # from snudda.plotting import PlotInput
    # input_file = os.path.join(network_path_pd0, "input-spikes.hdf5")
    # spi = PlotInput(input_file)
    # spi.plot_input("iSPN", 10, fig_size=(15,5))


    ###################################################################
    #Create PD2 "network"
    ###################################################################
    SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
    network_name_pd2=network_name+"_pd2"
    network_path_pd2 = os.path.join("networks",network_name_pd2)
    condition='PD2'
    subtype=subtype_cap#.lower()
    #subsubtype = 'str-ispn-e151123_c1_D2-mWT-P270-09'
    single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
    neurons_path_pd2 = os.path.join(SNUDDA_DATA_pd2, "neurons", "striatum")
    input_tuning_pd2 = InputTuning(network_path_pd2,snudda_data=SNUDDA_DATA_pd2)
    input_tuning_pd2.setup_network(neurons_path=neurons_path_pd2,
                               num_replicas=numrep,
                               neuron_types=subtype,
                               single_neuron_path=single_neuron_path,
                               parameter_key=para_key,
                               morphology_key=morph_key_pd2
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


    ###################################################################
    #Create PD2 "network" with strenghtened synapses (increased nInputs)
    ###################################################################
    SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
    network_name_pd2_str=network_name+"_pd2_str"
    network_path_pd2_str = os.path.join("networks",network_name_pd2_str)
    condition='PD2'
    subtype=subtype_cap#.lower()
    #subsubtype = 'str-ispn-e151123_c1_D2-mWT-P270-09'
    single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
    neurons_path_pd2 = os.path.join(SNUDDA_DATA_pd2, "neurons", "striatum")
    input_tuning_pd2 = InputTuning(network_path_pd2_str,snudda_data=SNUDDA_DATA_pd2)
    input_tuning_pd2.setup_network(neurons_path=neurons_path_pd2,
                               num_replicas=numrep,
                               neuron_types=subtype,
                               single_neuron_path=single_neuron_path,
                               parameter_key=para_key,
                               morphology_key=morph_key_pd2
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
    factor=1.4#1.66
    input_spikes = h5py.File(os.path.join(network_path_pd2_str,'input-spikes.hdf5'),'r+')
    nrnids = input_spikes['input'].keys()
    condbefore = input_spikes['input']['0']['CorticalSignal']['conductance'][()]
    for nrnid in nrnids:
        cond = input_spikes['input'][nrnid]['CorticalSignal']['conductance'][()]
        input_spikes['input'][nrnid]['CorticalSignal']['conductance'][()] = cond*factor
    input_spikes.close()

if __name__ == "__main__":
    # Neuron type
    nrn_types = {}
    # nrn_types['subtype_cap'] = "iSPN"
    # nrn_types['para_key'] = "p3ce323bd"  # "p3ce323bd"
    # nrn_types['morph_key_pd0'] = "m688a16a9"
    # nrn_types['morph_key_pd2'] = "me4daa6bd"
    # nrn_types['subsubtype'] = 'str-ispn-e151123_c1_D2-mWT-P270-09'

    # nrn_types['subtype_cap'] = "dSPN"
    # nrn_types['para_key'] = "p09b0945d"  # "p3ce323bd"
    # nrn_types['morph_key_pd0'] = "mbb8e5b24"
    # nrn_types['morph_key_pd2'] = "m0c6d8690"
    # nrn_types['subsubtype'] = 'str-dspn-e150602_c1_D1-mWT-0728MSN01'

    nrn_types['subtype_cap'] = "dSPN"
    nrn_types['para_key'] = "p105d4ed9"
    nrn_types['morph_key_pd0'] = "m93f282f3"
    nrn_types['morph_key_pd2'] = "mf2a73197"
    nrn_types['subsubtype'] = 'str-dspn-e150917_c6_D1-m21-6-DE'

    # Settings
    ninputs = 100
    network_name = "inpspk_vs_output_1.4_"+nrn_types['subtype_cap']+"_c6"
    numrep = 3  # number of neurons to simulate in parallel

    #Time
    stim_period = 1#0.5  # 1
    no_stim_period = 0.5#0.5  # 0.5
    start_time = 0.2#0.2  # 0.2
    n_steps = 5#10
    time=n_steps*(stim_period+no_stim_period)+start_time+0.2

    create_networks(ninputs, network_name, numrep,stim_period,no_stim_period,start_time,n_steps,nrn_types)

    os.system("mpiexec -n 3 snudda simulate networks/" + network_name + "_pd0" + " --time " + str(time))
    os.system("mpiexec -n 3 snudda simulate networks/" + network_name + "_pd2" + " --time " + str(time))
    os.system("mpiexec -n 3 snudda simulate networks/" + network_name + "_pd2_str" + " --time " + str(time))

    plot_inpspk_outp.plot_inpspk_outp_str(network_name,nrn_types['subtype_cap'])