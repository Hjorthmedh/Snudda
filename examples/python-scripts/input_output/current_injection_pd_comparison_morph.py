import os
import json
from collections import OrderedDict
import numpy as np
from snudda.input.input_tuning import InputTuning
from snudda.simulate.pair_recording import PairRecording

from snudda.plotting.plot_input_vs_output_cur import PlotInputVsOutput
import matplotlib.pyplot as plt
from snudda.plotting import PlotTraces

def setting_up_experiment(network_name, parameter_key, morphology_key):

    print('Setting up experiment')
    
    output_name='output_recording'+ '_' + morphology_key + '_' + parameter_key
    network_path = os.path.join("networks",network_name)
    print(f'network path: {network_path}')
    experiment_config_file_template = os.path.join("..", "..", "..", "snudda", "data", "experiment_config", "experiment-config-7-bobek-templ.json")
    experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "experiment_config", output_name+'_'+network_name+"_experiment-config.json")
    print(f'experiment config file: {experiment_config_file}')

    # defining the square pulse current injection protocol
    stim_period = 1         # stimulation period
    no_stim_period = 0.5    # recovery period
    start_time = 0.1        # start time
    n_steps = 30            # number of steps
    min_cur = 2e-10         # min current injected
    max_cur = 6e-10         # max current injected
    end_time = start_time + (stim_period + no_stim_period) * n_steps # end time
    time_step = (end_time - start_time) / n_steps
    
   
    with open(experiment_config_file_template, "r") as f:
        exp_proto = json.load(f, object_pairs_hook=OrderedDict)

    default_inj_proto = exp_proto["currentInjection"][0]        
    num_cell=1    
    for i in np.arange(0, num_cell):
        exp_proto["currentInjection"][i]["neuronID"] = int(i)  # 3
        exp_proto["currentInjection"][i]["start"] = np.arange(start_time, end_time, time_step).tolist()
        exp_proto["currentInjection"][i]["end"] = (np.arange(start_time, end_time, time_step) + stim_period).tolist()
        exp_proto["currentInjection"][i]["amplitude"] = np.linspace(min_cur, max_cur, n_steps).tolist()
    exp_proto["meta"]["simulationDuration"] = end_time
    exp_proto["meta"]["pairRecordingOutputFile"] = output_name + '.hdf5'
    
    with open(experiment_config_file, "w") as f:
        print("Writing to file: " + str(experiment_config_file))
        json.dump(exp_proto, f, indent=2)
    return experiment_config_file, output_name

def setting_up_specific_PD0_cell(network_name, neuron_type):
    network_name_pd0=network_name + "_pd0"
    SNUDDA_DATA_pd0 = os.path.join("..","..","20220930_3","PD0")
    network_path_pd0 = os.path.join("networks",network_name_pd0)
    condition='PD0'
    neurons_path_pd0 = os.path.join(SNUDDA_DATA_pd0, "neurons", "striatum")
    input_tuning_pd0 = InputTuning(network_path_pd0,snudda_data=SNUDDA_DATA_pd0)
    if neuron_type=='ispn':
        single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', neuron_type, neuron_folder)
        input_tuning_pd0.setup_network(neurons_path=neurons_path_pd0,
                           num_replicas=1,
                           neuron_types=neuron_type,
                           single_neuron_path=single_neuron_path,
                           parameter_key=parameter_key_pd0,
                           morphology_key=morphology_key_pd0
                           )
        print(f'####### Selected pd0 cell: {morphology_key_pd0}, {parameter_key_pd0}')
        return network_path_pd0
                           
def setting_up_specific_PD2_cell(network_name, neuron_type):
    network_name_pd2=network_name + "_pd2"
    SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
    network_path_pd2 = os.path.join("networks",network_name_pd2)
    condition='PD2'
    neurons_path_pd2 = os.path.join(SNUDDA_DATA_pd2, "neurons", "striatum")
    input_tuning_pd2 = InputTuning(network_path_pd2,snudda_data=SNUDDA_DATA_pd2)
    
    if neuron_type=='ispn':
        single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', neuron_type, neuron_folder)
        input_tuning_pd2.setup_network(neurons_path=neurons_path_pd2,
                           num_replicas=1,
                           neuron_types=neuron_type,
                           single_neuron_path=single_neuron_path,
                           parameter_key=parameter_key_pd2,
                           morphology_key=morphology_key_pd2
                           )
        print(f'####### Selected pd2 cell: {morphology_key_pd2}, {parameter_key_pd2}')    
    return network_path_pd2
            
def simulate_PD0_cell(network_path_pd0,experiment_config_file_pd0):
    pr = PairRecording(network_path=network_path_pd0,
                      experiment_config_file=experiment_config_file_pd0)
    pr.run()
    '''
    from neuron import h
    for sec in h.allsec():
        h.delete_section(sec=sec)
    import gc
    del pr
    gc.collect()
    '''
    pr.clear_neuron()    
    del pr
    
    print('PD0 simulations: done!')
    
def simulate_PD2_cell(network_path_pd2,experiment_config_file_pd2):
    pr2 = PairRecording(network_path=network_path_pd2, 
                        experiment_config_file=experiment_config_file_pd2)
    pr2.run()
    '''
    from neuron import h
    for sec in h.allsec():
        h.delete_section(sec=sec)
    import gc
    del pr2
    gc.collect()
    '''
    pr2.clear_neuron()    
    del pr2
    
    print('PD2 simulations: done!')    
    
def plot_pd0_sim_out(network_path_pd0, output_name_pd0):    
    network_file = os.path.join(network_path_pd0,"network-synapses.hdf5")
    network_output = os.path.join(network_path_pd0, "simulation", output_name_pd0+'.hdf5')
    pio_pd0 = PlotInputVsOutput(network_path_pd0, network_output, network_file=network_file)
    pio_pd0.plot_cur_input_vs_output(experiment_config_file=experiment_config_file_pd0,
                                     morphology_key= morphology_key_pd0,
                                     parameter_key=parameter_key_pd0)
    pt = PlotTraces(network_output, network_file=network_file)
    pt.plot_traces_sep()
    pt.plot_traces(fig_name= morphology_key_pd0 + '_' + parameter_key_pd0 + '_pd0_traces.pdf')
    return pio_pd0
        
    
def plot_pd2_sim_out(network_path_pd2, output_name_pd2):    
    network_file_pd2 = os.path.join(network_path_pd2,"network-synapses.hdf5")
    network_output_pd2 = os.path.join(network_path_pd2, "simulation", output_name_pd2+'.hdf5')
    pio_pd2 = PlotInputVsOutput(network_path_pd2, network_output_pd2, network_file=network_file_pd2)
    pio_pd2.plot_cur_input_vs_output(experiment_config_file=experiment_config_file_pd2,
                                     morphology_key= morphology_key_pd2,
                                     parameter_key=parameter_key_pd2)
    pt_pd2 = PlotTraces(network_output_pd2, network_file=network_file_pd2)
    #pt_pd2.plot_traces_sep()
    pt_pd2.plot_traces(fig_name= morphology_key_pd2 + '_' + parameter_key_pd2 + '_pd2_traces.pdf')    
    return pio_pd2
    
def plot_pd0_pd2_sim_out():   #plot pd0 and pd2 combined
    cur_steps_pd0, spike_rate_list_pd0 = pio_pd0.get_spike_rate_list(experiment_config_file=experiment_config_file_pd0)
    cur_steps_pd2, spike_rate_list_pd2 = pio_pd2.get_spike_rate_list(experiment_config_file=experiment_config_file_pd2)
    avg_spike_rate_list_pd0 = np.mean(spike_rate_list_pd0,axis=0)
    avg_spike_rate_list_pd2 = np.mean(spike_rate_list_pd2,axis=0)

    fig, ax = plt.subplots()
    ax.plot(cur_steps_pd0, avg_spike_rate_list_pd0, marker='.', label='PD0', color='blue')#,color='black'
    ax.plot(cur_steps_pd2, avg_spike_rate_list_pd2, marker='.', label='PD2', color='red')#,color='black'
    plt.legend()
    ax.set_ylabel('Spikes/s')
    ax.set_xlabel('Current injection (A)')
    ax.set_title(f'{neuron_type} \n pd0: {morphology_key_pd0}, {parameter_key_pd0} \n pd2: {morphology_key_pd2}, {parameter_key_pd2} ')
    figure_path = os.path.join(network_path_pd0, "figures", morphology_key_pd0 + '_' + morphology_key_pd2 + '_' + parameter_key_pd0 + '_' +  "curinj_vs_spikerate_pd0_pd2.png")
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)  
    return cur_steps_pd0, avg_spike_rate_list_pd0, cur_steps_pd2, avg_spike_rate_list_pd2    
    

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Current injection vs output experiment: PD comparisons") 
    parser.add_argument("experiment_name", help="Network name")
    parser.add_argument("neuron_type", help="Neuron type to analyze", choices=['dspn','ispn'])
    
    args = parser.parse_args()

    network_name=args.experiment_name
    neuron_type=args.neuron_type

    ### ispn ###
    #neuron_folder= 'str-ispn-e150908_c4_D2-m51-5-DE'
    neuron_folder= 'str-ispn-e150917_c11_D2-mWT-MSN1'
    #neuron_folder= 'str-ispn-e151123_c1_D2-mWT-P270-09'
    #neuron_folder= 'str-ispn-e160118_c10_D2-m46-3-DE'
    ###
    ### dspn ###
    #neuron_folder= 'str-dspn-e150602_c1_D1-mWT-0728MSN01'
    #neuron_folder= 'str-dspn-e150917_c6_D1-m21-6-DE'
    #neuron_folder= 'str-dspn-e150917_c9_D1-mWT-1215MSN03'
    #neuron_folder= 'str-dspn-e150917_c10_D1-mWT-P270-20'
    ###
    
    # collecting info from combined pd0 and pd2 meta
    meta_path=os.path.join('..', '..', '20220930_3','meta_json_files')
    
    with open(os.path.join(meta_path, neuron_folder+'_merged_meta.json'), 'r') as f:
        meta = json.load(f, object_pairs_hook=OrderedDict)
    
    #'''       
    for par_key in meta.keys():
        parameter_key_pd0=par_key
        parameter_key_pd2=par_key

        morphology_key_pd0_list=[]   
        morphology_key_pd2_list=[]

        for var in meta[par_key].keys():
            morphology_key_pd0_list.append(meta[par_key][var][0])
            morphology_key_pd2_list.append(meta[par_key][var][1])
                  
        morphologies = list(zip(morphology_key_pd0_list,morphology_key_pd2_list))
    
        for i in range(len(morphologies)):
            morphology_key_pd0=morphologies[i][0]
            morphology_key_pd2=morphologies[i][1]

            experiment_config_file_pd0,output_name_pd0=setting_up_experiment(network_name,
                                                                        parameter_key_pd0,
                                                                        morphology_key_pd0)
            experiment_config_file_pd2, output_name_pd2=setting_up_experiment(network_name,
                                                                        parameter_key_pd2,
                                                                        morphology_key_pd2)
            network_path_pd0=setting_up_specific_PD0_cell(network_name, neuron_type)
            network_path_pd2=setting_up_specific_PD2_cell(network_name, neuron_type)
            
            print("### Simulations are starting, remember to compile the mechanisms ###")
            print("### $nrnivmodl ../../20220930_3/mechanisms/ ###")

            simulate_PD0_cell(network_path_pd0, experiment_config_file_pd0)
            simulate_PD2_cell(network_path_pd2, experiment_config_file_pd2)
                                
        print(f'!!!!!!!!!! SIMULATIONS of {par_key} DONE !!!!!!!!!!')
    print('!!!!!!!!!! SIMULATIONS DONE !!!!!!!!!!')
    #'''
    #'''       

    fig_sum, ax_sum = plt.subplots()
    for par_key in meta.keys():
        parameter_key_pd0=par_key
        parameter_key_pd2=par_key

        morphology_key_pd0_list=[]   
        morphology_key_pd2_list=[]

        for var in meta[par_key].keys():
            morphology_key_pd0_list.append(meta[par_key][var][0])
            morphology_key_pd2_list.append(meta[par_key][var][1])
                  
        morphologies = list(zip(morphology_key_pd0_list,morphology_key_pd2_list))
        
        
        for i in range(len(morphologies)):
            morphology_key_pd0=morphologies[i][0]
            morphology_key_pd2=morphologies[i][1]

            experiment_config_file_pd0,output_name_pd0=setting_up_experiment(network_name,
                                                                        parameter_key_pd0,
                                                                        morphology_key_pd0)
            experiment_config_file_pd2, output_name_pd2=setting_up_experiment(network_name,
                                                                        parameter_key_pd2,
                                                                        morphology_key_pd2)
            network_path_pd0=setting_up_specific_PD0_cell(network_name, neuron_type)
            network_path_pd2=setting_up_specific_PD2_cell(network_name, neuron_type)
                           
            pio_pd0=plot_pd0_sim_out(network_path_pd0, output_name_pd0)
            pio_pd2=plot_pd2_sim_out(network_path_pd2, output_name_pd2)
            
            cur_steps_pd0, avg_spike_rate_list_pd0, cur_steps_pd2, avg_spike_rate_list_pd2   = plot_pd0_pd2_sim_out()  
            ax_sum.plot(cur_steps_pd0, avg_spike_rate_list_pd0, marker='.', label='PD0', color='blue', alpha=0.2)
            ax_sum.plot(cur_steps_pd2, avg_spike_rate_list_pd2, marker='.', label='PD2', color='red', alpha=0.2)
    ax_sum.set_ylabel('Spikes/s')
    ax_sum.set_xlabel('Current injection (A)')
    ax_sum.set_title(f'{neuron_type} \n neuron folder: {neuron_folder} ')
    figure_path = os.path.join(network_path_pd0, "figures", neuron_folder + '_' +  "curinj_vs_spikerate_pd0_pd2.png")
    fig_sum.savefig(figure_path, dpi=300)
    plt.close(fig_sum)       
    #'''

    

