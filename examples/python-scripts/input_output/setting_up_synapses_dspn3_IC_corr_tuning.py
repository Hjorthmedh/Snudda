import os
import json
from collections import OrderedDict
from snudda.input.input_tuning import InputTuning
import numpy as np

experiment_config_file_template=os.path.join("..","..","..","snudda","data","input_config","input_output_dspn_template_IC.json")

PD0_mkeys=["m615b0265", "mb5fabc3d", "maa9c90e2", "m5089d2b2", "me0f4fb29", "maabff8a9", "m58455f0c",
"m33883bbb", "m5d4aac41"]

PD0_pkeys=["p03c5181a", "p0ffc96fb", "p53c1927c", "pa0f43cf5", "pac70f661", "pd6f5bb48"]


#correlations=[0, 0.25, 0.5]
correlations=[0]

#number_synapses=[210, 220, 235, 240, 225, 235] # "m615b0265" 
#freq           =[6,   1  , 6  ,   6,   15,   5] 
#number_synapses=[215, 230, 240, 245, 205,  240] 
#freq           =[8,   3  , 7  ,   7,   4,   6] 
#number_synapses=[220, 250, 245, 255, 215,  240] 
#freq           =[9,   10  , 8  ,   9,   8,   6] 
#number_synapses=[225, 250, 255, 265, 220,  255]
#freq           =[15,   10,  11, 12.5,  10,   9] 
#number_synapses=[222, 250, 250, 260, 220,  260] # done
#####################################################################################################
#number_synapses=[160, 230, 230, 260, 210, 260] # "mb5fabc3d"
#freq           =[0,   1  , 2  ,   8,   4,   6] 
#number_synapses=[170, 240, 240, 265, 220, 265] 
#freq           =[0,   2  , 4  ,   9,   7,   8] 
#number_synapses=[230, 260, 250, 270, 230, 270] 
#freq           =[10,   8  , 8  , 10,   10,   9] 
#number_synapses=[230, 270, 260, 270, 230, 275] 
#freq           =[10+,  12,  10,  10,   11,   10+] 
#number_synapses=[230, 265, 260, 270, 230, 275] # #done

#####################################################################################################

#number_synapses=[220, 220, 220, 220, 220, 220] # "maa9c90e2"
#freq          =[12 ,   4,   4,   4,   12,   5] 
#number_synapses=[210, 250, 250, 250, 210, 240]
#freq           =[8  ,  20,   18,12.5, 10,  10] 
#number_synapses=[215, 235, 235, 240, 210, 240] #done


#####################################################################################################
#number_synapses=[250, 250, 250, 250, 250, 250] # "m5089d2b2"
#freq           =[20  ,  2,   4,   3, 20,  2] 
#number_synapses=[240, 270, 270, 270, 240, 270]
#freq           =[15  , 5,   8,   6, 15,  5] 
#number_synapses=[225, 285, 280, 290, 230, 285]
#freq           =[ 9  , 10,   10,   10, 9,  7] 
#number_synapses=[227, 285, 280, 290, 235, 295]
#freq           =[ 6  , 10,   10,   10, 10,  10] 
#####################################################################################################

#number_synapses=[220, 220, 220, 220, 220, 220] # "me0f4fb29"
#freq           =[15 ,  7,   7,   6, 15,  7] 
#number_synapses=[205, 230, 230, 230, 205, 235]
#freq           =[8 ,  12,   10,   8, 9,  12.5]
#number_synapses=[210, 225, 230, 235, 205, 230]
#freq           =[10 ,  12,   10,   9, 9,  10]


from snudda.input import SnuddaInput

#for i in range(len(PD0_mkeys)):
for i in [4]:
    morphology_key_pd0=PD0_mkeys[i]

    for j in range(len(PD0_pkeys)):
    #for j in range(1):
        parameter_key_pd0=PD0_pkeys[j]
        
        for corr in correlations:

            network_name_pd0="dspn3_corr_"+str(corr)+'_tuning_'+str(morphology_key_pd0)+"_"+str(parameter_key_pd0)+"_pd0"
            experiment_config_file=os.path.join("..","..","..","snudda","data","input_config","dspn3_corr_"+str(corr)+'_tuning_'+str(morphology_key_pd0)+"_"+str(parameter_key_pd0)+".json")

            ###################################################################
            #Make input file
            ###################################################################
            with open(experiment_config_file_template, "r") as f:
                exp_proto = json.load(f, object_pairs_hook=OrderedDict)
            stim_period = 2
            no_stim_period = 0.5
            start_time = 0.5
            n_steps = 10
            min_freq = 0.5
            max_freq = 5
            end_time = start_time + (stim_period + no_stim_period) * n_steps
            print(f'################################# end time: {end_time}')
            time_step = (end_time - start_time) / n_steps

            exp_proto['dspn']["CorticalSignal"]["start"] = np.arange(start_time, end_time, time_step).tolist()
            exp_proto['dspn']["CorticalSignal"]["end"] = (np.arange(start_time, end_time, time_step) + stim_period).tolist()
            exp_proto['dspn']["CorticalSignal"]["frequency"] = np.linspace(min_freq, max_freq, n_steps).tolist()
            exp_proto['dspn']["CorticalSignal"]["nInputs"]["dspn"] = number_synapses[j]
            exp_proto['dspn']["CorticalSignal"]["populationUnitCorrelation"] = corr
            with open(experiment_config_file, "w") as f:
                print("Writing to file: " + str(experiment_config_file))
                json.dump(exp_proto, f, indent=2)


            SNUDDA_DATA_pd0 = os.path.join("..","..","20220930_3","PD0")
            network_path_pd0 = os.path.join("networks","dspn3",network_name_pd0)
            
            
            ###################################################################
            #Create PD0 "network"
            ###################################################################

            SNUDDA_DATA_pd0 = os.path.join("..","..","20220930_3","PD0")
            network_path_pd0 = os.path.join("networks","dspn3",network_name_pd0)
            condition='PD0'
            subtype='dspn'
            subsubtype = 'str-dspn-e150917_c9_D1-mWT-1215MSN03'
            single_neuron_path = os.path.join('..', '..', '20220930_3', condition, 'neurons', 'striatum', subtype, subsubtype)
            neurons_path_pd0 = os.path.join(SNUDDA_DATA_pd0, "neurons", "striatum")
            input_tuning_pd0 = InputTuning(network_path_pd0,snudda_data=SNUDDA_DATA_pd0)
            input_tuning_pd0.setup_network(neurons_path=neurons_path_pd0,
                                       num_replicas=10,
                                       neuron_types=subtype,
                                       single_neuron_path=single_neuron_path,
                                       parameter_key=parameter_key_pd0,
                                       morphology_key=morphology_key_pd0
                                       )
            from snudda.input import SnuddaInput
            sii = SnuddaInput(network_path=network_path_pd0,input_config_file=experiment_config_file,verbose=True,random_seed=1234)
            sii.generate()

