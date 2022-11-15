from snudda.simulate.pair_recording import PairRecording
import os

network_name="curinj_network03"
output_name='pair_recording_output'#network-output-6'
network_path = os.path.join("networks",network_name)
#input_config_path=os.path.join("..","..","..","snudda","data","input_config","input-v10-scaled.json")
neurons_dir=os.path.join("..","..","20220930_3","PD0","neurons")
experiment_config_file_template = os.path.join("..", "..", "..", "snudda", "data", "experiment_config", "experiment-config-7-bobek-templ.json")
experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "experiment_config", output_name+'_'+network_name+"_experiment-config-7-bobek.json")

network_name_pd0=network_name + "_pd0"
SNUDDA_DATA_pd0 = os.path.join("..","..","20220930_3","PD0")
network_path_pd0 = os.path.join("networks",network_name_pd0)

SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
network_name_pd2=network_name + "_pd2"
network_path_pd2 = os.path.join("networks",network_name_pd2)

pr = PairRecording(network_path=network_path_pd0, experiment_config_file=experiment_config_file)
pr.run()

#For some reason it did not work to run directly after
# pr_pd2 = PairRecording(network_path=network_path_pd2, experiment_config_file=experiment_config_file)
# pr_pd2.run()