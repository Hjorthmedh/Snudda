from snudda.simulate.pair_recording import PairRecording
import os

network_name="curinj_network04"
output_name='pair_recording_output'#network-output-6'
network_path = os.path.join("networks",network_name)
#input_config_path=os.path.join("..","..","..","snudda","data","input_config","input-v10-scaled.json")
neurons_dir=os.path.join("..","..","20220930_3","PD0","neurons")
experiment_config_file_norm = os.path.join("..", "..", "..", "snudda", "data", "experiment_config", output_name+'_'+network_name+"_pd2_norm"+"_experiment-config-7-bobek.json")

SNUDDA_DATA_pd2 = os.path.join("..","..","20220930_3","PD2")
network_name_pd2=network_name+ "_pd2"
network_path_pd2 = os.path.join("networks",network_name_pd2)

output_path_pd2_norm = os.path.join(network_path_pd2, "simulation", "output_pd2_norm.hdf5")
pr_pd2_norm = PairRecording(network_path=network_path_pd2, experiment_config_file=experiment_config_file_norm,output_file=output_path_pd2_norm)
pr_pd2_norm.run()