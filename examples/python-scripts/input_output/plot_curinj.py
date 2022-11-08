import os
import numpy as np

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

#plot pd0
from snudda.plotting.plot_input_vs_output_cur import PlotInputVsOutput
import matplotlib.pyplot as plt
from snudda.plotting import PlotTraces

network_file = os.path.join(network_path_pd0,"network-synapses.hdf5")
network_output = os.path.join(network_path_pd0, "simulation", output_name+'.hdf5')
pio = PlotInputVsOutput(network_path_pd0, network_output, network_file=network_file)
pio.plot_cur_input_vs_output(experiment_config_file=experiment_config_file)
pt = PlotTraces(network_output, network_file=network_file)
pt.plot_traces_sep()
pt.plot_traces(fig_name="traces.pdf")

#plot pd2
network_file_pd2 = os.path.join(network_path_pd2,"network-synapses.hdf5")
network_output_pd2 = os.path.join(network_path_pd2, "simulation", output_name+'.hdf5')
pio_pd2 = PlotInputVsOutput(network_path_pd2, network_output_pd2, network_file=network_file_pd2)
pio_pd2.plot_cur_input_vs_output(experiment_config_file=experiment_config_file)
pt_pd2 = PlotTraces(network_output, network_file=network_file_pd2)
pt_pd2.plot_traces_sep()
pt_pd2.plot_traces(fig_name="traces.pdf")

cur_steps_pd0, spike_rate_list_pd0 = pio.get_spike_rate_list(experiment_config_file=experiment_config_file)
cur_steps_pd2, spike_rate_list_pd2 = pio_pd2.get_spike_rate_list(experiment_config_file=experiment_config_file)
avg_spike_rate_list_pd0 = np.mean(spike_rate_list_pd0,axis=0)
avg_spike_rate_list_pd2 = np.mean(spike_rate_list_pd2,axis=0)

#plot pd0 and pd2 combined
fig, ax = plt.subplots()
ax.plot(cur_steps_pd0, avg_spike_rate_list_pd0, marker='.', label='PD0', color='blue')#,color='black'
ax.plot(cur_steps_pd2, avg_spike_rate_list_pd2, marker='.', label='PD2', color='red')#,color='black'
plt.legend()
ax.set_ylabel('Spikes/s')
ax.set_xlabel('Current injection (A)')
figure_path = os.path.join(network_path_pd0, "figures", "curinj_vs_spikerate_pd0_pd2.png")
fig.savefig(figure_path, dpi=300)
plt.close(fig)