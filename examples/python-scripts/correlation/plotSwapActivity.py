import os
import numpy as np
from snudda.plotting import PlotTraces
from snudda.plotting.plot_input_vs_output_corr import PlotInputVsOutput
import matplotlib.pyplot as plt

##############################################
network_name_pd0="testnetwork6_pd0"
network_name_pd2="testnetwork6_pd2"
##############################################

output_name='output'
network_path_pd0 = os.path.join("networks",network_name_pd0)
network_file = os.path.join(network_path_pd0,"network-synapses.hdf5")
network_output = os.path.join(network_path_pd0, "simulation", output_name+'.hdf5')
experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "input_config",
                                      "corr_exp_input.json")

pt = PlotTraces(network_output, network_file=network_file)
pt.plot_traces_sep()
pt.plot_traces(fig_name="traces.pdf")

#Example of actual correlation when really high 0.5-0.99 correlation in the generator was used
#correlation2 = [0.051015826468973596, 0.03691469787465566, 0.026828019362451105, 0.046538034551528204, 0.04777386458376677, 0.022515168296948822, 0.08859144682424376, 0.12762489009798608, 0.11462237035617154, 0.11473182094690294, 0.18053817567478206, 0.2960403803861346, 0.27949973668304595, 0.09650867392460889, 0.37150008781460636, 0.18582919336387593, 0.48256222545329397, 0.7975479645262546, 0.8928468068038496, 0.9985615171893567]

pio = PlotInputVsOutput(network_path_pd0, network_output, network_file=network_file)
cors_pd0, correlation = pio.plot_input_output_corr(experiment_config_file=experiment_config_file)


output_name='output'
network_path_pd2 = os.path.join("networks",network_name_pd2)
network_file = os.path.join(network_path_pd2,"network-synapses.hdf5")
network_output = os.path.join(network_path_pd2, "simulation", output_name+'.hdf5')
experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "input_config",
                                      "corr_exp_input.json")

pt = PlotTraces(network_output, network_file=network_file)
pt.plot_traces_sep()
pt.plot_traces(fig_name="traces.pdf")

pio = PlotInputVsOutput(network_path_pd2, network_output, network_file=network_file)
cors_pd2, correlation = pio.plot_input_output_corr(experiment_config_file=experiment_config_file)


fig, ax = plt.subplots()
ax.plot(correlation,cors_pd0, marker='.', label='PD0')
ax.plot(correlation,cors_pd2, marker='.', label='PD2')
#ax.plot(np.array([correlation,correlation]).transpose(),np.array([cors_pd0,cors_pd2]).transpose(), marker='.')
plt.legend()
ax.set_ylabel('Output correlation')
ax.set_xlabel('Input correlation')
figure_path = os.path.join(network_path_pd0, "figures", "input_output_correlation_pd0_vs_pd2.png")
fig.savefig(figure_path, dpi=300)
plt.close(fig)
#There is a bug, such that that sometimes the last correlation value becomes nan. Not sure why at the moment.
# Some misstake in my plot code I think.

print("end")