import os
import numpy as np
from snudda.plotting import PlotTraces
from snudda.plotting.plot_input_vs_output import PlotInputVsOutput
import matplotlib.pyplot as plt
network_name="inpspk_vs_output_03_pd0"
output_name='output'
network_path_pd0 = os.path.join("networks",network_name)
network_file = os.path.join(network_path_pd0,"network-synapses.hdf5")
network_output = os.path.join(network_path_pd0, "simulation", output_name+'.hdf5')
experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "input_config",
                                      "spk_input_vs_out.json")

pt = PlotTraces(network_output, network_file=network_file)
pt.plot_traces_sep()
pt.plot_traces(fig_name="traces.pdf")
pio = PlotInputVsOutput(network_path_pd0, network_output, network_file=network_file)
freqs_pd0, out_freq_pd0 = pio.get_input_vs_output(experiment_config_file=experiment_config_file)

network_name="inpspk_vs_output_03_pd2"
output_name='output'
network_path_pd2 = os.path.join("networks",network_name)
network_file = os.path.join(network_path_pd2,"network-synapses.hdf5")
network_output = os.path.join(network_path_pd2, "simulation", output_name+'.hdf5')
experiment_config_file = os.path.join("..", "..", "..", "snudda", "data", "input_config",
                                      "spk_input_vs_out.json")

pt = PlotTraces(network_output, network_file=network_file)
pt.plot_traces_sep()
pt.plot_traces(fig_name="traces.pdf")

pio_pd2 = PlotInputVsOutput(network_path_pd2, network_output, network_file=network_file)
freqs_pd2, out_freq_pd2 = pio_pd2.get_input_vs_output(experiment_config_file=experiment_config_file)


fig, ax = plt.subplots()
ax.plot(freqs_pd0, out_freq_pd0, marker='.', label='PD0', color='blue', alpha=0.2)
ax.plot(freqs_pd0, np.mean(out_freq_pd0, axis=1), marker='.', label='PD0', color='blue', linewidth=2)
ax.plot(freqs_pd2, out_freq_pd2, marker='.', label='PD2', color='red', alpha=0.2)
ax.plot(freqs_pd2, np.mean(out_freq_pd2, axis=1), marker='.', label='PD2', color='red', linewidth=2)

#ax.plot(np.array([correlation,correlation]).transpose(),np.array([cors_pd0,cors_pd2]).transpose(), marker='.')
plt.legend()
ax.set_ylabel('Output frequency (Hz)')
ax.set_xlabel('Input frequency (Hz)')
figure_path = os.path.join(network_path_pd0, "figures", "input_output_freq_pd0_pd2.png")
fig.savefig(figure_path, dpi=300)
plt.close(fig)

print("end")