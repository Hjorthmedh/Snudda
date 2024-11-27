import os
import numpy as np
from snudda.plotting.plot_input_locations import SnuddaPlotInputLocations


class PlotDegeneration:

    def __init__(self, original_network_path, degenerated_network_path,
                 original_input=None, degenerated_input=None,
                 original_snudda_data=None, degenerated_snudda_data=None):

        if os.path.isfile(original_network_path):
            original_network_file = original_network_path
            original_network_path = os.path.dirname(original_network_path)
        else:
            original_network_file = None

        if os.path.isfile(degenerated_network_path):
            degenerated_network_file = degenerated_network_path
            degenerated_network_path = os.path.dirname(degenerated_network_path)
        else:
            degenerated_network_file = None

        self.original_plot = SnuddaPlotInputLocations(network_path=original_network_path,
                                                      network_file=original_network_file,
                                                      input_file=original_input,
                                                      snudda_data=original_snudda_data)
        self.degenerated_plot = SnuddaPlotInputLocations(network_path=degenerated_network_path,
                                                         network_file=degenerated_network_file,
                                                         input_file=degenerated_input,
                                                         snudda_data=degenerated_snudda_data)

        self.fig_path = os.path.join(degenerated_network_path, "figures")
        if not os.path.exists(self.fig_path):
            os.mkdir(self.fig_path)

    def plot_neuron(self, neuron_id, figure_size=None, show_internal_synapses=True, hide_axis=False):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=figure_size)
        ax = fig.add_subplot(111, projection='3d')
        ax = self.original_plot.plot_neuron_inputs(neuron_id=neuron_id,
                                                   neuron_colour=np.array([0.6, 0.6, 0.6]),
                                                   external_colour=np.array([1, 0.5, 0]),
                                                   internal_colour=np.array([0, 0.5, 1]),
                                                   show_internal_synapses=show_internal_synapses,
                                                   ax=ax,
                                                   size=2,
                                                   save_fig=False,
                                                   show_figure=False,
                                                   hide_axis=hide_axis,
                                                   figure_size=figure_size)

        ax = self.degenerated_plot.plot_neuron_inputs(neuron_id=neuron_id,
                                                      neuron_colour=np.array([0, 0, 0]),
                                                      show_internal_synapses=show_internal_synapses,
                                                      ax=ax,
                                                      size=50,
                                                      save_fig=True,
                                                      hide_axis=hide_axis,
                                                      show_figure=True)


def cli():

    import argparse
    parser = argparse.ArgumentParser("plot_degeneration")
    parser.add_argument("original_network_path", help="Path to original network directory")
    parser.add_argument("degenerated_network_path", help="Path to degenerated network directory")
    parser.add_argument("neuron_id", help="Neuron ID to inspect", type=int)
    args = parser.parse_args()

    pd = PlotDegeneration(original_network_path=args.original_network_path,
                          degenerated_network_path=args.degenerated_network_path)
    pd.plot_neuron(neuron_id=args.neuron_id)


if __name__ == "__main__":
    cli()
