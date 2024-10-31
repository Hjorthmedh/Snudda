import plotly.graph_objects as go
from plotly.subplots import make_subplots

import numpy as np
from snudda.utils import SnuddaLoadSimulation
from snudda.utils.load import SnuddaLoad

# To get clean plots:
#
# import plotly.io as pio
# pio.templates.default = "simple_white"
#


class PlotTrace:

    def __init__(self, snudda_load_simulation=None, network_simulation_file=None):

        self.snudda_load_simulation = None
        self.colour_map = dict()
        self.traces = dict()
        self.default_colour = None
        self.default_width = 1

        if snudda_load_simulation is not None:
            self.snudda_load_simulation = snudda_load_simulation

        if network_simulation_file is not None:
            self.load_data(network_simulation_file=network_simulation_file)

    def load_data(self, network_simulation_file, do_depol_block_test=False):
        self.snudda_load_simulation = SnuddaLoadSimulation(network_simulation_output_file=network_simulation_file,
                                                           quiet_load=True, do_test=do_depol_block_test)

    def define_colour_by_neuron_id(self, colour_by_neuron_id):
        self.colour_map |= colour_by_neuron_id

    def define_colour_by_neuron_type(self, colour_by_neuron_type):
        neuron_names = [SnuddaLoad.to_str(x) for x
                        in self.snudda_load_simulation.network_simulation_file["meta_data"]["name"]]
        neuron_types = [x.split("_")[0] for x in neuron_names]

        for neuron_id, nt in enumerate(neuron_types):
            if nt in colour_by_neuron_type:
                self.colour_map[neuron_id] = colour_by_neuron_type[nt]

    def define_colour_by_neuron_name(self, colour_by_name):
        neuron_names = self.snudda_load_simulation.get_neuron_name()
        for neuron_id, name in enumerate(neuron_names):
            if name in colour_by_name:
                self.colour_map[neuron_id] = colour_by_name[name]

    def clear_colour(self):
        self.colour_map = dict()

    def get_instant_frequency(self, neuron_id):
        # neuron_id is assumed to be scalar int here
        spike_times = self.snudda_load_simulation.get_spikes(neuron_id=neuron_id).flatten()

        isi = 1.0 / np.diff(spike_times)
        isi_time = spike_times[1:]

        return isi_time, isi

    def get_voltage_trace(self, neuron_id):

        volt = self.snudda_load_simulation.get_voltage(neuron_id=neuron_id)
        time = self.snudda_load_simulation.get_time()

        return time, volt[:, 0]

    def create_trace(self, neuron_id):

        time, volt = self.get_voltage_trace(neuron_id=neuron_id)
        isi_time, isi = self.get_instant_frequency(neuron_id=neuron_id)

        colour = self.colour_map.get(neuron_id, self.default_colour)
        neuron_name = self.snudda_load_simulation.get_neuron_name(neuron_id=[neuron_id])[0]
        label = f"{neuron_name} ({neuron_id})"

        isi_sct = go.Scatter(x=isi_time.T, y=isi.T, mode="lines",
                             line={"color": colour, "width": self.default_width},
                             name=label, legendgroup=label, showlegend=False)

        volt_sct = go.Scatter(x=time.T, y=volt.T, mode="lines",
                              line={"color": colour, "width": self.default_width},
                              name=label, legendgroup=label)

        if neuron_id not in self.traces:
            self.traces[neuron_id] = dict()

        self.traces[neuron_id]["isi"] = isi_sct
        self.traces[neuron_id]["volt"] = volt_sct

    def create_traces(self, neuron_id_list=None):
        if neuron_id_list is None:
            neuron_id_list = [int(x) for x
                              in self.snudda_load_simulation.network_simulation_file["neurons"].keys()]

        for neuron_id in neuron_id_list:
            self.create_trace(neuron_id=neuron_id)

    def make_figure(self, show=False):

        fig = make_subplots(cols=1, rows=2, shared_xaxes=True)

        for neuron_id, plots in self.traces.items():
            if neuron_id not in self.colour_map:
                raise KeyError(f"No colour defined for {neuron_id = } or its neuron type")
            fig.add_trace(plots["volt"], row=1, col=1)
            fig.add_trace(plots["isi"], row=2, col=1)

        fig.update_layout(xaxis={"title": "Time (s)"},
                          yaxis1={"title": "Volt (V)"},
                          yaxis2={"title": "Frequency (Hz)"})

        if show:
            fig.show()

        return fig

    def plot_traces(self, neuron_id=None, show=False):
        if isinstance(neuron_id, (int, np.integer)):
            neuron_id = [neuron_id]

        self.create_traces(neuron_id_list=neuron_id)
        fig = self.make_figure(show=show)
        return fig

    # To save: fig.write_image("yourfigname.png")