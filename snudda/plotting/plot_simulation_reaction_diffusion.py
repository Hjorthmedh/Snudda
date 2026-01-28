import os
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import numpy as np

from snudda.utils import SnuddaLoadSimulation

class PlotReactionDiffusion:

    def __init__(self, network_path, simulation_file, snudda_simulation=None):

        if snudda_simulation is None:
            self.sls = SnuddaLoadSimulation(network_path=network_path,
                                            network_simulation_output_file=simulation_file)
        else:
            self.sls = snudda_simulation

    def list_neuron_info(self, neuron_id):
        return self.sls.list_data_types(neuron_id=neuron_id)

    def plot(self, neuron_id, species=None, species_label=None,
             ylabel=None,
             compartment_id = 0,
             normalise=False,
             fig_name=None, fig_path="figures", title=None,
             width=800, height=700):

        """ compartment_id is based on the order the compartments are added, ie. 0 is first one added, usually soma"""

        if ylabel is None:
            if normalise:
                ylabel = "Concentration (normalised)"
            else:
                ylabel = "Concentration (mM)"

        # pio.renderers.default = "iframe"  # Do not save plots in the notebook, they can get BIG
        pio.renderers.default = "plotly_mimetype"

        palette = px.colors.qualitative.Set2

        fig = go.Figure()

        if species is None:
            species = self.sls.list_data_types(neuron_id=neuron_id)

        time = self.sls.get_time()
        all_data = self.sls.get_all_data(neuron_id=neuron_id, exclude=["spikes", "voltage"])

        for i, s in enumerate(species):
            idx = time >= 0.0

            data = all_data[s]

            color = palette[i % len(palette)]

            if species_label is None:
                s_label = s
            else:
                s_label = species_label[i]

            try:
                # data variable contains 'data', 'sec_id_x', 'syninfo'
                if normalise:
                    fig.add_trace(go.Scatter(x=time[idx],
                                             y=data[0][neuron_id].T[compartment_id][idx]/np.max(data[0][neuron_id].T[compartment_id][idx]),
                                             name=s_label, line={"width": 4, "color": color}))
                else:
                    fig.add_trace(go.Scatter(x=time[idx], y=data[0][neuron_id].T[compartment_id][idx],
                                             name=s_label, line={"width": 4, "color": color}))
            except Exception as e:
                import traceback
                print(traceback.format_exc())
                print(e)
                import pdb
                pdb.set_trace()

        fig.update_layout(xaxis_title="Time (s)", yaxis_title=ylabel, width=1000, height=800,
                          paper_bgcolor="white", plot_bgcolor="white",
                          font={"size": 18},  # General font size for all elements
                          title={"text": title, "font": {"size": 60}, "x": 0.5, "xanchor": "center", "y": 0.9},
                          legend={"font": {"size": 40}},  # Specific font size for legend
                          xaxis={"title": {"font": {"size": 40}}, "tickfont": {"size": 40}},
                          yaxis={"title": {"font": {"size": 40}},
                                 "tickfont": {"size": 40}}, # Y-axis title and tick labels
                          margin=dict(l=100, r=100, t=160, b=100)) # Margin


        if fig_name is not None:
            if fig_path is not None:
                os.makedirs(fig_path, exist_ok=True)
                fig_name = os.path.join(fig_path, fig_name)

            fig.write_image(fig_name, width=width, height=height)

        fig.show()

        return fig