import plotly
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import numpy as np
import os


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
        # pio.renderers.default = "plotly_mimetype"
        pio.renderers.default = "png"

        palette = px.colors.qualitative.Set2

        fig = go.Figure()

        if species is None:
            species = self.sls.list_data_types(neuron_id=neuron_id, exclude=["spikes", "voltage"])

        time = self.sls.get_time()
        all_data = self.sls.get_all_data(neuron_id=neuron_id, exclude=["spikes", "voltage"])

        for i, s in enumerate(species):
            idx = time >= 0.0

            data = all_data[s]

            color = palette[i % len(palette)]

            if compartment_id is None:
                comp_ofs = 0
            else:
                comp_ofs = np.where(data[1][neuron_id][0] == compartment_id)[0]
                if len(comp_ofs) == 0:
                    print(f"Species {s}, not recorded for section_id {compartment_id}")
                    continue
                else:

                    if len(comp_ofs) > 1:
                        print(
                            f"Warning, species {s} recorded at multiple points in section_id {compartment_id}, using first occurence")

                    comp_ofs = comp_ofs[0]

            if species_label is None:
                s_label = s

                if compartment_id != 0:

                    try:
                        section_id = data[1][neuron_id][0][comp_ofs]
                        section_x = data[1][neuron_id][1][comp_ofs]
                        s_label = f"{s} ({section_id}:{section_x})"

                    except Exception as e:
                        import traceback
                        print(traceback.format_exc())
                        print(e)
                        import pdb
                        pdb.set_trace()

            else:
                s_label = species_label[i]

            try:
                # data variable contains 'data', 'sec_id_x', 'syninfo'
                if normalise:
                    fig.add_trace(go.Scatter(x=time[idx],
                                             y=data[0][neuron_id].T[comp_ofs][idx]/np.max(data[0][neuron_id].T[comp_ofs][idx]),
                                             name=s_label, line={"width": 4, "color": color}))
                else:
                    fig.add_trace(go.Scatter(x=time[idx], y=data[0][neuron_id].T[comp_ofs][idx],
                                             name=s_label, line={"width": 4, "color": color}))
            except Exception as e:
                import traceback
                print(traceback.format_exc())
                print(e)
                import pdb
                pdb.set_trace()

        fig.update_layout(xaxis_title="Time (s)", yaxis_title=ylabel, width=1000, height=800,
                          paper_bgcolor="white", plot_bgcolor="white",
                          font={"size": 15},  # General font size for all elements
                          title={"text": title, "font": {"size": 20}, "x": 0.5, "xanchor": "center", "y": 0.9},
                          legend={"font": {"size": 15}},  # Specific font size for legend
                          xaxis={"title": {"font": {"size": 15}}, "tickfont": {"size": 15}},
                          yaxis={"title": {"font": {"size": 15}},
                                 "tickfont": {"size": 15}}, # Y-axis title and tick labels
                          margin=dict(l=100, r=100, t=160, b=100)) # Margin


        if fig_name is not None:
            if fig_path is not None:
                os.makedirs(fig_path, exist_ok=True)
                fig_name = os.path.join(fig_path, fig_name)

            if ".html" in fig_name:
                fig.write_html(fig_name)
            else:
                fig.write_image(fig_name, width=width, height=height)

        fig.show()

        return fig

def plot_cli():
    import argparse

    parser = argparse.ArgumentParser(description="Plot reaction diffusion data from Snudda simulation")
    parser.add_argument("network_path", type=str, help="Path to the network directory")
    parser.add_argument("simulation_file", type=str, help="Path to the simulation output file")
    parser.add_argument("--neuron_id", type=int, required=True, help="Neuron ID to plot")
    parser.add_argument("--species", type=str, nargs="+", default=None, help="Species to plot")
    parser.add_argument("--species_label", type=str, nargs="+", default=None, help="Labels for species")
    parser.add_argument("--ylabel", type=str, default=None, help="Y-axis label")
    parser.add_argument("--compartment_id", type=int, default=0, help="Compartment ID (default: 0, soma)")
    parser.add_argument("--normalise", action="store_true", help="Normalise the concentration data")
    parser.add_argument("--fig_name", type=str, default=None, help="Output figure filename")
    parser.add_argument("--fig_path", type=str, default="figures",
                        help="Output figure directory (default: figures)")
    parser.add_argument("--title", type=str, default=None, help="Plot title")
    parser.add_argument("--width", type=int, default=800, help="Figure width in pixels (default: 800)")
    parser.add_argument("--height", type=int, default=700, help="Figure height in pixels (default: 700)")

    args = parser.parse_args()

    plotter = PlotReactionDiffusion(
        network_path=args.network_path,
        simulation_file=args.simulation_file
    )

    if args.fig_name is None:
        fig_name = f"plot-{args.neuron_id}.html"
    else:
        fig_name = args.fig_name

    plotter.plot(
        neuron_id=args.neuron_id,
        species=args.species,
        species_label=args.species_label,
        ylabel=args.ylabel,
        compartment_id=args.compartment_id,
        normalise=args.normalise,
        fig_name=fig_name,
        fig_path=args.fig_path,
        title=args.title,
        width=args.width,
        height=args.height
    )

if __name__ == "__main__":
    plot_cli()
