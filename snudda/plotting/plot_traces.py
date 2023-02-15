# python3 plot_traces.py save/simulation/network-output.hdf5 save/network-synapses.hdf5


import sys
import os

import h5py
import numpy as np
from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation
import matplotlib.pyplot as plt

import re
import ntpath
import time


class PlotTraces:

    ############################################################################

    def __init__(self, output_file, network_file=None, input_file=None, experiment_name=None):

        self.output_file = output_file
        self.network_file = network_file
        self.input_file = input_file

        self.time = []
        self.voltage = dict([])

        self.neuron_name_remap = {"FSN": "FS"}

        if experiment_name is not None:
            self.experiment_name = experiment_name
        else:
            self.experiment_name = ""

        if network_file is None and "simulation" in output_file:
            network_path = os.path.dirname(os.path.dirname(output_file))
            network_file = os.path.join(network_path, "network-synapses.hdf5")
            if os.path.exists(network_file):
                self.network_file = network_file

        if network_file is not None:
            network_path = os.path.dirname(os.path.dirname(network_file))
        else:
            network_path = None

        if self.network_file is not None:
            print(f"Loading network info from {self.network_file}")
            self.network_info = SnuddaLoad(self.network_file)
        else:
            self.network_info = None

        if self.input_file is not None:
            print(f"Loading input info from {self.input_file}")
            self.input_info = h5py.File(self.input_file, "r")
        else:
            network_path = os.path.dirname(os.path.dirname(output_file))
            input_file = os.path.join(network_path, "input-spikes.hdf5")
            if os.path.exists(input_file):
                self.input_file = input_file
                print(f"Loading input info from {self.input_file}")
                self.input_info = h5py.File(self.input_file, "r")
            else:
                self.input_info = None

        self.output_load = SnuddaLoadNetworkSimulation(network_simulation_output_file=output_file,
                                                       network_path=network_path)

        self.voltage = self.output_load.get_voltage()
        self.time = self.output_load.get_time()

    ############################################################################

    def neuron_name(self, neuron_type):

        if neuron_type in self.neuron_name_remap:
            return self.neuron_name_remap[neuron_type]
        else:
            return neuron_type

    ############################################################################

    def plot_traces(self, trace_id=None, offset=150e-3, colours=None, skip_time=None, time_range=None,
                    line_width=1, fig_size=None,
                    mark_current=None, mark_current_y=None,
                    title=None, fig_name=None, mark_depolarisation_block=True):

        """
            Plot the traces of neuron trace_id

            Args:
                trace_id (int or list) : ID of trace to show, can be integer or list
                offset (float) : Offset between multiple traces, float or None
                colours : What colour to plot
                skip_time (float) : Skip portion of the start, modifies time shown
                time_range (float, float) : Range to plot
                mark_current (list) : List of tuples of start, end time
                mark_current_y (float) : Y-coordinate of where to mark the current
                title (str) : Plot title
                fig_name (str) : Figure file to save to

        """

        depol_dict = self.output_load.get_depolarisation_dictionary()

        if skip_time is not None:
            print(f"!!! Excluding first {skip_time} s from the plot")

        assert time_range is None or skip_time is None, f"Only specify one of skip_time and time_range"

        if trace_id is None:
            if self.network_info:
                trace_id = [x["neuronID"] for x in self.network_info.data["neurons"]]
            else:
                trace_id = [x for x in self.voltage]
        elif isinstance(trace_id, (int, np.integer)):
            trace_id = [trace_id]

        if colours is None:
            colours = {"dSPN": (77. / 255, 151. / 255, 1.0),
                       "iSPN": (67. / 255, 55. / 255, 181. / 255),
                       "FS": (6. / 255, 31. / 255, 85. / 255),
                       "ChIN": [252. / 255, 102. / 255, 0],
                       "LTS": [150. / 255, 63. / 255, 212. / 255]}

        print(f"Plotting traces: {trace_id}")
        print(f"Plotted {len(trace_id)} traces (total {len(self.voltage)})")

        types_in_plot = set()

        if self.network_info is not None:
            cell_types = [n["type"] for n in self.network_info.data["neurons"]]
            cell_id_check = [n["neuronID"] for n in self.network_info.data["neurons"]]
            try:
                assert (np.array([cell_id_check[x] == x for x in trace_id])).all(), \
                    "Internal error, assume IDs ordered"
            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)
                print("This is strange...")
                import pdb
                pdb.set_trace()

            cols = [colours[c] if c in colours else [0, 0, 0] for c in cell_types]

        fig = plt.figure(figsize=fig_size)

        ofs = 0

        if skip_time is not None:
            time_idx = np.where(self.time >= skip_time)[0]
        else:
            skip_time = 0.0
            time_idx = range(0, len(self.time))

        if time_range is not None:
            time_idx = np.where(np.logical_and(time_range[0] <= self.time, self.time <= time_range[1]))[0]

        plot_count = 0

        for r in trace_id:

            if r not in self.voltage:
                if self.network_info is not None and self.network_info.data["neurons"][r]["virtualNeuron"]:
                    # It was a virtual neuron that lacked voltage, skip warning
                    continue

                print(f"Missing data for trace {r}")
                continue

            plot_count += 1
            types_in_plot.add(self.output_load.network_simulation_file["metaData"]["type"][r].decode())

            if colours is None or self.network_info is None:
                colour = "black"
            else:
                try:
                    colour = cols[r]
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    import pdb
                    pdb.set_trace()

            plt.plot(self.time[time_idx] - skip_time,
                     self.voltage[r][time_idx] + ofs,
                     color=colour, linewidth=line_width)

            if mark_depolarisation_block and r in depol_dict:
                for depol_start_t, depol_end_t in depol_dict[r]:
                    if time_range:
                        if depol_end_t < time_range[0] or time_range[1] < depol_start_t:
                            # Both outside, skip it in plot
                            continue

                        start_in_range = time_range[0] <= depol_start_t and depol_start_t <= time_range[1]
                        end_in_range = time_range[0] <= depol_end_t and depol_end_t <= time_range[1]

                        if start_in_range and not end_in_range:
                            depol_end_t = time_range[1]
                        elif end_in_range and not start_in_range:
                            depol_start_t = time_range[0]

                        if depol_start_t < time_range[0] and time_range[1] < depol_end_t:
                            depol_start_t = time_range[0]
                            depol_end_t = time_range[1]

                    plt.plot([depol_start_t-skip_time, depol_end_t-skip_time],
                             [ofs, ofs], color="red", linewidth=line_width)

            if offset:
                ofs += offset

        if mark_current is not None:
            for t_start, t_end in mark_current:
                plt.plot([t_start - skip_time, t_end - skip_time], [mark_current_y, mark_current_y], 'r-', linewidth=5)

        if plot_count == 0:
            plt.close()
            return

        plt.xlabel('Time')
        plt.ylabel('Voltage')

        if title is None and self.input_info is not None and len(trace_id) == 1:
            n_inputs = 0
            for input_type in self.input_info["input"][str(trace_id[0])]:
                n_inputs += self.input_info["input"][str(trace_id[0])][input_type]["spikes"].shape[0]

            title = f"{self.network_info.data['neurons'][trace_id[0]]['name']} receiving {n_inputs} inputs"

        if title is not None:
            plt.title(title)

        if offset != 0 and offset is not None and len(trace_id) > 1:
            ax = fig.axes[0]
            ax.set_yticklabels([])

        plt.tight_layout()

        fig_path = self.get_figure_path()

        if fig_name is None:
            if len(types_in_plot) > 1:
                fig_name = f"Network-voltage-trace-{self.experiment_name}-{'-'.join(types_in_plot)}.pdf"
            elif len(trace_id) <= 10:
                trace_id_str = '-'.join([str(x) for x in trace_id])
                fig_name = f"Network-voltage-trace-{self.experiment_name}-{types_in_plot.pop()}-{trace_id_str}.pdf"
            else:
                fig_name = f"Network-voltage-trace-{self.experiment_name}-{types_in_plot.pop()}-traces.pdf"

        plt.savefig(os.path.join(fig_path, fig_name), dpi=600)

        print(f"Saving to figure {os.path.join(fig_path, fig_name)}")

        plt.ion()
        plt.show()
        # plt.draw()
        # plt.pause(0.5)  # Show interactive plot (that user can interact with for a short period of time)

        return fig

    def get_figure_path(self):

        if self.network_file:
            fig_path = os.path.join(os.path.dirname(os.path.realpath(self.network_file)), "figures")
        else:
            fig_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(self.output_file))), "figures")

        if not os.path.exists(fig_path):
            os.makedirs(fig_path)

        return fig_path

    ############################################################################

    def plot_trace_neuron_name(self, neuron_name, num_traces=1, skip_time=0, plot_offset=0, fig_name=None,
                               fig_size=None, num_offset=0):
        assert self.network_info is not None, "You need to specify networkInfo file"
        neuron_names = [x["name"] for x in self.network_info.data["neurons"]]
        traceID = [x[0] for x in enumerate(neuron_names) if x[1].lower() == neuron_name.lower()]
        num_traces = min(len(traceID), num_traces)

        if num_traces <= 0:
            print(f"No traces of neuron(s) {neuron_name} to show")
            return

        fig = self.plot_traces(offset=plot_offset, trace_id=traceID[num_offset:num_offset + num_traces],
                               skip_time=skip_time, fig_size=fig_size,
                               title=neuron_names[traceID[0]], fig_name=fig_name)

        time.sleep(1)
        return fig

    def plot_traces_sep(self, trace_id=None, offset=150e-3, colours=None, skip_time=None,
                        title=None, fig_name=None, fig_size=None, folder_name=None):

        if folder_name is None:
            folder_name = ""

        # Plot traces and save as separate images
        if skip_time is not None:
            print(f"!!! Excluding first {skip_time} s from the plot")

        if not trace_id:
            if self.network_info:
                trace_id = [x["neuronID"] for x in self.network_info.data["neurons"]]
            else:
                trace_id = [x for x in self.voltage]
        elif isinstance(trace_id, (int, np.integer)):
            trace_id = [trace_id]

        if colours is None:
            colours = {"dSPN": (77. / 255, 151. / 255, 1.0),
                       "iSPN": (67. / 255, 55. / 255, 181. / 255),
                       "FS": (6. / 255, 31. / 255, 85. / 255),
                       "ChIN": [252. / 255, 102. / 255, 0],
                       "LTS": [150. / 255, 63. / 255, 212. / 255]}

        print(f"Plotting traces: {trace_id}")
        print(f"Plotted {len(trace_id)} traces (total {len(self.voltage)})")

        types_in_plot = set()

        if self.network_info is not None:
            cell_types = [n["type"] for n in self.network_info.data["neurons"]]
            cell_id_check = [n["neuronID"] for n in self.network_info.data["neurons"]]
            try:
                assert (np.array([cell_id_check[x] == x for x in trace_id])).all(), \
                    "Internal error, assume IDs ordered"
            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)
                print("This is strange...")
                import pdb
                pdb.set_trace()

            cols = [colours[c] if c in colours else [0, 0, 0] for c in cell_types]

        if not fig_size:
            fig_size = (10, 5)

        if skip_time is not None:
            time_idx = np.where(self.time >= skip_time)[0]
        else:
            skip_time = 0.0
            time_idx = range(0, len(self.time))

        fig_path = self.get_figure_path()

        if not os.path.exists(os.path.join(fig_path, str(self.experiment_name) + folder_name)):
            os.makedirs(os.path.join(fig_path, str(self.experiment_name) + folder_name))
        plot_count = 0
        for r in trace_id:
            fig = plt.figure(figsize=fig_size)
            if r not in self.voltage:
                print(f"Missing data for trace {r}")
                continue

            plot_count += 1
            types_in_plot.add(self.network_info.data["neurons"][r]["type"])

            plt.plot(self.time[time_idx] - skip_time,1000*self.voltage[r][time_idx], color='black')

            plt.xlabel('Time (s)')
            plt.ylabel('Membrane potential (mV)')

            if title is None and self.input_info is not None and len(trace_id) == 1:
                n_inputs = 0
                for input_type in self.input_info["input"][str(trace_id[0])]:
                    n_inputs += self.input_info["input"][str(trace_id[0])][input_type]["spikes"].shape[0]

                title = f"{self.network_info.data['neurons'][trace_id[0]]['name']} receiving {n_inputs} synaptic inputs"
            title = f"{self.network_info.data['neurons'][trace_id[r]]['name']}"
            plt.title(title)
            plt.tight_layout()
            
            plt.savefig(os.path.join(fig_path, self.experiment_name + folder_name,
                                     f"Network-spikes-{self.experiment_name}-{r}-{title}.png"))
            plt.close(fig)

    ############################################################################

    def plot_trace_neuron_type(self, neuron_type, num_traces=10, offset=0, skip_time=0.0, fig_size=None,
                               mark_depolarisation_block=True):

        assert self.network_info is not None, "You need to specify networkInfo file"

        neuron_types = [x["type"] for x in self.network_info.data["neurons"]]

        # Find numbers of the relevant neurons

        trace_id = [x[0] for x in enumerate(neuron_types) if x[1].lower() == neuron_type.lower()]

        if num_traces is None:
            num_traces = len(trace_id)
        else:
            num_traces = min(len(trace_id), num_traces)

        if num_traces <= 0:
            print(f"No traces of {neuron_type} to show")
            return

        fig = self.plot_traces(offset=offset, trace_id=trace_id[:num_traces], skip_time=skip_time,
                               title=self.neuron_name(neuron_type), fig_size=fig_size,
                               mark_depolarisation_block=mark_depolarisation_block)

        time.sleep(1)
        return fig

    ############################################################################

    def plot_synaptic_currents(self, post_id, fig_size=None):

        """
            Plot synaptic currents impinging on neuron post_id

            Args: post_id (int) : Neuron ID of post synaptic neuron
        """

        data, sec_id_x, syn_info = self.output_load.get_data("synaptic_current", neuron_id=[post_id])
        time = self.output_load.get_time()

        plt.figure(figsize=fig_size)
        line_id_list = []
        for trace, pre_id in zip(data[post_id].T, syn_info[post_id][1]):
            plt.plot(time, trace, label=pre_id)

        plt.legend()
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")
        plt.title(f"Synaptic currents on {post_id} ({self.network_info.data['neurons'][post_id]['name']})")
        plt.ion()
        plt.show()

        fig_path = self.get_figure_path()
        fig_name = os.path.join(fig_path, f"{self.experiment_name}-synaptic-currents-{post_id}.png")
        plt.savefig(fig_name, dpi=300)

    ############################################################################


def snudda_plot_traces_cli():

    import argparse
    parser = argparse.ArgumentParser("Plot traces")
    parser.add_argument("output_file")
    parser.add_argument("--input_file")
    parser.add_argument("--network_file")
    parser.add_argument("--plot_offset", type=float, default=0)
    parser.add_argument("--skip_time", type=float, default=0)
    parser.add_argument("--max_num_traces", type=int, default=None)
    parser.add_argument("--traceID", type=str, default=None, help="Trace ID to plot, separated by comma: e.g. 1,3,14")
    args = parser.parse_args()

    npt = PlotTraces(output_file=args.output_file, network_file=args.network_file, input_file=args.input_file)

    if args.traceID is not None:
        trace_id = [int(x) for x in args.traceID.split(",")]
        npt.plot_traces(offset=args.plot_offset, trace_id=trace_id, skip_time=args.skip_time,
                        title=f"Traces {args.traceID}")
    else:
        for neuron_type in npt.output_load.iter_neuron_type():
            npt.plot_trace_neuron_type(neuron_type=neuron_type, num_traces=args.max_num_traces,
                                       offset=args.plot_offset, skip_time=args.skip_time,
                                       fig_size=(10, 5))


if __name__ == "__main__":
    snudda_plot_traces_cli()