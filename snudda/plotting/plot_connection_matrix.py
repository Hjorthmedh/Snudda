import matplotlib.pyplot as plt
import numpy as np
from snudda.utils import SnuddaLoad


class PlotConnectionMatrix:

    def __init__(self, network_path):

        self.network_path = network_path
        self.snudda_load = SnuddaLoad(network_path)

        print("Creating connection matrix...")
        self.con_mat = self.snudda_load.create_connection_matrix(sparse_matrix=False)

    def plot_sorted_connection_matrix(self, fig_path, neuron_order=None):

        all_neuron_id = []

        if neuron_order is None:
            neuron_order = sorted(list(set(self.snudda_load.get_neuron_types())))

        print(f"neuron_order = {neuron_order}")

        label_edges = [0]
        label_centers = []
        labels = []

        for neuron_type in neuron_order:
            print(f"Plotting {neuron_type}")
            neuron_id = self.snudda_load.get_neuron_id_of_type(neuron_type=neuron_type)
            all_neuron_id.append(neuron_id)
            label_centers.append(label_edges[-1] + len(neuron_id)/2)
            label_edges.append(len(neuron_id)+label_edges[-1])
            labels.append(neuron_type)

        all_neuron_id = np.concatenate(all_neuron_id)

        sorted_con_mat = self.con_mat[all_neuron_id, :][:, all_neuron_id]

        cbar_kw = {}
        cbarlabel = ""

        neuron_types = self.snudda_load.get_neuron_types(neuron_id=all_neuron_id)

        fig, ax = plt.subplots()
        # cmap: 'Spectral', 'coolwarm', 'jet
        im = ax.imshow(sorted_con_mat, cmap='Spectral')

        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

        ax.set_xticks(label_centers, labels=labels)
        ax.set_yticks(label_centers, labels=labels)

        #ax.set_xticks(label_edges, minor=True)
        #ax.set_yticks(label_edges, minor=True)
        
        ax.set_xlabel('Postsynaptic neurons')
        ax.set_ylabel('Presynaptic neurons')


        fig.tight_layout()
        print(f"Saving figure to {fig_path}")
        plt.savefig(fig_path)
        plt.show()


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser("Plot connection matrix")
    parser.add_argument("network_path", help="Path to network folder", type=str)
    parser.add_argument("output_file", help="Figure file", type=str)
    parser.add_argument("--order", help="e.g. dSPN,iSPN,FS (no spaces)", type=str)

    args = parser.parse_args()

    if args.order is not None:
        type_order = args.order.split(",")
    else:
        type_order = None

    print(f"type_order: {type_order}")

    pc = PlotConnectionMatrix(network_path=args.network_path)
    pc.plot_sorted_connection_matrix(neuron_order=type_order, fig_path=args.output_file)

