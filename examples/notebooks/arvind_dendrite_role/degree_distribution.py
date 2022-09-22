import os
import numpy as np
from snudda.utils import SnuddaLoad


class DegreeDistribution:

    def __init__(self, network_file):

        self.network_file = network_file
        self.network_loader = SnuddaLoad(network_file=network_file, load_synapses=True)
        self.connection_matrix = self.network_loader.create_connection_matrix(sparse_matrix=False)
        self.degree_distribution = dict()

    def get_degree_distribution(self, pre_neuron, post_neuron):

        pre_neuron_id = self.network_loader.get_neuron_id_of_type(neuron_type=pre_neuron)
        post_neuron_id = self.network_loader.get_neuron_id_of_type(neuron_type=post_neuron)

        degree_distribution = np.zeros((len(post_neuron_id) + 1), dtype=np.uint)

        sub_matrix = self.connection_matrix[pre_neuron_id, :][:, post_neuron_id]
        for pre_id_A in pre_neuron_id:
            for pre_id_B in pre_neuron_id:
                if pre_id_A == pre_id_B:
                    continue

                # If connection matrix was sparse
                # con_vect_A = (self.connection_matrix[pre_id_A, :][:, post_neuron_id] > 0).toarray()
                # con_vect_B = (self.connection_matrix[pre_id_B, :][:, post_neuron_id] > 0).toarray()

                con_vect_A = (self.connection_matrix[pre_id_A, :][post_neuron_id] > 0)
                con_vect_B = (self.connection_matrix[pre_id_B, :][post_neuron_id] > 0)

                degree = np.sum(np.logical_and(con_vect_A, con_vect_B))
                degree_distribution[degree] += 1

        degree_distribution = np.trim_zeros(degree_distribution, trim="b")

        return degree_distribution

    def get_all_degree_distributions(self):

        neuron_types = self.network_loader.get_neuron_types(return_set=True)
        degree_distributions = dict()

        for pre_neuron in neuron_types:
            for post_neuron in neuron_types:
                print(f"Processing {pre_neuron} -> {post_neuron}")
                deg_dist = self.get_degree_distribution(pre_neuron=pre_neuron, post_neuron=post_neuron)

                if deg_dist.size > 0:
                    degree_distributions[pre_neuron, post_neuron] = deg_dist

        return degree_distributions

    def write_all_degree_distributions(self, plot=True):

        file_name = os.path.join(os.path.dirname(self.network_file),
                                 f"degree-distribution.csv")

        neuron_types = self.network_loader.get_neuron_types(return_set=True)

        print(f"Writing degree distribution to {file_name}")

        degree_distributions = self.get_all_degree_distributions()

        with open(file_name, "w") as f:
            for (pre_neuron, post_neuron), distribution in degree_distributions.items():
                f.write(f"{pre_neuron}, {post_neuron}, ")
                f.write(", ".join([f'{x}' for x in distribution]))
                f.write("\n")

        self.plot_degree_distributions(degree_distributions=degree_distributions)

    def plot_degree_distributions(self, degree_distributions=None):

        if degree_distributions is None:
            degree_distributions = self.get_all_degree_distributions()

        import matplotlib.pyplot as plt
        fig = plt.figure()

        plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(degree_distributions)))))

        for (pre_neuron, post_neuron), distribution in degree_distributions.items():
            degree = np.arange(0, len(distribution))
            plt.plot(degree, distribution, label=f"{pre_neuron}->{post_neuron}")

        plt.legend(prop={'size': 6})
        plt.xlabel("Degree")
        plt.ylabel("Count")

        fig_name = os.path.join(os.path.dirname(self.network_file), "degree_distribution.png")
        plt.savefig(fig_name, dpi=300)
        print(f"Writing figure to {fig_name}")


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser("Extract degree distributions")
    parser.add_argument("network_file")

    args = parser.parse_args()

    dd = DegreeDistribution(args.network_file)
    dd.write_all_degree_distributions()

