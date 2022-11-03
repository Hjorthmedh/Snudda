import os
import sys
import numpy as np
import h5py
import pandas
import timeit


class FilterFlagserData:

    def __init__(self, flagser_data_file_name, network_meta_file_name, filtered_file_name=None):

        self.data_file_name = flagser_data_file_name
        self.network_meta_file_name = network_meta_file_name

        self.data_file = None
        self.meta_data = None
        self.filtered_data_file = None

        if filtered_file_name:
            self.filtered_file_name = filtered_file_name
        else:
            self.filtered_file_name = f"{flagser_data_file_name}-filtered.hdf5"

    def close(self):
        if self.data_file:
            self.data_file.close()
            self.data_file = None

        if self.filtered_data_file:
            self.filtered_data_file.close()
            self.filtered_data_file = None

    def load_data(self):

        self.data_file = open(self.data_file_name, "r")
        # self.meta_data = np.genfromtxt(self.network_meta_file_name, delimiter=",")
        self.meta_data = pandas.read_csv(self.network_meta_file_name,
                                         names=["neuron_id", "neuron_type", "neuron_name", "x", "y", "z", "morphology", "population_id"])

    def find_core_neurons(self, population_fraction=None, distance_to_centre=None):

        neuron_id = self.meta_data["neuron_id"].to_numpy()
        positions = self.meta_data[["x", "y", "z"]].to_numpy()

        centre_point = np.mean(positions, axis=0)
        dist_to_centre = np.linalg.norm(positions - centre_point, axis=1)

        if distance_to_centre:
            print(f"Core filtering, maximal distance to centre {distance_to_centre} meters")
            distance_threshold = distance_to_centre
        else:
            assert population_fraction is not None, f"Either distance_to_centre or population_fraction must be set"

            print(f"Filtering simplices with at least one member in core "
                  f"(core fraction {population_fraction} of all neurons)")
            distance_threshold = np.percentile(dist_to_centre, population_fraction)

        centre_idx = np.where(dist_to_centre <= distance_threshold)[0]

        lookup_table = np.zeros((len(neuron_id),), dtype=bool)
        lookup_table[centre_idx] = True

        self.filtered_data_file = open(self.filtered_file_name, "w")

        start_time = timeit.default_timer()

        keep_count = 0
        total_count = 0

        self.filtered_data_file.write("# Centre IDs:" + ", ".join([str(x) for x in neuron_id[centre_idx]]) + "\n")

        for row in self.data_file:

            neuron_id = np.array([int(x) for x in row.strip().split(" ")])

            keep_flag = np.sum(np.take(lookup_table, neuron_id)) > 0
            if keep_flag:
                self.filtered_data_file.write(",".join([f"{x}" for x in neuron_id]) + "\n")
                keep_count += 1

            total_count += 1

        self.close()

        print(f"Total time used {(timeit.default_timer() - start_time):0.1f} seconds.")
        print(f"Kept {keep_count} out of {total_count} rows")

    @staticmethod
    def get_dim_from_name(dim_name):
        return int(dim_name[6:])

    def check_if_exists(self, dim_name):
        return dim_name in self.filtered_data_file and \
               self.filtered_data_file["meta/num_simplices"][self.get_dim_from_name(dim_name)] >= 0

    def save_dim_data(self, dim_name, dim_data):
        print(f"Writing {dim_name} ({len(dim_data)} rows)")
        self.filtered_data_file.create_dataset(dim_name, data=dim_data)
        self.filtered_data_file["meta/num_simplices"][self.get_dim_from_name(dim_name)] = len(dim_data)


def filter_flagser_cli():

    import argparse
    parser = argparse.ArgumentParser("Filter flagser")
    parser.add_argument("flagser_file", help="Flagser file")
    parser.add_argument("meta_data_file", help="Meta data file")
    parser.add_argument("output_file", help="Output file name")
    parser.add_argument("--percentile", help="Percentile of population kept (0-100)", type=int, default=None)
    parser.add_argument("--distance", help="Maximum distance to centre for core neurons (micrometers)",
                        type=float, default=None)

    args = parser.parse_args()

    ff = FilterFlagserData(flagser_data_file_name=args.flagser_file,
                           network_meta_file_name=args.meta_data_file,
                           filtered_file_name=args.output_file)
    ff.load_data()
    ff.find_core_neurons(distance_to_centre=args.distance, population_fraction=args.percentile)


if __name__ == "__main__":
    filter_flagser_cli()