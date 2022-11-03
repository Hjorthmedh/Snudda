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

    def load_data(self):

        self.data_file = h5py.File(self.data_file_name, "r")
        # self.meta_data = np.genfromtxt(self.network_meta_file_name, delimiter=",")
        self.meta_data = pandas.read_csv(self.network_meta_file_name,
                                         names=["neuron_id", "neuron_type", "neuron_name", "x", "y", "z", "morphology"])

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

        row_num = np.array([len(x) for x in self.data_file.values()])
        row_dim = np.array([self.get_dim_from_name(x) for x in self.data_file.keys()])

        # Dimensions might not be in order, reorder them
        sort_idx = np.argsort(row_dim)
        work_load = np.multiply(row_num[sort_idx], row_dim[sort_idx])
        work_load_done = 0
        total_work_load = np.sum(work_load)

        self.open_output()

        for dim, all_rows_on_file in self.data_file.items():

            all_rows = all_rows_on_file[()].copy()

            if self.check_if_exists(dim):
                print(f"{dim} already stored in {self.filtered_file_name}, skipping.")
                continue

            start_time = timeit.default_timer()

            print(f"Processing {dim} with {len(all_rows)}")
            keep_rows = np.zeros((len(all_rows,)), dtype=bool)
            filtered_data = []

            all_rows_flag = np.take(lookup_table, all_rows)
            keep_rows = np.where(np.sum(all_rows_flag, axis=1))[0]

            filtered_data = [all_rows[x] for x in keep_rows]
            print(f"Keeping {len(filtered_data)} out of {len(all_rows)} simplices")
            print(f"Time used for {dim}: {(timeit.default_timer() - start_time):0.1f} seconds")

            self.save_dim_data(dim_name=dim, dim_data=filtered_data)

            dim_value = self.get_dim_from_name(dim)
            work_load_done += work_load[dim_value-1]

            time_used = timeit.default_timer() - start_time
            fraction_done = work_load_done / total_work_load
            total_time_estimated = time_used / fraction_done
            time_left = total_time_estimated - time_used
            sys.stdout.write(f"\rTime used: {time_used:0.1f}s, estimated left {time_left:0.1f}s  "
                             f"({total_time_estimated:0.1f}s)  ")
            sys.stdout.flush()

        self.close()

        print(f"Total time used {(timeit.default_timer() - start_time):0.1f} seconds.")

    def open_output(self):

        if os.path.isfile(self.filtered_file_name):
            self.filtered_data_file = h5py.File(self.filtered_file_name, "a")
        else:
            self.filtered_data_file = h5py.File(self.filtered_file_name, "w")

            meta_info = {"neuron_id": self.meta_data["neuron_id"].to_numpy(),
                         "type": self.meta_data["neuron_type"].to_list(),
                         "name": self.meta_data["neuron_name"].to_list(),
                         "position": self.meta_data[["x", "y", "z"]].to_numpy(),
                         "morphology": self.meta_data["morphology"].to_list(),
                         "num_simplices": np.full((len(self.data_file)+1,), -1, dtype=int)}

            meta = self.filtered_data_file.create_group("meta")
            for key, value in meta_info.items():
                meta.create_dataset(key, data=value)

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
    parser.add_argument("--distance", help="Maximum distance to centre for core neurons (micrometers)", type=float, default=None)

    args = parser.parse_args()

    ff = FilterFlagserData(flagser_data_file_name=args.flagser_file,
                           network_meta_file_name=args.meta_data_file,
                           filtered_file_name=args.output_file)
    ff.load_data()
    ff.find_core_neurons(distance_to_centre=args.distance, population_fraction=args.percentile)


if __name__ == "__main__":
    filter_flagser_cli()