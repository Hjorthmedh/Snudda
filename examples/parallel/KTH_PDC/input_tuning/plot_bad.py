import os
import glob
import json
import numpy as np
from snudda.utils import snudda_parse_path
import matplotlib.pyplot as plt
        

class PlotBad:

    def __init__(self, network_path, figure_path, snudda_data=None):
        self.network_path = network_path
        self.figure_path = figure_path
        self.snudda_data = snudda_data
        self.meta = None
        self.parameters = None

        self.load_config(network_path)
        self.bad_key_list = self.get_keys_in_dir(figure_path)
        self.bad_param_list = [x[1] for x in self.bad_key_list]
        
    def get_keys_in_dir(self, path, file_extension=".png"):

        file_list = glob.glob(os.path.join(path, f"*{file_extension}"))

        key_list = []
        
        for f in file_list:
            f_parts = os.path.basename(f).split("-")
            morph_key = f_parts[0]
            param_key = f_parts[1]

            key_list.append((morph_key, param_key))

        return key_list

    def load_config(self, network_path):

        network_config_path = os.path.join(network_path, "network-config.json")
    
        with open(network_config_path, "r") as f:
            self.network_config = json.load(f)

        if self.snudda_data is None:
            self.snudda_data = self.network_config["snudda_data"]
            
        neuron_paths = []
            
        for region_name, region_data in self.network_config["regions"].items():
            for neuron_name, neuron_data in region_data["neurons"].items():
                for name, path in neuron_data["neuron_path"].items():
                    neuron_paths.append(snudda_parse_path(path, self.snudda_data))
            
        self.load_meta(neuron_paths)
        self.load_parameters(neuron_paths)
        self.extract_parameters()

    def load_meta(self, neuron_paths):

        all_meta = dict()
        read_meta_files = []
        
        for neuron_path in neuron_paths:
            meta_file = os.path.join(neuron_path, "meta.json")

            if meta_file in read_meta_files:
                # Already loaded this meta.json file
                continue
            
            with open(meta_file, "r") as f:
                new_meta = json.load(f)
                read_meta_files.append(meta_file)

                for key in new_meta.keys():
                    if key in all_meta:
                        raise KeyError("Parameter key {key} already exists!")
                
                all_meta |= new_meta

        self.meta = all_meta

    def load_parameters(self, neuron_paths):

        all_params = dict()
        read_param_files = []

        for neuron_path in neuron_paths:
            param_file = os.path.join(neuron_path, "parameters.json")

            if param_file in read_param_files:
                # Already loaded the parameters.json
                continue

            with open(param_file, "r") as f:
                new_params = json.load(f)
                read_param_files.append(param_file)

                for key in new_params.keys():
                    if key in all_params:
                        raise KeyError("Parameter key {key} already exists!")

                all_params |= new_params

        self.parameters = all_params

    def extract_parameters(self):

        self.params = dict()

        for param_key, param_list in self.parameters.items():
            self.params[param_key] = dict()
            
            for param_data in param_list:
                if param_data["type"] == "global":
                    continue
                
                try:
                    param_name = param_data["param_name"]
                    param_loc = param_data["sectionlist"]
                    param_value = param_data["value"]
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()
                    
                self.params[param_key][f"{param_name}_{param_loc}"] = param_value

    def get_param(self, param_name):

        par_val = []
        par_status = []

        for param_key, param_data in self.params.items():
            if param_key in self.bad_param_list:
                status = 0
            else:
                status = 1

            if param_name in param_data:
                par_status.append(status)
                par_val.append(param_data[param_name])
            else:
                print(f"{param_name = } missing for {param_key = }")
                
        return par_val, par_status

    def get_all_param_names(self):
        param_names = []
        for par_key, par_data in self.params.items():
            for par_name in par_data.keys():                
                param_names.append(par_name)

        return np.unique(param_names)
    
    def plot_params(self, param_name):

        param_val, param_status = self.get_param(param_name)
        color = [(0,0,0) if x == 1 else (1,0,0) for x in param_status]
        plt.figure()
        plt.scatter(param_status, param_val, c=color)
        plt.title(param_name)
        plt.ion()
        plt.show()
        fig_path = os.path.join(self.network_path, "figures", f"{param_name}-bad-summary.png")
        plt.savefig(fig_path)

    def summary_plot(self, normalise=True):
        param_names = self.get_all_param_names()

        plot_x = []
        plot_y = []

        plot_ctr = 0

        label_x = []
        label_str = []
        spacing = 6
        
        for param_name in param_names:

            param_val, param_status = self.get_param(param_name)

            if normalise:
                min_val = np.min(param_val)
                max_val = np.max(param_val)

                param_val = [(x-min_val)/(max_val-min_val) for x in param_val]
            
            if len(np.unique(param_val)) == 1:
                # All values same, skip.
                continue

            plot_x += [x+plot_ctr*spacing for x in param_status]
            plot_y += param_val

            label_x.append(plot_ctr*spacing)
            label_str.append(param_name)
            
            plot_ctr += 1

            
        color = [(0,0,0) if x % 2 == 1 else (1,0,0) for x in plot_x]
        sizes = [5 if x % 2 == 1 else 20 for x in plot_x]
        
        plt.figure(figsize=(15,8))
        plt.scatter(plot_x, plot_y, c=color, s=sizes)
        plt.xticks(label_x, label_str, rotation=90)
        plt.subplots_adjust(bottom=0.4)
        title_str = self.network_path
        if normalise:
            title_str += " (normalised values plotted)"
        plt.title(title_str)
        plt.ion()
        plt.show()

        fig_path = os.path.join(self.network_path, "figures", f"bad-summary.png")        
        plt.savefig(fig_path)
            
            
    def plot(self):
        param_names = self.get_all_param_names()

        for pn in param_names:
            self.plot_params(pn)

        
def cli():

    import argparse
    parser = argparse.ArgumentParser(description="Investigate parameters for plots")
    parser.add_argument("network_path", help="Path to network folder")
    parser.add_argument("figure_path", help="Path to figure folder")
    parser.add_argument("--snudda_data", type=str, default=None)
    args = parser.parse_args()
    
    pb = PlotBad(network_path=args.network_path, figure_path=args.figure_path,
                 snudda_data=args.snudda_data)
    # pb.plot()
    pb.summary_plot()


if __name__ == "__main__":
    cli()
    input("Press a key to exit")
