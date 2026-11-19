# This code loads all the setup from a JSON file, and the simulates it
# import os
import json
import shutil
import subprocess
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from snudda.simulate.network_pair_pulse_simulation import SnuddaNetworkPairPulseSimulation


class RunMultiplePairRecordingExperiments:

    def __init__(self, experiment_config):

        """ Running pair recording experiments defined in experiment_config.


        """

        if isinstance(experiment_config, dict):
            self.experiment_config = experiment_config
        else:
            print(f"Loading {experiment_config}")
            with open(experiment_config, "rt") as f:
                self.experiment_config = json.load(f)

        self.pps = None

    def process(self, exp_key=None):

        if exp_key is None:
            for exp_key in self.experiment_config.keys():
                self.process(exp_key)
            return

        self.setup_slice(exp_key)
        self.run_slice(exp_key)
        self.analyse(exp_key)

    def setup_slice(self, exp_key):

        exp_info = self.experiment_config[exp_key]
        exp_setup = exp_info["setup"]
        if isinstance(exp_setup["post_type"], list):
            post_type = "ALL"
        else:
            post_type = exp_setup["post_type"]

        self.pps = SnuddaNetworkPairPulseSimulation(network_path=exp_info["network_path"],
                                               snudda_data=exp_info["snudda_data"],
                                               exp_type=exp_info["exp_name"],
                                               pre_type=exp_setup["pre_type"],
                                               post_type=post_type,
                                               n_stimulated_neurons=exp_setup["n_stimulated_neurons"],
                                               max_dist=exp_setup["max_dist"],
                                               hold_voltage=exp_setup["hold_v"],
                                               current_injection=exp_setup["cur_inj"],
                                               random_seed=exp_setup["random_seed"])

        neuron_names, neuron_paths, number_of_neurons = zip(*[(k,v["neuron_path"], v["num"]) for k, v in exp_info["neurons"].items()])

        self.pps.setup(neuron_paths=neuron_paths, neuron_names=neuron_names,
                       volume_type="slice",
                       slice_side_len=exp_info["slice"]["length"],
                       slice_depth=exp_info["slice"]["depth"],
                       number_of_neurons=number_of_neurons,
                       connection_config=exp_info["connection_config"], parallel=True)

        self.pps.ablate_post(remove_pre_without_targets=True)


    def run_slice(self, exp_key):

        exp_info = self.experiment_config[exp_key]
        exp_setup = exp_info["setup"]


        # import subprocess

        run_str = self.pps.write_simulation_config(gaba_rev=exp_setup["GABA_rev"],
                                                   clamp_mode=exp_setup["clamp_mode"],
                                                   return_run_str=True)

        # On Cray/Dardel mpirun is not available; use srun instead.
        if shutil.which("mpirun") is None and shutil.which("srun") is not None:
            run_str = run_str.replace("mpirun", "srun", 1)

        print(f"Running:\n   {run_str}")
        result = subprocess.run(run_str, shell=True)
        if result.returncode != 0:
            # srun on Cray/Dardel returns non-zero at MPI finalize even on success;
            # proceed to analyse() so data is not silently dropped.
            print(f"Warning: srun returned exit code {result.returncode}, attempting to analyse anyway.")

    def analyse(self, exp_key, post_type=None):
        exp_info = self.experiment_config[exp_key]
        post_type = post_type or exp_info["setup"]["post_type"]

        if isinstance(post_type, list):
            return all([self.analyse(exp_key=exp_key, post_type=pt) for pt in post_type])

        exp_anal = exp_info.get("analyse", {})

        filter_small_amps = exp_anal.get("filter_small_amps", False)
        filter_spiking = exp_anal.get("filter_spiking", False)
        spike_threshold = exp_anal.get("spike_threshold", 0)
        time_window = exp_anal.get("time_window", 0.050)

        model_mean, model_sd, trace_fig, hist_fig = \
            self.pps.analyse(post_type=post_type,
                             exp_data_file=exp_info["experiment"],
                             clamp_mode=exp_info["setup"]["clamp_mode"],
                             time_window = time_window,
                             spike_threshold = spike_threshold,
                             filter_spiking = filter_spiking,
                             filter_small_amps=filter_small_amps)


        pre_type = exp_info["setup"]["pre_type"]
        exp_name = exp_info["exp_name"]
        exp_text = ""
        is_ok = False

        if exp_name in exp_info["experiment"]:
            exp_text += f"{pre_type} -> {post_type}: "
            exp_mean, exp_std = exp_info["experiment"][exp_name][pre_type][post_type]
            if model_mean is None or exp_mean is None:
                exp_text += "NO DATA.\n"
                is_ok = False
            elif abs(model_mean - exp_mean) < exp_std:
                exp_text += "Success.\n"
                is_ok = True
            else:
                exp_text += "FAILED.\n"
                is_ok = False

        with open("pair_summary.txt", "a") as f:
            f.write(exp_text)

        plt.close("all")

        return is_ok

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("Slice pair pulse experiment")
    parser.add_argument("exp_file", type=str)
    parser.add_argument("--exp_key", type=str, default=None)
    args = parser.parse_args()

    mpr = RunMultiplePairRecordingExperiments(experiment_config=args.exp_file)

    mpr.process(exp_key=args.exp_key)