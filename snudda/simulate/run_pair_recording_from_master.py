"""Run pair/population recording experiments using a pre-built master network.

Each experiment sub-slices the master (via SnuddaCut), ablates irrelevant neurons,
optionally marks pre-neurons as virtual (spike injection instead of morphology sim),
then runs the simulation.  This skips the expensive place/detect/prune phase.

The JSON config format is the same as pair_recording_experiment.json and
pop_recording_experiment.json, with these extra top-level fields per entry:
  - "slice": {"length": <m>, "depth": <m>}   — required for master workflow
  - "setup": {
      ...existing fields...,
      "make_pre_virtual": false,   # optional (default false)
      "stim_all_at_once": false,   # optional (default false)
    }
  - Entries that have "volume": {"type": "default"} (no slice) fall back to
    the original setup() flow automatically.

Usage:
    python run_pair_recording_from_master.py <config.json> <master_network_path>
         [--exp_key KEY] [--make_pre_virtual]
"""

import argparse
import json
import shutil
import subprocess
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from snudda.simulate.network_pair_pulse_simulation import SnuddaNetworkPairPulseSimulation


class RunPairRecordingFromMaster:

    def __init__(self, experiment_config, master_network_path):
        if isinstance(experiment_config, dict):
            self.experiment_config = experiment_config
        else:
            print(f"Loading {experiment_config}")
            with open(experiment_config, "rt") as f:
                self.experiment_config = json.load(f)

        self.master_network_path = master_network_path
        self.pps = None

    def process(self, exp_key=None, make_pre_virtual=False):
        if exp_key is None:
            for key in self.experiment_config.keys():
                if key.startswith("_"):
                    continue
                self.process(exp_key=key, make_pre_virtual=make_pre_virtual)
            return

        self.setup_slice(exp_key, make_pre_virtual=make_pre_virtual)
        self.run_slice(exp_key)
        self.analyse(exp_key)

    def setup_slice(self, exp_key, make_pre_virtual=False):
        exp_info = self.experiment_config[exp_key]
        exp_setup = exp_info["setup"]

        post_type = exp_setup["post_type"]
        if isinstance(post_type, list):
            post_type_pps = "ALL"
        else:
            post_type_pps = post_type

        self.pps = SnuddaNetworkPairPulseSimulation(
            network_path=exp_info["network_path"],
            snudda_data=exp_info.get("snudda_data"),
            exp_type=exp_info["exp_name"],
            pre_type=exp_setup["pre_type"],
            post_type=post_type_pps,
            n_stimulated_neurons=exp_setup.get("n_stimulated_neurons"),
            max_dist=exp_setup.get("max_dist"),
            hold_voltage=exp_setup["hold_v"],
            current_injection=exp_setup["cur_inj"],
            random_seed=exp_setup.get("random_seed"))

        use_master = "slice" in exp_info and "volume" not in exp_info

        if use_master:
            self.pps.setup_from_master(
                master_network_path=self.master_network_path,
                slice_side_len=exp_info["slice"]["length"],
                slice_depth=exp_info["slice"]["depth"])
        else:
            # Fallback: original flow for default-volume experiments
            neuron_names, neuron_paths, number_of_neurons = zip(
                *[(k, v["neuron_path"], v["num"])
                  for k, v in exp_info["neurons"].items()])

            volume_type = exp_info.get("volume", {}).get("type", "slice")
            if volume_type == "slice":
                self.pps.setup(neuron_paths=neuron_paths,
                               neuron_names=neuron_names,
                               volume_type="slice",
                               slice_side_len=exp_info["slice"]["length"],
                               slice_depth=exp_info["slice"]["depth"],
                               number_of_neurons=number_of_neurons,
                               connection_config=exp_info["connection_config"],
                               parallel=True)
            else:
                self.pps.setup(neuron_paths=neuron_paths,
                               neuron_names=neuron_names,
                               number_of_neurons=number_of_neurons,
                               connection_config=exp_info["connection_config"],
                               parallel=True)

        # make_pre_virtual: allow JSON override or CLI flag
        make_virt = make_pre_virtual or exp_setup.get("make_pre_virtual", False)
        self.pps.ablate_post(remove_pre_without_targets=True,
                             make_pre_virtual=make_virt)

    def run_slice(self, exp_key):
        exp_info = self.experiment_config[exp_key]
        exp_setup = exp_info["setup"]

        stim_all_at_once = exp_setup.get("stim_all_at_once", False)
        reversal_potential = exp_setup.get("reversal_potential")
        gaba_rev = exp_setup.get("GABA_rev")

        run_str = self.pps.write_simulation_config(
            gaba_rev=gaba_rev,
            reversal_potential=reversal_potential,
            clamp_mode=exp_setup["clamp_mode"],
            stim_all_at_once=stim_all_at_once,
            return_run_str=True)

        if shutil.which("mpirun") is None and shutil.which("srun") is not None:
            run_str = run_str.replace("mpirun", "srun", 1)

        print(f"Running:\n   {run_str}")
        subprocess.run(run_str, shell=True, check=True)

    def analyse(self, exp_key, post_type=None):
        exp_info = self.experiment_config[exp_key]
        post_type = post_type or exp_info["setup"]["post_type"]

        if isinstance(post_type, list):
            return all([self.analyse(exp_key=exp_key, post_type=pt)
                        for pt in post_type])

        exp_anal = exp_info.get("analyse", {})
        time_window = exp_anal.get("time_window", 0.050)
        filter_small_amps = exp_anal.get("filter_small_amps", False)
        filter_spiking = exp_anal.get("filter_spiking", False)
        spike_threshold = exp_anal.get("spike_threshold", 0)

        model_mean, model_sd, trace_fig, hist_fig = self.pps.analyse(
            post_type=post_type,
            exp_data_file=exp_info["experiment"],
            clamp_mode=exp_info["setup"]["clamp_mode"],
            time_window=time_window,
            spike_threshold=spike_threshold,
            filter_spiking=filter_spiking,
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
    parser = argparse.ArgumentParser(
        description="Run pair/population recording experiments from a master network")
    parser.add_argument("exp_file", type=str,
                        help="JSON experiment config (pair_recording_experiment.json or pop_recording_experiment.json)")
    parser.add_argument("master_network_path", type=str,
                        help="Path to pre-built master network (contains network-synapses.hdf5)")
    parser.add_argument("--exp_key", type=str, default=None,
                        help="Run only this experiment key (default: all)")
    parser.add_argument("--make_pre_virtual", action="store_true",
                        help="Mark pre-synaptic neurons as virtual (spike injection, no morphology sim)")
    args = parser.parse_args()

    runner = RunPairRecordingFromMaster(
        experiment_config=args.exp_file,
        master_network_path=args.master_network_path)

    runner.process(exp_key=args.exp_key, make_pre_virtual=args.make_pre_virtual)
