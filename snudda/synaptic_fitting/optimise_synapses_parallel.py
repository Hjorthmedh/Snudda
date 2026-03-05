# TODO: 2026-02-26 -- Dear future humans.
#      -- We have copied optimise_synapse_full.py (but not all of it)
#      -- we need to write code to run on self.pc.id() == 0 that sets up the optimisation
#          and performs all the goodie stuff
#
#      We got to 'get_refined_parameters' function, probably not needed?
#      We are going to use mpirun and neuron parallel, to run multiple instances at the same time.

# We want to run synapse optimisation in parallel
# 1. Start using mpirun (to get multiple instances)
# 2. Worker id 0 (master) sets up optimisation, bounds etc.
# 3. Setup models (with default parameters, will be overwritten each iteration)
# 4. Start optimisation loop on 0.
# 5. Get the paramater sets to investigate
# 6. Distribute parameters to all the workers, set them to the synapses
# 7. Simulate the simulation
# 8. Gather the results master node
# 9. Update the opt object with results
# 10. Repeat 5-10

import os
import sys
import shutil
import timeit

import numpy as np
import json
import copy
import time

import scipy

from mpi4py import MPI  # This must be imported before neuron, to run parallel
import neuron
from neuron import h  # , gui
from snudda.utils.snudda_path import snudda_parse_path, get_snudda_data
from snudda.synaptic_fitting.parameter_bookkeeper import ParameterBookkeeper

from run_synapse_run import RunSynapseRun


class SynapseOptimiser:

    def __init__(self, data_file,
                 entropy=1023456734529028340264793840,
                 synapse_type="glut",
                 load_parameters=True,
                 snudda_data=None,
                 neuron_set_file="neuron_set.json",
                 synapse_parameter_file=None,
                 verbose=True):

        self.log_file = None
        self.verbose = verbose
        self.rng = None


        self.data_file = data_file
        self.neuron_set_file = neuron_set_file
        self.seed = None
        self.entropy=entropy

        self.snudda_data = get_snudda_data(snudda_data=snudda_data)

        self.rsr_synapse_model = None
        self.synapse_type = synapse_type
        self.data = None
        self.volt = None
        self.time = None
        self.sample_freq = None
        self.trace_holding_voltage = None
        self.stim_time = None
        self.cell_type = None
        self.neuron_set_file = None

        self.synapse_parameter_data = None
        self.synapse_section_id = None
        self.synapse_section_x = None

        self.cell_properties = None

        self.load_trace_data()

        # if load_parameters:
        #     self.load_parameter_data()

        self.synapse_parameter_file = synapse_parameter_file

        if synapse_parameter_file:
            with open(synapse_parameter_file, 'r') as f:
                self.write_log(f"Reading synapse parameters from {synapse_parameter_file}")
                self.synapse_parameters = json.load(f)["data"]
        else:
            self.synapse_parameters = {}

        self.parameter_data_file_name = f"{self.data_file}-parameters-full.json"



        self.pc = h.ParallelContext()
        self.n_workers = self.pc.nhost()

        print(f"I am {self.pc.id()}/{self.pc.nhost()}")
        if self.pc.id() == 0:
            rng = np.random.default_rng(seed=0)
            param_list = [rng.random(size=(3, 1)) for x in range(self.pc.nhost())]
            print(f"Before sending: {param_list = }")
        else:
            param_list = []

        self.param = self.pc.py_scatter(param_list)

        self.pc.barrier()
        print(f"Done {self.pc.id()} {self.param.T = }")

    def setup_rng(self):

        if self.rng is not None:
            print(f"setup_rng: rng already setup, skipping.")
            return

        seeds = []

        if self.pc.id() == 0:
            # Setup and distribute random seeds to all workers
            seed_sequence = np.random.SeedSequence(entropy=self.entropy)
            seeds=seed_sequence.generate_state(self.n_workers)

        self.seed = self.pc.py_scatter(seeds)
        self.rng = np.random.default_rng(seed=self.seed)

        self.pc.barrier()


    def prepare_models(self):

        if self.pc.id() == 0:
            # Setup the model on master node, this sets self.synapse_section_id (and _x)

            self.synapse_model = self.setup_model(synapse_density_override=None,
                                                  n_synapses_override=None,
                                                  synapse_position_override=None)




        self.pc.barrier()
        # Distribute section id, section x that was picked by master, and any other needed parameters
        self.synapse_section_id, self.synapse_section_x \
            = self.pc.broadcast((self.synapse_section_id, self.synapse_section_x))

        if self.pc.id() != 0:
            # Setup models on all other nodes (but not master)
            synapse_position_override = (self.synapse_section_id, self.synapse_section_x)

            self.synapse_model = self.setup_model(synapse_position_override=synapse_position_override)

        self.pc.barrier()

    def run_models(self, model_parameter_list):

        model_parameters = self.pc.py_scatter(model_parameter_list)



        # we need model parameters, and position of synapses (section_id, section_x)

        self.run_model()

        # TODO: 2026-03-05 WE ARE HERE, WORKING ON THIS FUNCTION!! SciLifeLab rulez!


    def optimise(self):

        self.setup_rng()

        # synapse_density_override
        # n_synapess_override
        # synapse_position_override





    def write_log(self, text, flush=True):  # Change flush to False in future, debug
        if self.log_file is not None:
            self.log_file.write(text + "\n")

            if self.verbose:
                print(text)

            if flush:
                self.log_file.flush()
        else:
            if self.verbose:
                print(text)


    def load_trace_data(self, data_file):

        self.write_log(f"Loading {data_file}")

        with open(data_file, "r") as f:
            self.data = json.load(f)

            self.volt = np.array(self.data["data"]["mean_norm_trace"]).flatten()
            self.sample_freq = self.data["meta_data"]["sample_frequency"]

            if "holding_voltage" in self.data["meta_data"]:
                self.trace_holding_voltage = self.data["meta_data"]["trace_holding_voltage"]
            else:
                self.trace_holding_voltage = None

            dt = 1 / self.sample_freq
            self.time = 0 + dt * np.arange(0, len(self.volt))

            self.stim_time = np.array(self.data["meta_data"]["stim_time"])

            self.cell_type = self.data["meta_data"]["cell_type"]

    def save_parameter_data(self):

        if self.pc.id() != 0:
            self.write_log("No servants are allowed to write output to json, ignoring call.")
            return

        self.write_log(f"Saving data to {self.parameter_data_file_name}")
        self.synapse_parameter_data.save(self.parameter_data_file_name)

    def load_parameter_data(self):

        # TODO: How should this be done with new architecture?

        self.synapse_parameter_data = ParameterBookkeeper(old_book_file=self.parameter_data_file_name, n_max=100)
        self.synapse_parameter_data.check_integrity()

        best_dataset = self.synapse_parameter_data.get_best_dataset()

        if best_dataset is not None:
            self.synapse_section_id = best_dataset["section_id"]
            self.synapse_section_x = best_dataset["section_x"]

        if self.pc.id() != 0:
            # This is to prevent duplicating entries
            self.synapse_parameter_data.clear()



    def get_cell_properties(self):

        if self.cell_properties is None:
            with open(self.neuron_set_file, 'r') as f:
                self.cell_properties = json.load(f)

        cell_type = self.data["meta_data"]["cell_type"]

        return copy.deepcopy(self.cell_properties[cell_type])


    def update_cell_properties(self, holding_current):

        cell_type = self.data["meta_data"]["cell_type"]

        with open(self.neuron_set_file, 'r') as f:
            self.cell_properties = json.load(f)

        self.cell_properties[cell_type]["holding_current"] = holding_current

        with open(self.neuron_set_file, 'w') as f:
            json.dump(self.cell_properties, f, indent=4)


    def setup_model(self,
                    synapse_density_override=None,
                    n_synapses_override=None,
                    synapse_position_override=None):

        self.write_log(f"setup_model: synapse_position-override: {synapse_position_override}")

        params = {}

        t_stim = self.stim_time

        # Read the info needed to setup the neuron hosting the synapses
        c_prop = self.get_cell_properties()

        if synapse_position_override is not None:
            synapse_section_id, synapse_section_x = synapse_position_override
        else:
            synapse_section_id, synapse_section_x = None, None

        if synapse_density_override is not None:
            synapse_density = synapse_density_override
        else:
            synapse_density = c_prop["synapse_density"]

        if n_synapses_override is not None:
            n_synapses = n_synapses_override
        else:
            if "num_synapses" in c_prop:
                n_synapses = c_prop["num_synapses"]
            elif 'nSynapses' in c_prop:
                n_synapses = c_prop['nSynapses']
            else:
                raise Exception('Setup error: number of synapses not no specified in setup file (which ever that is?)')

        if "holding_current" in c_prop:
            holding_current = c_prop["holding_current"]
        else:
            holding_current = None

        neuron_morphology_key = c_prop["neuron_morphology_key"]
        neuron_parameter_key = c_prop["neuron_parameter_key"]

        # Use the trace holding voltage if it exists, otherwise use the holding voltage in the neuronSet json file.
        if self.trace_holding_voltage is not None:
            trace_holding_voltage = self.trace_holding_voltage
        elif "baseline_voltage" in c_prop:
            trace_holding_voltage = c_prop["baseline_voltage"]
            self.trace_holding_voltage = trace_holding_voltage
        else:
            trace_holding_voltage = None
            assert f"You need to specify either a trace_holding_voltage in {self.data_file}" \
                   f"or specify baselineVoltage in neuronSet.json for the neuron type in question."

        # Temporarily force regeneration of holding current
        holding_current = None

        print(f"Using random seed {self.seed}")

        self.rsr_synapse_model = \
            RunSynapseRun(neuron_path=snudda_parse_path(c_prop["neuron_path"], self.snudda_data),
                          neuron_morphology_key=neuron_morphology_key,
                          neuron_parameter_key=neuron_parameter_key,
                          stim_times=t_stim,
                          num_synapses=n_synapses,
                          synapse_density=synapse_density,
                          holding_voltage=trace_holding_voltage,
                          holding_current=holding_current,
                          synapse_type=self.synapse_type,
                          params=params,
                          time=self.sim_time,
                          log_file=self.log_file,
                          synapse_section_id=synapse_section_id,
                          synapse_section_x=synapse_section_x,
                          random_seed=self.seed,
                          verbose=True)

        if self.rsr_synapse_model.holding_current != holding_current:
            self.update_cell_properties(holding_current=self.rsr_synapse_model.holding_current)

        return self.rsr_synapse_model


    def get_peak_idx(self, stim_time, time, volt):

        freq = 1.0 / (stim_time[1] - stim_time[0])

        p_window = 1.0 / (2 * freq) * np.ones(stim_time.shape)
        p_window[-1] *= 5

        peak_info = self.find_peaks_helper(p_time=stim_time,
                                           p_window=p_window,
                                           time=time,
                                           volt=volt)

        return peak_info["peakIdx"]


    # Find peaks within pStart[i] and pStart[i]+pWindow[i]
    # The value is not the amplitude of the peak, just the voltage at the peak

    def find_peaks_helper(self, p_time, p_window, time=None, volt=None):

        peak_idx = []
        peak_time = []
        peak_volt = []

        for pt, pw in zip(p_time, p_window):
            t_start = pt
            t_end = pt + pw

            t_idx = np.where(np.logical_and(t_start <= time, time <= t_end))[0]

            if len(t_idx) == 0:
                self.write_log(f"No time points within {t_start} and {t_end}", flush=True)
                import pdb
                pdb.set_trace()

            assert len(t_idx) > 0, f"No time points within {t_start} and {t_end}"

            if self.synapse_type == "glut":
                p_idx = t_idx[np.argmax(volt[t_idx])]
            elif self.synapse_type == "gaba":
                # We assume that neuron is more depolarised than -65, ie gaba is
                # also depolarising
                p_idx = t_idx[np.argmax(volt[t_idx])]
            else:
                self.write_log("Unknown synapse type : " + str(self.synapse_type), flush=True)
                import pdb
                pdb.set_trace()

            peak_idx.append(int(p_idx))
            peak_time.append(time[p_idx])
            peak_volt.append(volt[p_idx])

        # Save to cache -- obs peakVolt is NOT amplitude of peak, just volt

        peak_dict = {"peakIdx": np.array(peak_idx),
                    "peakTime": np.array(peak_time),
                    "peakVolt": np.array(peak_volt)}  # NOT AMPLITUDE

        return peak_dict

    def find_trace_heights(self, time, volt, peak_idx):

        decay_func = lambda x, a, b, c: a * np.exp(-x / b) + c

        v_base = np.mean(volt[int(0.3 * peak_idx[0]):int(0.8 * peak_idx[0])])

        peak_height = np.zeros((len(peak_idx)))
        peak_height[0] = volt[peak_idx[0]] - v_base

        decay_fits = []

        for idx_b in range(1, len(peak_idx)):

            if peak_height[0] > 0:
                if idx_b < len(peak_idx) - 1:
                    p0d = [0.06, 0.05, self.trace_holding_voltage]
                else:
                    p0d = [1e-5, 100, self.trace_holding_voltage]

                    if self.synapse_type == "gaba":
                        p0d = [1e-8, 10000, self.trace_holding_voltage]
            else:
                # In some cases for GABA we had really fast decay back
                if idx_b < len(peak_idx) - 1:
                    p0d = [-0.06, 0.05, self.trace_holding_voltage]
                else:
                    p0d = [-1e-5, 1e5, self.trace_holding_voltage]

            peak_idx_a = peak_idx[idx_b - 1]  # Prior peak
            peak_idx_b = peak_idx[idx_b]  # Next peak

            if idx_b < len(peak_idx) - 1:
                # Not the last spike
                idx_start = int(peak_idx_a * 0.9 + peak_idx_b * 0.1)
                idx_end = int(peak_idx_a * 0.1 + peak_idx_b * 0.9)
            else:
                # Last spike, use only last half of decay trace
                idx_start = int(peak_idx_a * 0.5 + peak_idx_b * 0.5)
                idx_end = int(peak_idx_a * 0.1 + peak_idx_b * 0.9)  # might need 0.85 as last

            try:
                assert idx_start < idx_end
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr, flush=True)

                import matplotlib.pyplot as plt

                plt.figure()
                plt.plot(time, volt)
                plt.xlabel("Time (error plot)")
                plt.ylabel("Volt (error plot)")
                plt.ion()
                plt.show()
                plt.title("ERROR!!!")
                import pdb
                pdb.set_trace()

            t_ab = time[idx_start:idx_end]
            v_ab = volt[idx_start:idx_end]

            t_ab_fit = t_ab - t_ab[0]
            v_ab_fit = v_ab

            try:

                try:
                    fit_params, pcov = scipy.optimize.curve_fit(decay_func, t_ab_fit, v_ab_fit, p0=p0d)
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    ### DEBUGGING START
                    import pickle
                    from datetime import datetime

                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    filename = f"curve_fit_args_{timestamp}.pkl"

                    dump_data = {
                        "t_ab_fit": t_ab_fit,
                        "v_ab_fit": v_ab_fit,
                        "p0d": p0d,
                    }

                    with open(filename, "wb") as f:
                        pickle.dump(dump_data, f)

                    print(f"Saved to {filename}")

                    ### DEBIGGING END

                    self.write_log("!!! Failed to converge, trying with smaller decay constant")
                    p0d[1] *= 0.01
                    fit_params, pcov = scipy.optimize.curve_fit(decay_func, t_ab_fit, v_ab_fit, p0=p0d)

                t_b = time[peak_idx_b] - t_ab[0]
                v_base_b = decay_func(t_b, fit_params[0], fit_params[1], fit_params[2])

                peak_height[idx_b] = volt[peak_idx_b] - v_base_b

                v_fit = decay_func(t_ab - t_ab[0], fit_params[0], fit_params[1], fit_params[2])
                decay_fits.append((t_ab, v_fit))

            except:
                self.write_log("Check that the threshold in the peak detection before is OK")
                # self.plot(name)
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr, flush=True)

                if True:
                    import matplotlib.pyplot as plt

                    plt.figure()
                    plt.plot(t_ab, v_ab, 'r')
                    plt.title("Error in findTraceHeights")
                    plt.xlabel("time")
                    plt.ylabel("volt")
                    # plt.plot(tAB,vFit,'k-')
                    plt.ion()
                    plt.show()

                import pdb
                pdb.set_trace()

        return peak_height.copy(), decay_fits, v_base

    def neuron_synapse_helper_glut(self, t_spike,
                                   u, tau_r, tau_f, tau_ratio, cond,
                                   smooth_exp_trace8, smooth_exp_trace9, exp_peak_height,
                                   return_type="peaks"):

        if self.debug_pars_flag:
            self.debug_pars.append([u, tau_r, tau_f, tau_ratio, cond])

        params = self.synapse_parameters
        tau = tau_r * tau_ratio

        peak_h, t_sim, v_sim = self.run_model(t_spike, u,
                                              tau_r, tau_f, cond, tau,
                                              params=params,
                                              return_trace=True)

        # Calculate error in decay fit
        sim_trace8, sim_time8 = self.smoothing_trace(v_sim, self.num_smoothing,
                                                     time=t_sim,
                                                     start_time=self.decay_start_fit8,
                                                     end_time=self.decay_end_fit8)

        sim_trace9, sim_time9 = self.smoothing_trace(v_sim, self.num_smoothing,
                                                     time=t_sim,
                                                     start_time=self.decay_start_fit9,
                                                     end_time=self.decay_end_fit9)

        # We only want to use the bit of the trace after max
        idx_max8 = np.argmax(smooth_exp_trace8)
        idx_max9 = np.argmax(smooth_exp_trace9)

        # Calculating error in peak height
        if self.normalise_trace:
            h_diff = np.abs(peak_h/peak_h[0] - exp_peak_height/exp_peak_height[0])
        else:
            h_diff = np.abs(peak_h - exp_peak_height)

        h_diff[0] *= 3
        h_diff[-2] *= 2
        h_diff[-1] *= 3

        # This is to prevent the model spiking
        spike_penalty = np.sum(peak_h > 0.03) * 20

        h_error = np.sum(h_diff) / len(h_diff)

        decay_error8 = np.mean((smooth_exp_trace8[idx_max8:] - sim_trace8[idx_max8:]) ** 2)
        decay_error9 = np.mean((smooth_exp_trace9[idx_max9:] - sim_trace9[idx_max9:]) ** 2)

        fit_error = h_error + decay_error8 + decay_error9 + spike_penalty

        print(f"{fit_error = }")

        if spike_penalty > 0:
            self.write_log("Action potential detected in trace. Penalising!")

        if False:
            peak_base = v_sim[-1]
            plt.figure()
            plt.plot(t_sim, v_sim, 'k-')
            plt.plot(sim_time8, sim_trace8, 'y--')
            plt.plot(sim_time8, smooth_exp_trace8, 'r--')
            plt.plot(sim_time9, sim_trace9, 'y--')
            plt.plot(sim_time9, smooth_exp_trace9, 'r--')

            for tp, expH, modH in zip(t_spike, exp_peak_height, peak_h):
                plt.plot([tp, tp], [peak_base, expH + peak_base], 'r-', linewidth=3)
                plt.plot([tp, tp], [peak_base, modH + peak_base], 'b-')
            plt.title("hE = %g, dE8 = %g, dE9 = %g" \
                      % (h_error, decay_error8, decay_error9))

            plt.ion()
            plt.show()

        if return_type == "peaks":
            return peak_h
        elif return_type == "error":
            return fit_error
        elif return_type == "full":
            return fit_error, peak_h, t_sim, v_sim
        else:
            assert False, "Unknown return type: " + str(return_type)

    def run_model(self, t_spike, u, tau_r, tau_f, cond, tau,
                  params=None,
                  return_trace=False):

        if params is None:
            params = {}

        # self.writeLog("Running neuron model")

        assert self.rsr_synapse_model is not None, \
            "!!! Need to call setupModel first"

        # Should we make a copy of params, to not destroy it? ;)
        params["U"] = u
        params["tauR"] = tau_r
        params["tauF"] = tau_f
        params["cond"] = cond
        params["tau"] = tau

        # self.writeLog("params=" + str(params))

        self.pc.barrier()
        (t_sim, v_sim, i_sim) = self.rsr_synapse_model.run2(pars=params)

        if t_sim.shape != v_sim.shape:
            self.write_log("Shape are different, why?!", flush=True)
            import pdb
            pdb.set_trace()

        peak_idx = self.get_peak_idx(time=t_sim, volt=v_sim, stim_time=t_spike)
        peak_height, decay_fits, v_base = self.find_trace_heights(t_sim, v_sim, peak_idx)

        if return_trace:
            return peak_height, t_sim, v_sim
        else:
            return peak_height

    # This should read from a JSON file instead

    def get_model_bounds(self):

        mb = self.data["model_data"]

        param_list = ["U", "tauR", "tauF", "tauRatio", "cond"]
        lower_bound = [mb[x][0] for x in param_list]
        upper_bound = [mb[x][1] for x in param_list]

        return lower_bound, upper_bound


    # TODO: BUtcher this function
    def parallel_optimise_single_cell(self, n_trials=10000, post_opt=False):

        start_time = timeit.default_timer()

        if self.pc.id() != 0:
            self.write_log("parallel_optimise_single_cell should only be called on master node")
            return

        # 1. Setup workers
        params = self.synapse_parameters

        if self.synapse_section_id is not None:
            syn_override = self.synapse_section_id, self.synapse_section_x
        else:
            syn_override = None

        # 2. Setup one cell to optimise, randomise synapse positions
        synapse_model = self.setup_model(params=params,
                                         synapse_position_override=syn_override)

        # (volt,time) = self.getData(dataType,cellID)
        peak_idx = self.get_peak_idx(stim_time=self.stim_time,
                                     time=self.time,
                                     volt=self.volt)
        t_spikes = self.time[peak_idx]

        sigma = np.ones(len(peak_idx))
        sigma[-1] = 1. / 3

        peak_height, decay_fits, v_base = self.find_trace_heights(self.time, self.volt, peak_idx)

        if False:
            plt.figure()
            plt.plot(self.time, self.volt, color="black")
            for idx in range(len(decay_fits)):
                plt.plot(decay_fits[idx][0], decay_fits[idx][1], color="red")
            plt.show()

        skip_sets = self.synapse_parameter_data.old_iter
        print(f"Starting sobol sequence at position {skip_sets}")

        # 2b. Create list of all parameter points to investigate
        model_bounds = self.get_model_bounds()
        parameter_points = self.setup_parameter_set(model_bounds, n_trials, skip_sets=skip_sets)

        # 3. Send synapse positions to all workers, and split parameter points
        #    between workers

        if self.d_view is not None:
            self.setup_parallel(self.d_view)

            self.d_view.scatter("parameter_points", parameter_points, block=True)

            self.d_view.push({"params": params,
                              "synapse_section_id": synapse_model.synapse_section_id,
                              "synapse_section_x": synapse_model.synapse_section_x,
                              "model_bounds": model_bounds,
                              "stim_time": self.stim_time,
                              "peak_height": peak_height},
                             block=True)

            cmd_str_setup = \
                "ly.sobol_worker_setup(params=params," \
                + "synapse_position_override=(synapse_section_id,synapse_section_x))"

            self.d_view.execute(f"ly.write_log('TESTING LOG A')", block=True)

            self.write_log("Calling sobol_worker_setup")
            self.d_view.execute(cmd_str_setup, block=True)

            self.d_view.execute(f"ly.write_log('TESTING LOG B')", block=True)

            cmd_str = ("res = ly.sobol_scan(t_stim = stim_time,"
                       "                    h_peak = peak_height,"
                       "                    parameter_sets=parameter_points,"
                       "                    model_bounds=model_bounds,"
                       "                    smooth_exp_trace8=ly.smooth_exp_volt8,"
                       "                    smooth_exp_trace9=ly.smooth_exp_volt9)")

            self.write_log("Executing workers, bang bang")
            self.d_view.execute(cmd_str, block=True)

            # 5. Gather worker data
            self.write_log("Gathering results from workers")
            # res = self.d_view["res"]
            res = self.d_view.gather("res", block=True)
            self.write_log("Results gathered.")

            #  for r in res:
            self.synapse_parameter_data.merge(res)

            self.save_parameter_data()

        else:

            # No dView, run in serial mode...
            self.sobol_worker_setup(params=params,
                                    synapse_position_override=(synapse_model.synapse_section_id,
                                                               synapse_model.synapse_section_x))

            self.sobol_scan(parameter_sets=parameter_points,
                            t_stim=self.stim_time,
                            h_peak=peak_height,
                            model_bounds=model_bounds,
                            smooth_exp_trace8=self.smooth_exp_volt8,
                            smooth_exp_trace9=self.smooth_exp_volt9)

            self.save_parameter_data()


        self.write_log(f"Sobol search done. Best parameter {self.synapse_parameter_data.get_best_parameterset()}")

        if post_opt:
            # This updates parameters and saves new parameter cache
            self.get_refined_parameters()

            self.save_parameter_data()

        end_time = timeit.default_timer()

        self.write_log(f"Optimisation duration: {end_time - start_time}.1f s", flush=True)



if __name__ == "__main__":
    so = SynapseOptimiser()


"""

from skopt import gp_minimize
from skopt import Optimizer
from joblib import Parallel, delayed


## Parallell optimisation

opt = Optimizer(dimensions=m_bounds, random_state=42)

n_iterations = 10
batch_size = 8  # tune this to your number of CPU cores

# TODO: This can not handle the self reference, need to make func
#       self contained.

for i in range(n_iterations):
    print(f"Iteration {i}/{n_iterations}")
    x_batch = opt.ask(n_points=batch_size)

    # Here run simulation, get results and put it in y_batch

    # y_batch = Parallel(n_jobs=-1)(delayed(func)(x) for x in x_batch)  # couldnt pickle
    opt.tell(x_batch, y_batch)

best_idx = opt.yi.index(min(opt.yi))
print("Best value:", min(opt.yi))
print("Best params:", opt.Xi[best_idx])

fit_params = opt.Xi[best_idx]
min_error = opt.yi
"""