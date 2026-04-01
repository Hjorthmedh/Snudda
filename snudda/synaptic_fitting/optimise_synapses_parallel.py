# TODO: 2026-02-26 -- Dear future humans.
#      -- We have copied optimise_synapse_full.py (but not all of it)
#      -- we need to write code to run on self.pc.id() == 0 that sets up the optimisation
#          and performs all the goodie stuff
#
#      We got to 'get_refined_parameters' function, probably not needed?
#      We are going to use mpirun and neuron parallel, to run multiple instances at the same exp_time.

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
import lzma

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
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel

from run_synapse_run import RunSynapseRun

from skopt import gp_minimize
from skopt import Optimizer
from joblib import Parallel, delayed

# TODO: If we want to use tmGlut_double we have to make sure all variables are initialised
#       currently several of the variables are not declared, and thus 0...
#       factor_ampa, tau1_ampa, I2_ampa, nmda_ratio, ... etc

# TODO: 2026-03-31
#       X 1. Synapse conductance needs to be set to same values -- conductance was set by optimizer
#       2. Fix error function (weighting and include decay?)
#       X 3. Plot results, and verify it looks ok
#       4. Run on Dardel
#       5. tmGlutDouble
#       6. Celebrate!
#

class SynapseOptimiser:

    def __init__(self, data_file,
                 entropy=1023456734529028340264793840,
                 synapse_type="glut",  # Change to "glut2" for tmGlut_double
                 load_parameters=True,
                 snudda_data=None,
                 neuron_set_file="neuron_set.json",
                 synapse_parameter_file=None,
                 verbose=True):

        self.log_file = None
        self.verbose = verbose
        self.rng = None
        self.sim = None

        self.data_file = data_file
        self.parameter_data_file_name = f"{self.data_file}-parameters-optimised.json"
        self.opt_state_data_file_name = f"{self.data_file}-opt-state.json.xz"   # Maybe change format if files get too big...

        self.neuron_set_file = neuron_set_file
        self.seed = None
        self.entropy=entropy

        self.snudda_data = get_snudda_data(snudda_data=snudda_data)

        self.rsr_synapse_model = None
        self.synapse_type = synapse_type
        self.data = None
        self.exp_volt = None
        self.exp_time = None
        self.exp_volt_interpolated = {}
        self.sample_freq = None
        self.sim_time = 1.5
        self.trace_holding_voltage = None
        self.stim_time = None
        self.exp_peak_height = None
        self.cell_type = None

        self.last_run_volt = None
        self.last_run_time = None

        self.synapse_parameter_data = None
        self.synapse_section_id = None
        self.synapse_section_x = None

        self.cell_properties = None

        self.pc = h.ParallelContext()
        self.n_workers = self.pc.nhost()

        self.load_trace_data()

        self.load_parameters = load_parameters

        if load_parameters:
            self.load_parameter_data()

        self.synapse_parameter_file = synapse_parameter_file

        if synapse_parameter_file:
            with open(synapse_parameter_file, 'r') as f:
                self.write_log(f"Reading synapse parameters from {synapse_parameter_file}")
                self.synapse_parameters = json.load(f)["data"]
        else:
            self.synapse_parameters = {}

        self.setup_rng()

    def setup_rng(self):

        if self.rng is not None:
            print(f"setup_rng: rng already setup, skipping.")
            return

        seeds = []

        if self.pc.id() == 0:
            # Setup and distribute random seeds to all workers
            seed_sequence = np.random.SeedSequence(entropy=self.entropy)
            seeds=list(seed_sequence.generate_state(self.n_workers))

        self.seed = self.pc.py_scatter(seeds)
        self.rng = np.random.default_rng(seed=self.seed)

        print(f"Worker: {self.pc.id()} -- seed: {self.seed}")

        self.pc.barrier()


    # This function should only be called ones

    def prepare_models(self):

        if self.sim is None:
            self.sim = NrnSimulatorParallel(cvode_active=False)


        # This sets self.rsr_synapse_model
        print(f"Worker {self.pc.id()} calling setup_model")
        self.setup_model(synapse_density_override=None,
                         n_synapses_override=None,
                         synapse_position_override=(self.synapse_section_id, self.synapse_section_x),
                         init_synapses=self.pc.id() == 0)

        if self.pc.id() == 0:
            # Get synapse id and x from master node
            self.synapse_section_id = self.rsr_synapse_model.synapse_section_id
            self.synapse_section_x = self.rsr_synapse_model.synapse_section_x
            # TODO: Do we need to get the conductances also?!! !!!

        self.pc.barrier()
        # Distribute section id, section x that was picked by master, and any other needed parameters
        self.synapse_section_id, self.synapse_section_x \
            = self.pc.py_broadcast((self.synapse_section_id, self.synapse_section_x), 0)

        if self.pc.id() != 0:
            # Setup models on all other nodes (but not master)

            print(f"Worker {self.pc.id()} adding master nodes synapses.")
            self.rsr_synapse_model.setup_synapses(synapse_type=self.synapse_type,
                                                  num_synapses=len(self.synapse_section_id),
                                                  synapse_section_id=self.synapse_section_id,
                                                  synapse_section_x=self.synapse_section_x)

        self.pc.barrier()

    def run_models(self, model_parameter_list):

        # model_parameter_list should be the parameters for the master, and empty [] for the workers
        # master then distributes the parameters to the workers.

        # prepare_models should already be called, so that synapse position is fixed apriori

        model_parameters = self.pc.py_scatter(model_parameter_list, 0)

        # we need model parameters, and position of synapses (section_id, section_x)

        if len(model_parameters) != 5:
            raise ValueError(f"There should be five model parameters: {model_parameters}")

        m_params = { "U": model_parameters[0],
                     "tauR": model_parameters[1],
                     "tauF": model_parameters[2],
                     "tauRatio": model_parameters[3],
                     "cond": model_parameters[4] }

        t_sim, v_sim, i_sim = self.rsr_synapse_model.run2(pars=m_params)

        self.last_run_time = t_sim
        self.last_run_volt = v_sim

        # We use normalised voltage instead of v_sim
        v_norm = (v_sim - np.min(v_sim)) / (np.max(v_sim) - np.min(v_sim))

        peak_idx = self.get_peak_idx(time=t_sim, volt=v_norm, stim_time=self.stim_time)
        peak_height, decay_fits, v_base = self.find_trace_heights(t_sim, v_norm, peak_idx)

        # We need to take decay into accounts also for error, first version only uses peak heights
        error = self.error_calculation(peak_height=peak_height,
                                       decay_fits=decay_fits,
                                       time=t_sim,
                                       volt=v_norm,
                                       v_base=v_base)

        error = self.pc.py_gather(error, 0)

        return error

        # TODO: 2026-03-05 WE ARE HERE, WORKING ON THIS FUNCTION!! SciLifeLab rulez!


    def error_calculation(self, peak_height, decay_fits, time, volt, v_base):

        decay_window = [0.01, 0.045]

        try:
            # Error in peak heights
            peak_error = np.abs(peak_height - self.exp_peak_height)

            # Weight errors
            peak_error[0] *= 3
            peak_error[-2] *= 2
            peak_error[-1] *= 3

            decay_error = 0

            # Error in decay
            for st in self.stim_time:
                start_idx = np.argmin(np.abs(time - (st + decay_window[0])))
                end_idx = np.argmin(np.abs(time - (st + decay_window[1])))

                if st not in self.exp_volt_interpolated:
                    self.exp_volt_interpolated[st] = np.interp(time[start_idx:end_idx],
                                                               self.exp_time,
                                                               self.exp_volt)

                # TODO, interpolate points to match exp data!!

                decay_error += np.sum(np.abs(volt[start_idx:end_idx] - self.exp_volt_interpolated[st])) / (end_idx - start_idx)

            print(f"Peak error: {np.sum(peak_error)}, decay error: {decay_error}")

            error = np.sum(peak_error) + decay_error

        except Exception as e:
            import traceback
            print(traceback.format_exc())
            print(e)
            import pdb
            pdb.set_trace()

        return error



    def error_calculation_peaks_only(self, peak_height, decay_fits, v_base):

        try:
            peak_error = np.sum(np.abs(peak_height - self.exp_peak_height))
        except Exception as e:
            import traceback
            print(traceback.format_exc())
            print(e)
            import pdb
            pdb.set_trace()

        return peak_error

    def load_opt_state(self, opt):

        if self.pc.id() != 0:
            return

        if os.path.isfile(self.opt_state_data_file_name):
            print(f"Loading optmisation state from {self.opt_state_data_file_name}")
            with lzma.open(self.opt_state_data_file_name, "rt") as f:
                state = json.load(f)

            print(f"Found {len(state['yi'])} previous data points.")

            # Instruct the optimizer about previous evaluations
            opt.tell(state["xi"], state["yi"])

    def save_opt_state(self, opt):

        if self.pc.id() != 0:
            return

        state = { "xi": opt.Xi,
                  "yi": opt.yi }

        print(f"Saving optmisation state to {self.opt_state_data_file_name}")
        with lzma.open(self.opt_state_data_file_name, "wt") as f:
            json.dump(state, f, indent=4)



    def optimise(self, n_iterations=10, load_state=True):

        error_list = []
        start_time = time.perf_counter()

        if self.seed is None:
            self.setup_rng()

        self.prepare_models()

        if self.pc.id() == 0:
            model_bounds = self.get_model_bounds()
            model_bounds = [x for x in zip(*model_bounds)]
            opt = Optimizer(dimensions=model_bounds, random_state=42)

            if self.load_parameters:
                self.load_opt_state(opt)

        for i in range(n_iterations):

            if self.pc.id() == 0:
                print(f"Iteration {i}/{n_iterations}")
                model_parameter_list = opt.ask(n_points=self.n_workers)
                # TODO: Should we round model_parameter_list to N decimals before proceeding?
            else:
                model_parameter_list = []

            error = self.run_models(model_parameter_list)

            if self.pc.id() == 0:
                opt.tell(model_parameter_list, error)
                print(f"Error: {error}")

                if i % 100 == 0 and i > 0:
                    # Just for safety let's save every 100 iterations...
                    print(f"Iteration {i}: Saving state to {self.opt_state_data_file_name}")
                    self.save_opt_state(opt)

                error_list.append(np.min(opt.yi))

        if self.pc.id() == 0:
            best_idx = opt.yi.index(min(opt.yi))
            print("Best value:", opt.yi[best_idx])
            print("Best params:", opt.Xi[best_idx])
            fit_params = opt.Xi[best_idx]
            min_error = opt.yi[best_idx]

            self.synapse_parameter_data.add_parameters(parameter_set=fit_params,
                                                       section_id=self.rsr_synapse_model.synapse_section_id,
                                                       section_x=self.rsr_synapse_model.synapse_section_x,
                                                       error=min_error)

            self.save_parameter_data()
            self.save_opt_state(opt)

        self.pc.barrier()

        # This reruns the best run, then plots it
        self.run_best_run()
        self.plot_last_run()
        self.plot_error(error_list)

        duration = time.perf_counter() - start_time
        print(f"Duration: {duration} seconds")

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


    def load_trace_data(self, data_file=None):

        if data_file is None:
            data_file = self.data_file
        else:
            self.data_file = data_file

        self.write_log(f"Loading {data_file}")

        with open(data_file, "r") as f:
            self.data = json.load(f)

        self.exp_volt = np.array(self.data["data"]["mean_norm_trace"]).flatten()
        self.sample_freq = self.data["meta_data"]["sample_frequency"]

        if "holding_voltage" in self.data["meta_data"]:
            self.trace_holding_voltage = self.data["meta_data"]["holding_voltage"]
        else:
            self.trace_holding_voltage = np.mean(self.data["data"]["mean_norm_trace"][:10])
            print(f"Guessing holding voltage: {self.trace_holding_voltage}")

        if self.trace_holding_voltage > 0:
            raise ValueError(f"Your holding voltage is probably wrong: {self.trace_holding_voltage} V")

        dt = 1 / self.sample_freq
        self.exp_time = 0 + dt * np.arange(0, len(self.exp_volt))

        self.stim_time = np.array(self.data["meta_data"]["stim_time"])

        self.cell_type = self.data["meta_data"]["cell_type"]

        peak_idx = self.get_peak_idx(time=self.exp_time, volt=self.exp_volt, stim_time=self.stim_time)

        self.exp_peak_height, _, _ = self.find_trace_heights(time=self.exp_time,
                                                             volt=self.exp_volt,
                                                             peak_idx=peak_idx)

    def save_parameter_data(self):

        if self.pc.id() != 0:
            self.write_log("No servants are allowed to write output to json, ignoring call.")
            return

        self.write_log(f"Saving data to {self.parameter_data_file_name}")
        self.synapse_parameter_data.save(self.parameter_data_file_name)

    def load_parameter_data(self):

        if self.pc.id() != 0:
            return

        print(f"Loading parameters from {self.parameter_data_file_name}")

        self.synapse_parameter_data = ParameterBookkeeper(old_book_file=self.parameter_data_file_name, n_max=100)
        self.synapse_parameter_data.check_integrity()

        best_dataset = self.synapse_parameter_data.get_best_dataset()

        if best_dataset is not None:
            self.synapse_section_id = best_dataset["section_id"]
            self.synapse_section_x = best_dataset["section_x"]


    def get_cell_properties(self):

        if self.cell_properties is None:
            with open(self.neuron_set_file, 'r') as f:
                self.cell_properties = json.load(f)

        cell_type = self.data["meta_data"]["cell_type"]

        return copy.deepcopy(self.cell_properties[cell_type])


    def update_cell_properties(self, holding_current):

        if self.pc.id() != 0:
            return

        cell_type = self.data["meta_data"]["cell_type"]

        with open(self.neuron_set_file, 'r') as f:
            self.cell_properties = json.load(f)

        self.cell_properties[cell_type]["holding_current"] = holding_current

        with open(self.neuron_set_file, 'w') as f:
            json.dump(self.cell_properties, f, indent=4)


    def setup_model(self,
                    synapse_density_override=None,
                    n_synapses_override=None,
                    synapse_position_override=None,
                    init_synapses=True):

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

        print(f"t_stim = {t_stim}")

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
                          sim=self.sim,
                          random_seed=self.seed,
                          init_synapses=init_synapses,
                          verbose=True,
                          pc=self.pc)

        self.pc.barrier()

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
                self.write_log(f"No exp_time points within {t_start} and {t_end}", flush=True)
                import pdb
                pdb.set_trace()

            assert len(t_idx) > 0, f"No exp_time points within {t_start} and {t_end}"

            if self.synapse_type in ("glut", "glut2"):
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

        # Save to cache -- obs peakVolt is NOT amplitude of peak, just exp_volt

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
                    plt.xlabel("exp_time")
                    plt.ylabel("exp_volt")
                    # plt.plot(tAB,vFit,'k-')
                    plt.ion()
                    plt.show()

                import pdb
                pdb.set_trace()

        return peak_height.copy(), decay_fits, v_base

    # This should read from a JSON file instead

    def get_model_bounds(self):

        mb = self.data["model_data"]

        param_list = ["U", "tauR", "tauF", "tauRatio", "cond"]
        lower_bound = [mb[x][0] for x in param_list]
        upper_bound = [mb[x][1] for x in param_list]

        return lower_bound, upper_bound

    def plot_last_run(self, fig_name=None):

        if self.pc.id() != 0:
            return

        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(self.last_run_time, (self.last_run_volt  - np.min(self.last_run_volt))/ (np.max(self.last_run_volt) - np.min(self.last_run_volt)), color='black', label="model")
        plt.plot(self.exp_time, self.exp_volt , color='red', label="experiment")
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage")
        plt.legend()

        if fig_name is None:
            fig_name = os.path.join("figures", os.path.basename(self.data_file).split(".")[0] + ".png")

        os.makedirs("figures", exist_ok=True)

        plt.savefig(fig_name, dpi=300)

    def run_best_run(self):

        # Get the best parameters
        # Distribute the best parameters to workers (all will be identical), wasteful (run just one worker)
        # Run

        if self.pc.id() == 0:
            best_param = self.synapse_parameter_data.get_best_parameterset()
            best_param_list = [best_param for x in range(self.n_workers)]
        else:
            best_param_list = []

        self.run_models(best_param_list)

    def plot_error(self, error_list):
        if self.pc.id() != 0:
            return

        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(error_list)
        plt.ylabel("Error")

        fig_name = os.path.join("figures", os.path.basename(self.data_file).split(".")[0] + "-error.png")

        plt.savefig(fig_name, dpi=300)
        plt.close()



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Synapse optimisation using MPI parallel workers.")
    parser.add_argument("data_file", type=str,
                        help="Path to the data file (JSON) to optimise against.")
    parser.add_argument("--iterations", type=int, default=100,
                        help="Number of optimisation iterations to run (default: 100).")
    parser.add_argument("--snudda_data", type=str, default=None,
                        help="Path to the Snudda data directory.")
    args = parser.parse_args()

    so = SynapseOptimiser(data_file=args.data_file,
                          snudda_data=args.snudda_data)
    so.optimise(n_iterations=args.iterations)

    # mpirun -n 5 python optimise_synapses_parallel.py ../data/synapses/example_data/10_MSN12_GBZ_CC_H20.json --iterations 50 --snudda_data /home/hjorth/HBP/BasalGangliaData/data/

    # Stored:  /media/psf/KTH/2025-11-13-Yvonne-data-teanalysing/Yvonne2019/CategorisedSTP/
    # mpirun -n 5 python optimise_synapses_parallel.py ../data/synapses/example_data/mixed_test_data.json  --iterations 3 --snudda_data /home/hjorth/HBP/BasalGangliaData/data/

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