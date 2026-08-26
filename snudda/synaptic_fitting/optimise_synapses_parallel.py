#
# TO RUN THIS CODE YOU ALSO NEED:
# pip install scikit-optimize

# TODO: 2026-02-26 -- Dear future humans.
#      -- We have copied optimise_synapse_full.py (but not all of it)
#      -- we need to write code to run on self.pc.id() == 0 that sets up the optimisation
#          and performs all the goodie stuff
#
# TODO: 2026-04-09: tmglut_double (glut2) is running on Dardel, hopefully the result will look good.
#       We need to get the data for all the synapses, for optimisation.
#       Then run and verify.
#       We need a way to setup the runs for multiple different neurons of the same type.
#       So we want to be able to specify neuron-set file, and random seed.
#       * We need to make a valiation step, where we randomly pick neuron (morph and parameter keys) and check
#       that it looks good.
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
from scipy.signal import find_peaks

import datetime

from mpi4py import MPI  # This must be imported before neuron, to run parallel
import neuron
from neuron import h  # , gui

from snudda.utils.snudda_path import snudda_parse_path, get_snudda_data
from snudda.synaptic_fitting.parameter_bookkeeper import ParameterBookkeeper
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel
from snudda.utils.debug import can_debug

from run_synapse_run import RunSynapseRun

from skopt import gp_minimize
from skopt import Optimizer
from skopt.acquisition import _gaussian_acquisition

from joblib import Parallel, delayed
from skopt.learning import RandomForestRegressor as SkoptRandomForestRegressor



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

def ask_batch_fast(opt, n_points, candidate_pool_size=2000, min_dist_frac=0.03):
    """Replacement for opt.ask(n_points=..., strategy='cl_min').

    skopt's built-in batch ask() copies the optimizer and refits the
    surrogate model n_points+1 times (once per fake 'lie' point told to
    the copy). This function instead reuses the model already fit by the
    last opt.tell() call and ranks a pool of random candidates by
    acquisition value -- one model, zero extra refits.
    """
    if not opt.models:
        # Still in the initial random-sampling phase, no surrogate to reuse yet.
        return opt.ask(n_points=n_points)

    est = opt.models[-1]

    X_pool = opt.space.transform(
        opt.space.rvs(n_samples=candidate_pool_size, random_state=opt.rng)
    )

    values = _gaussian_acquisition(
        X=X_pool,
        model=est,
        y_opt=np.min(opt.yi),
        acq_func="EI",
        acq_func_kwargs=opt.acq_func_kwargs,
    )

    # Normalise each dimension to [0, 1] just for the diversity check below --
    # your params span wildly different scales (cond ~1e-9, U ~0-1), so raw
    # distance would be dominated by whichever parameter has the biggest units.
    bounds = np.array(opt.space.transformed_bounds)
    ranges = np.maximum(bounds[:, 1] - bounds[:, 0], 1e-12)
    X_norm = (X_pool - bounds[:, 0]) / ranges

    order = np.argsort(values)  # best (lowest) acquisition value first
    chosen, chosen_norm = [], []

    for idx in order:
        x_n = X_norm[idx]
        if chosen_norm:
            d = np.linalg.norm(np.array(chosen_norm) - x_n, axis=1).min()
            if d < min_dist_frac:
                continue  # too close to an already-picked point this batch, skip
        chosen.append(idx)
        chosen_norm.append(x_n)
        if len(chosen) == n_points:
            break

    if len(chosen) < n_points:
        # Diversity filter was too strict to fill the batch -- top up with
        # the next-best points regardless of spacing.
        for idx in order:
            if idx not in chosen:
                chosen.append(idx)
            if len(chosen) == n_points:
                break

    return opt.space.inverse_transform(X_pool[chosen])

class SynapseOptimiser:

    def __init__(self, data_file,
                 entropy=1023456734529028340264793840,
                 synapse_type="glut",  # Change to "glut2" for tmGlut_double
                 load_parameters=True,
                 snudda_data=None,
                 neuron_set_file="neuron_set.json",
                 n_synapses = None,
                 n_soma_synapses = None,
                 synapse_density = None,
                 name_tag=None,
                 synapse_parameter_file=None,
                 update_neuronset_file=True,
                 verbose=True):

        os.makedirs("log", exist_ok=True)
        self.pc = h.ParallelContext()

        self.name_tag = name_tag
        slurm_job_id = os.environ.get("SLURM_JOB_ID", os.environ.get("SLURM_JOBID"))

        if slurm_job_id is not None:
            log_file_name = f"log/opt_log_job_{slurm_job_id}_rank-{self.pc.id()}.txt"
        else:
            log_file_name = f"log/opt_log_rank-{self.pc.id()}.txt"

        if self.name_tag is not None:
            log_file_name = log_file_name.replace(".txt", f"-{self.name_tag}.txt")

        self.log_file = open(log_file_name, "w")
        print(f"Writing to log file: {log_file_name}")

        self.verbose = verbose
        self.rng = None
        self.sim = None

        # This is useful if running multiple optimisations at the same time
        # the first run will set this, subsequent runs will cause trouble if they run in parallel
        # and all try to update it at same time.
        self.update_neuronset_file = update_neuronset_file

        self.parameter_list = ["U", "tauR", "tauF", "tauRatio", "nmda_ratio"]

        self.data_file = data_file
        new_path, data_file_name = os.path.split(data_file)
        new_path = os.path.join(new_path, "fitted")
        os.makedirs(new_path, exist_ok=True)
        new_file_name = data_file_name.replace(".json", "")

        if self.name_tag is not None:
            new_file_name += f"-{self.name_tag}"

        self.parameter_data_file_name = os.path.join(new_path, f"{new_file_name}-parameters-optimised-{synapse_type}.json")
        self.opt_state_data_file_name = os.path.join(new_path, f"{new_file_name}-opt-state-{synapse_type}.json.xz")   # Maybe change format if files get too big...
        self.parameter_export_file_name = os.path.join(new_path, f"{new_file_name}-synapse-parameters-{synapse_type}.json")

        self.neuron_set_file = neuron_set_file
        self.seed = None
        self.entropy=entropy

        # These are by default None, normally read from neuron_set file,
        # but we have option to override them when calling
        self.n_synapses_override = n_synapses
        self.synapse_density_override = synapse_density
        self.n_soma_synapses_override = n_soma_synapses

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
        self.last_run_parameters = None
        self.syn_recording = None

        self.synapse_parameter_data = None
        self.synapse_section_id = None
        self.synapse_section_x = None

        self.cell_properties = None

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
            if self.verbose:
                self.write_log(f"setup_rng: rng already setup, skipping.")
            return

        seeds = []

        if self.pc.id() == 0:
            # Setup and distribute random seeds to all workers
            seed_sequence = np.random.SeedSequence(entropy=self.entropy)
            seeds=list(seed_sequence.generate_state(self.n_workers))

        self.seed = self.pc.py_scatter(seeds)
        self.rng = np.random.default_rng(seed=self.seed)

        if self.verbose:
            self.write_log(f"Worker: {self.pc.id()} -- seed: {self.seed}")

        self.pc.barrier()


    # This function should only be called ones

    def prepare_models(self):

        if self.sim is None:
            self.sim = NrnSimulatorParallel(cvode_active=False)

        if self.verbose:
            self.write_log(f"synapse_parameters = {self.synapse_parameters}")
            self.write_log(f"Worker {self.pc.id()} calling setup_model")

        # This sets self.rsr_synapse_model
        self.setup_model(synapse_density_override=self.synapse_density_override,
                         n_synapses_override=self.n_synapses_override,
                         n_soma_synapses_override=self.n_soma_synapses_override,
                         synapse_params=self.synapse_parameters,
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

            if self.verbose:
                self.write_log(f"Worker {self.pc.id()} adding master nodes synapses.")

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

        # print(f"Worker {self.pc.id()} received: {model_parameters}")

        try:
            m_params = dict(zip(self.parameter_list, model_parameters))

            t_sim, v_sim, i_sim = self.rsr_synapse_model.run2(pars=m_params)

            self.last_run_time = t_sim
            self.last_run_volt = v_sim

            # We use normalised voltage instead of v_sim
            v_norm = (v_sim - np.min(v_sim)) / (np.max(v_sim) - np.min(v_sim))

            max_volt = np.max(v_sim)

            peak_idx = self.get_peak_idx(time=t_sim, volt=v_norm, stim_time=self.stim_time)
            peak_height, decay_fits, v_base = self.find_trace_heights(t_sim, v_norm, peak_idx)

            # We need to take decay into accounts also for error, first version only uses peak heights
            error = self.error_calculation(peak_height=peak_height,
                                           decay_fits=decay_fits,
                                           time=t_sim,
                                           volt=v_norm,
                                           v_base=v_base,
                                           max_volt=max_volt)

            self.last_run_parameters = m_params

        except Exception:
            import traceback
            error_str = traceback.format_exc()
            self.write_log(f"Error during model evaluation, and error calculation: {error_str}")
            # We set a high error to mark this as bad.
            error = 1e4

        if self.verbose:
            self.write_log(f"Worker {self.pc.id()} error: {error}")

        error = self.pc.py_gather(error, 0)

        return error

        # TODO: 2026-03-05 WE ARE HERE, WORKING ON THIS FUNCTION!! SciLifeLab rulez!


    def error_calculation(self, peak_height, decay_fits, time, volt, v_base, max_volt):

        decay_windows = [(0.01, 0.045)] * (len(self.stim_time) - 2) + [(0.01, 0.21)] * 2

        try:
            # Error in peak heights
            peak_error = np.abs(peak_height - self.exp_peak_height)

            # Weight errors
            peak_error[0] *= 3
            peak_error[-2] *= 2
            peak_error[-1] *= 3

            decay_error = 0

            # Error in decay
            for st, decay_window  in zip(self.stim_time, decay_windows):
                start_idx = np.argmin(np.abs(time - (st + decay_window[0])))
                end_idx = np.argmin(np.abs(time - (st + decay_window[1])))

                # Interpolate points, since the sampling rates might differ between exp and  simulation
                if st not in self.exp_volt_interpolated:
                    self.exp_volt_interpolated[st] = np.interp(time[start_idx:end_idx],
                                                               self.exp_time,
                                                               self.exp_volt)

                decay_error += np.sum(np.abs(volt[start_idx:end_idx] - self.exp_volt_interpolated[st])) / (end_idx - start_idx)

            peak_data, peak_prop = find_peaks(volt, threshold=np.median(volt))
            n_peaks = len(peak_data)
            n_peak_error = np.abs(9 - n_peaks) * 10

            spike_threshold = -20e-3

            if max_volt > spike_threshold:
                max_volt_error = (max_volt - spike_threshold) * 1e4
            else:
                max_volt_error = 0

            if self.verbose:
                self.write_log(f"Peak error: {np.sum(peak_error)}, decay error: {decay_error}, num peak error: {n_peak_error}, max_volt_error: {max_volt_error}")

            error = np.sum(peak_error) + decay_error + n_peak_error + max_volt_error

        except Exception as e:
            import traceback
            t_str = traceback.format_exc()
            self.write_log(f"error_calculation: {t_str}")
            print(e)

            if can_debug():
                import pdb
                pdb.set_trace()

            # Raise the error onwards, outer loop has to handle it.
            raise e

        return error



    def error_calculation_peaks_only(self, peak_height, decay_fits, v_base):

        try:
            peak_error = np.sum(np.abs(peak_height - self.exp_peak_height))
        except Exception as e:
            import traceback
            self.write_log(traceback.format_exc())
            print(e)

            if can_debug():
                import pdb
                pdb.set_trace()

            raise e

        return peak_error

    def load_opt_state(self, opt):

        if self.pc.id() != 0:
            return

        if os.path.isfile(self.opt_state_data_file_name):
            if self.verbose:
                self.write_log(f"Loading optmisation state from {self.opt_state_data_file_name}")
            with lzma.open(self.opt_state_data_file_name, "rt") as f:
                state = json.load(f)

            if self.verbose:
                self.write_log(f"Found {len(state['yi'])} previous data points.")

            # Instruct the optimizer about previous evaluations
            opt.tell(state["xi"], state["yi"])

    def save_opt_state(self, opt):

        if self.pc.id() != 0:
            return

        if self.verbose:
            self.write_log(f"Saving opt state: {len(opt.Xi)} xi points, {len(opt.yi)} yi points")
            # self.write_log(f"yi = {opt.yi}")

        state = { "xi": opt.Xi,
                  "yi": opt.yi }

        self.write_log(f"Saving optmisation state to {self.opt_state_data_file_name}")
        with lzma.open(self.opt_state_data_file_name, "wt") as f:
            json.dump(state, f)



    def optimise(self, n_iterations=10, load_state=True):

        if n_iterations < 1:
            return

        error_list = []
        start_time = time.perf_counter()

        if self.seed is None:
            self.setup_rng()

        self.prepare_models()

        if self.pc.id() == 0:
            model_bounds = self.get_model_bounds()
            model_bounds = [x for x in zip(*model_bounds)]

            # base_estimator = SkoptRandomForestRegressor(n_estimators=20, n_jobs=-1, random_state=42)
            # opt = Optimizer(dimensions=model_bounds, random_state=42, base_estimator=base_estimator)

            # 2026-07-31, speedup suggested by Claude
            base_estimator = SkoptRandomForestRegressor(n_estimators=10, n_jobs=1, random_state=42)
            opt = Optimizer(dimensions=model_bounds, random_state=42, base_estimator=base_estimator,
                            acq_func="EI",
                            acq_optimizer_kwargs={"n_points": 2000})

            # opt = Optimizer(dimensions=model_bounds, random_state=42, base_estimator="RF")

            if self.load_parameters:
                self.load_opt_state(opt)

        for i in range(n_iterations):

            if self.pc.id() == 0:
                if self.verbose:
                    self.write_log(f"Iteration {i}/{n_iterations}")

                # model_parameter_list = opt.ask(n_points=self.n_workers)

                # 2026-07-31, speedup suggsted by Claude. Optimisation was running really slowly on big cluster
                model_parameter_list = ask_batch_fast(opt, n_points=self.n_workers)

            else:
                model_parameter_list = []

            error = self.run_models(model_parameter_list)

            if self.verbose:
                self.write_log(f"Worker {self.pc.id()} has neuron = {id(self.rsr_synapse_model.neuron)}")

            if self.pc.id() == 0:
                opt.tell(model_parameter_list, error)

                if i % 50 == 0 and i > 0:
                    # Just for safety let's save every 50 iterations...
                    elapsed_time = time.perf_counter() - start_time
                    self.write_log(f"Iteration {i}: Saving state to {self.opt_state_data_file_name} (elapsed time: {elapsed_time:.0f} seconds)")
                    self.save_opt_state(opt)

                error_list.append(np.min(opt.yi))

        if self.pc.id() == 0:
            best_idx = opt.yi.index(min(opt.yi))
            self.write_log(f"Best value: {opt.yi[best_idx]}")
            self.write_log(f"Best params: {opt.Xi[best_idx]}")
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

        if self.pc.id() == 0:
            self.plot_error(error_list)
            self.plot_error(opt.yi, fig_name_info="-ALL", linestyle="None")

        duration = time.perf_counter() - start_time
        self.write_log(f"Duration: {duration} seconds")

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

            if self.verbose:
                self.write_log(f"Guessing holding voltage: {self.trace_holding_voltage}")

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
            if self.verbose:
                self.write_log("No servants are allowed to write output to json, ignoring call.")
            return

        self.write_log(f"Saving data to {self.parameter_data_file_name}")
        self.synapse_parameter_data.save(self.parameter_data_file_name)

    def load_parameter_data(self):

        if self.pc.id() != 0:
            return

        self.write_log(f"Loading parameters from {self.parameter_data_file_name}")

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
                    n_soma_synapses_override=None,
                    synapse_position_override=None,
                    synapse_params=None,
                    init_synapses=True):

        self.write_log(f"setup_model: synapse_position-override: {synapse_position_override}")

        if synapse_params is None:
            synapse_params = {}

        t_stim = self.stim_time

        self.pc.barrier()

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

        if n_soma_synapses_override is not None:
            n_soma_synapses = n_soma_synapses_override
        else:
            n_soma_synapses = c_prop.get("num_soma_synapses", 0)

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
            raise ValueError(f"You need to specify either a trace_holding_voltage in {self.data_file} " \
                             f"or specify baselineVoltage in neuronSet.json for the neuron type in question.")

        # Temporarily force regeneration of holding current
        holding_current = None

        if self.verbose:
            self.write_log(f"Using random seed {self.seed}")
            self.write_log(f"t_stim = {t_stim}")

        self.rsr_synapse_model = \
            RunSynapseRun(neuron_path=snudda_parse_path(c_prop["neuron_path"], self.snudda_data),
                          neuron_morphology_key=neuron_morphology_key,
                          neuron_parameter_key=neuron_parameter_key,
                          stim_times=t_stim,
                          num_synapses=n_synapses,
                          num_soma_synapses=n_soma_synapses,
                          synapse_density=synapse_density,
                          holding_voltage=trace_holding_voltage,
                          holding_current=holding_current,
                          synapse_type=self.synapse_type,
                          params=synapse_params,
                          time=self.sim_time,
                          log_file=self.log_file,
                          synapse_section_id=synapse_section_id,
                          synapse_section_x=synapse_section_x,
                          sim=self.sim,
                          random_seed=self.seed,
                          init_synapses=init_synapses,
                          verbose=self.verbose,
                          pc=self.pc)

        self.pc.barrier()

        if self.rsr_synapse_model.holding_current != holding_current:
            self.log_file.write(f"Mismatch between passed holding current {holding_current}, and post-model detected holding current {self.rsr_synapse_model.holding_current}.")
            if self.update_neuronset_file:
                self.update_cell_properties(holding_current=self.rsr_synapse_model.holding_current)

        self.pc.barrier()

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

                if can_debug():
                    import pdb
                    pdb.set_trace()

                raise ValueError(f"No exp_time points within {t_start} and {t_end}")


            if self.synapse_type in ("glut", "glut2"):
                p_idx = t_idx[np.argmax(volt[t_idx])]
            elif self.synapse_type == "gaba":
                # We assume that neuron is more depolarised than -65, ie gaba is
                # also depolarising
                p_idx = t_idx[np.argmax(volt[t_idx])]
            else:
                self.write_log("Unknown synapse type : " + str(self.synapse_type), flush=True)
                if can_debug():
                    import pdb
                    pdb.set_trace()
                else:
                    raise ValueError(f"Unknown synapse type : {self.synapse_type}")

            peak_idx.append(int(p_idx))
            peak_time.append(time[p_idx])
            peak_volt.append(volt[p_idx])

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
            except Exception:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr, flush=True)

                import matplotlib.pyplot as plt

                plt.figure()
                plt.plot(time, volt)
                plt.xlabel("Time (error plot)")
                plt.ylabel("Volt (error plot)")
                plt.title(f"{idx_start =}, {idx_end =}")
                plt.ion()
                plt.savefig(f"error-plot-{datetime.datetime.now()}.png")
                plt.show()
                plt.title("ERROR!!!")
                raise ValueError(f"find_trace_heights: {idx_start =}, {idx_end =}")

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
                    if False:
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

            except Exception as e:
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
                    plt.savefig(f"error-plot-find-trace-heights{datetime.datetime.now()}.png")
                    plt.ion()
                    plt.show()

                if can_debug():
                    import pdb
                    pdb.set_trace()

                raise e

        return peak_height.copy(), decay_fits, v_base

    # This should read from a JSON file instead

    def get_model_bounds(self):

        mb = self.data["model_data"]

        try:
            lower_bound = [mb[x][0] for x in self.parameter_list]
            upper_bound = [mb[x][1] for x in self.parameter_list]
        except Exception as e:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr, flush=True)
            print(tstr)

            if can_debug():
                import pdb
                pdb.set_trace()
            else:
                raise ValueError(tstr)

        return lower_bound, upper_bound

    def plot_last_run(self, fig_name=None, normalise_volt=True):

        if self.pc.id() != 0:
            return

        import matplotlib.pyplot as plt

        plt.figure()

        # Skip artifacts
        t_idx = np.where(0.1 <= self.last_run_time)[0]

        if normalise_volt:
            # Skip the first 0.1s of the trace, sometimes artifacts, that we dont want to afffect normalisation
            volt = (self.last_run_volt  - np.min(self.last_run_volt[t_idx]))/ (np.max(self.last_run_volt[t_idx]) - np.min(self.last_run_volt[t_idx]))
        else:
            volt = self.last_run_volt


        plt.plot(self.last_run_time[t_idx], volt[t_idx], color='black', label="model")

        if normalise_volt:
            # Only plot experimental data if trace plotted was normalised, since data is normalised
            te_idx = np.where(0.1 <= self.exp_time)[0]

            plt.plot(self.exp_time[te_idx], self.exp_volt[te_idx] , color='red', label="experiment")

        plt.xlabel("Time (s)")
        plt.ylabel("Voltage")
        plt.legend()

        title_str = ", ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}"
                               for k, v in self.last_run_parameters.items())
        plt.title(title_str, fontsize=8, wrap=True)

        if fig_name is None:
            if self.name_tag is not None:
                name_tag = f"-{self.name_tag}"
            else:
                name_tag = ""

            fig_name = os.path.join("figures", os.path.basename(self.data_file).split(".")[0] + f"-{self.synapse_type}{name_tag}.png")

        os.makedirs("figures", exist_ok=True)

        plt.savefig(fig_name, dpi=300)

        if self.syn_recording is not None:
            plt.figure()
            for syn_record in self.syn_recording:
                plt.plot(syn_record)
            plt.title(title_str, fontsize=8, wrap=True)

            fig_name_cur = fig_name.replace(".png", "-current.png")
            plt.savefig(fig_name_cur, dpi=300)
            print(f"Plotting current to {fig_name_cur}")


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

    def export_best_parameters(self, output_file=None, parameter_key="1", overwrite=False):

        if self.pc.id() == 0:
            if output_file is None:
                output_file = self.parameter_export_file_name

            best_param = self.synapse_parameter_data.get_best_parameterset()

            if not overwrite and os.path.isfile(output_file):
                with open(output_file, "r") as f:
                    data = json.load(f)
            else:
                data = dict()

            # In case any extra parameters were specified
            syn_param = self.synapse_parameters.copy()
            syn_param |= { k:v for k,v in zip(self.parameter_list, best_param) }

            # Because tauR > tau we use tauRatio during optimisation, convert back to tau in parameter file
            # since MOD file does not know about tauRatio
            if "tauRatio" in syn_param:
                if "tau" in syn_param:
                    raise ValueError(f"tau = tauR * tauRatio, however tau already set in {syn_param = }")

                syn_param["tau"] = syn_param["tauR"] * syn_param["tauRatio"]
                del syn_param["tauRatio"]

            data[parameter_key] = {
                "meta": {"parameter_data_file": os.path.basename(self.parameter_data_file_name)},
                "synapse":  syn_param
            }

            os.makedirs(os.path.dirname(output_file), exist_ok=True)

            with open(output_file, "w") as f:
                json.dump(data, f)

            print(f"Wrote synapse parameters to {output_file}")


    def run_user_specified_parameters(self, user_parameter_list):

        if self.pc.id() != 0:
            return

        if len(user_parameter_list) != 5:
            raise ValueError(f"Expected 5 parameters (U, tauR,tauF, tauRatio,cond), got {len(user_parameter_list)}")

        m_params = dict(zip(self.parameter_list, user_parameter_list))

        # Let's also record currents from synapses
        self.syn_recording = []

        for syn in self.rsr_synapse_model.synapses:
            data = neuron.h.Vector()
            data.record(syn._ref_i)
            self.syn_recording.append(data)

        t_sim, v_sim, i_sim = self.rsr_synapse_model.run2(pars=m_params)

        self.last_run_time = t_sim
        self.last_run_volt = v_sim

        self.last_run_parameters = m_params


    def plot_error(self, error_list, fig_name_info="", marker=".", linestyle="-"):

        if self.pc.id() != 0:
            return

        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(error_list, marker=marker, linestyle=linestyle)
        plt.ylabel("Error")

        if self.name_tag is not None:
            name_tag = f"-{self.name_tag}"
        else:
            name_tag = ""

        fig_name = os.path.join("figures", os.path.basename(self.data_file).split(".")[0] + fig_name_info + f"-{self.synapse_type}{name_tag}error.png")

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
    parser.add_argument("--synapse_type", default="glut", help="Specify synapse ['glut', 'glut2']")
    parser.add_argument("--synapse_parameter_file", type=str, default=None)
    parser.add_argument("--neuron_set_file", type=str, default="neuron_set.json")
    parser.add_argument("--name_tag", type=str, default=None, help="Name tag for run")
    parser.add_argument("--user_parameters", default=None,
                        type=lambda s: [float(x) for x in s.split(",")],
                        help="Run user parameters: U,tauR,tauF,tauRatio,nmda_ratio")
    parser.add_argument("--n_synapses", type=int, default=None, help="Override number of synapses in config file")
    parser.add_argument("--n_soma_synapses", type=int, default=None, help="Override number of soma synapses in config file")
    parser.add_argument("--synapse_density", type=str, default=None, help="Override synapse density in config file")
    parser.add_argument("--entropy", type=int, default=1023456734529028340264793840, help="Entropy for random generator")
    parser.add_argument("--export", action="store_true", help="Export synapse parameters")
    parser.add_argument("--profile", action="store_true", default=False)
    parser.add_argument("--verbose", action="store_true", default=False)
    parser.add_argument("--dont-update-neuronset-file", dest="update_neuronset_file", action="store_false", default=True)

    args = parser.parse_args()

    if args.synapse_type == "glut2" and args.synapse_parameter_file is None:
        raise Exception("Synapse parameter file is required for glut2 synapse (tmGlut_double).")

    so = SynapseOptimiser(data_file=args.data_file,
                          snudda_data=args.snudda_data,
                          synapse_type=args.synapse_type,
                          synapse_parameter_file=args.synapse_parameter_file,
                          neuron_set_file=args.neuron_set_file,
                          n_synapses=args.n_synapses,
                          n_soma_synapses=args.n_soma_synapses,
                          synapse_density=args.synapse_density,
                          load_parameters=args.user_parameters is None,
                          name_tag=args.name_tag,
                          update_neuronset_file=args.update_neuronset_file,
                          entropy=args.entropy,
                          verbose=args.verbose)

    if args.user_parameters is not None:
        so.prepare_models()
        so.run_user_specified_parameters(args.user_parameters)
        so.plot_last_run(f"figures/user_specified_parameters-{args.name_tag}.png")
        sys.exit(0)


    if args.profile:
        import cProfile
        prof_file = f"profile-synaptic-opt.prof"
        cProfile.runctx("so.optimise(n_iterations=args.iterations)", globals(), locals(), filename=prof_file)

        if so.pc.id() == 0:
            # To analyse profile data:
            import pstats
            from pstats import SortKey
            p = pstats.Stats(prof_file)
            p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats(100)

    else:
        so.optimise(n_iterations=args.iterations)

    if args.export:
        so.export_best_parameters()


    # mpirun -n 5 python optimise_synapses_parallel.py ../data/synapses/example_data/10_MSN12_GBZ_CC_H20.json --iterations 50 --snudda_data /home/hjorth/HBP/BasalGangliaData/data/

    # Stored:  /media/psf/KTH/2025-11-13-Yvonne-data-teanalysing/Yvonne2019/CategorisedSTP/
    # mpirun -n 5 python optimise_synapses_parallel.py ../data/synapses/example_data/mixed_test_data.json  --iterations 3 --snudda_data /home/hjorth/HBP/BasalGangliaData/data/

    #
    # mpirun -n 5 python optimise_synapses_parallel.py ../data/synapses/example_data/mixed_test_data.json  --iterations 3 --snudda_data /home/hjorth/HBP/BasalGangliaData/data/ --synapse_type glut2 --synapse_parameter_file ../data/synapses/example_data/M1-ipsi_dSPN.json
