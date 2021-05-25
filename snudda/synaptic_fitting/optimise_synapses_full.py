import os
import numpy as np
import scipy
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib
import neuron
import json
import time

from snudda.synaptic_fitting.parameter_bookkeeper import ParameterBookkeeper


# TODO: 2021-05-12 -- Save more than the best parameter set in json file. Have one dictionary item per saved parameterset,
#                     that way we can remove all associated data easily in one go, when updating the json file.

# TODO: 2021-05-12 -- What value should u0 have? A way around it, repeat the stimulation train multiple times, use laster runs
# TODO: 2021-05-11 -- Set self.synapse_parameters
# TODO: 2021-05-11 -- Holding voltage

# TODO 2020-10-19
#
# We need to make sure params contains the nmda_ratio and other parameters
# that are important for mod file (but not optimised for at this stage)
# Please fix.
#
#
# python3 optimise_synapses_full.py DATA/Yvonne2020/M1RH-ipsi_MSN_D1_20Hz_depressing.json --synapseParameters ../data/synapses/v3/M1RH-ipsi_D1-MSN.json --st glut
#
#
# TODO 2020-10-09
#
# We are no longer reading data from Yvonne's HDF5 directly, instead we are
# reading data from a JSON file, that Ilaria created from Igor exports.
#
# Consequences and TODO:
# - We no longer need CellID to identify which cell we are optmising,
#   each JSON file only contains one dataset (and it is an averaged dataset)
# - Holding voltage is no longer extractable from data, we need to set it
# - Create a JSON optmisation parameter file which contains the holding voltage
#   to be used, as well as the modelbounds (currently in getModelBounds)
# - The JSON data is now loaded into self.data, go through all functions
#   remove references to CellID, and extract the data directly from the
#   self.data variable.
#
#

# TODO 2020-07-02
# -- We just wrote parallelOptimiseSingleCell -- need to make sure we call it
#    the function will optimise one cellID, using all workers available
#    need to debug the code, to make sure it works... have fun! :)
#

#
# TODO 2020-06-16
# -- We need to make sure one neuron can be optimised by all workers
#    collectively, right now one worker does one cell alone
#
# -- Determine synapse locations, pass it to all workers
# -- Best random is easy to parallelise, just do the work, then gather it at
#    the master node.
# -- Difficult: How to parallise scipy.optimize.minimize
#    Possible idea: let func to optimize handle vectors, and then each
#    position in vector is sent to one worker.


# !!! Add check that if the voltage is 0 or 5, then the trace is skipped entirely

from run_synapse_run import RunSynapseRun
from snudda.utils.numpy_encoder import NumpyEncoder


class OptimiseSynapsesFull(object):

    ############################################################################

    # optMethod not used anymore, potentially allow user to set sobol or refine

    # datafile is JSON file from Ilaria's Igor extraction

    def __init__(self, data_file, synapse_type="glut", load_parameters=True,
                 role="master", d_view=None, verbose=True, log_file_name=None,
                 opt_method="sobol", pretty_plot=False,
                 model_bounds="model_bounds.json",
                 neuron_set_file="neuronSet.json",
                 synapse_parameter_file=None):

        # Parallel execution role, "master" or "servant"
        self.role = role

        self.parallel_setup_flag = False  # Set to True when servants are done setup
        self.d_view = d_view
        self.verbose = verbose
        self.log_file_name = log_file_name
        self.opt_method = opt_method
        self.num_smoothing = 200  # How many smoothing points do we use?
        self.sim_time = 1.5
        self.neuron_set_file = neuron_set_file

        self.debug_pars_flag = False
        self.debug_pars = []
        self.cell_properties = None
        self.synapse_parameter_data = None

        self.data = None
        self.volt = None
        self.sample_freq = None
        self.time = None
        self.stim_time = None
        self.cell_type = None

        self.synapse_parameters = None
        self.synapse_section_id = None
        self.synapse_section_x = None

        self.pretty_plot = pretty_plot

        print("Init optMethod = " + str(opt_method))

        if self.log_file_name is not None and len(self.log_file_name) > 0:
            print(f"Log file: {self.log_file_name}")
            self.log_file = open(self.log_file_name, 'w')
        else:
            self.log_file = None

        self.fig_resolution = 300

        self.data_file = data_file
        self.load_trace_data(data_file=data_file)

        self.synapse_parameter_file = synapse_parameter_file

        if synapse_parameter_file:
            with open(synapse_parameter_file, 'r') as f:
                print(f"Reading synapse parameters from {synapse_parameter_file}")
                self.synapse_parameters = json.load(f)["data"]
        else:
            self.synapse_parameters = {}

        self.parameter_data_file_name = f"{self.data_file}-parameters-full.json"
        self.load_parameters = load_parameters
        self.synapse_type = synapse_type

        self.rsr_synapse_model = None
        self.rsr_delta_model = None

        self.model_info = None

        # Sobol params
        self.synapse_model = None
        self.smooth_exp_volt8 = None
        self.smooth_exp_time8 = None
        self.smooth_exp_volt9 = None
        self.smooth_exp_time9 = None

        self.decay_start_fit8 = 0.45
        self.decay_end_fit8 = 0.8
        self.decay_start_fit9 = 1.0
        self.decay_end_fit9 = 1.3

        with open(model_bounds, 'r') as f:
            print(f"Loading model bounds from {model_bounds}")
            self.model_bounds = json.load(f)

        if load_parameters:
            self.load_parameter_data()

        if self.role == "master":
            self.setup_parallel(d_view=d_view)

    ############################################################################

    def __delete__(self):

        # Save the parameter cache before closing
        if self.parameter_data_file_name:
            self.save_parameter_data()
        else:
            print("exiting: parameter_data_file_name not set, not saving parameter data")

        if self.log_file is not None:
            self.log_file.close()
            self.log_file = None

    ############################################################################

    def load_trace_data(self, data_file):

        self.write_log(f"Loading {data_file}")
        with open(data_file, "r") as f:
            self.data = json.load(f)

            self.volt = np.array(self.data["data"]["mean_norm_trace"])
            self.sample_freq = self.data["metadata"]["sample_frequency"]

            dt = 1 / self.sample_freq
            self.time = 0 + dt * np.arange(0, len(self.volt))

            self.stim_time = np.array(self.data["metadata"]["stim_time"])

            self.cell_type = self.data["metadata"]["cell_type"]

    ############################################################################

    def save_parameter_data(self):
        if self.role != "master":
            self.write_log("No servants are allowed to write output to json, ignoring call.")
            return

        print(f"Saving data to {self.parameter_data_file_name}")
        self.synapse_parameter_data.save(self.parameter_data_file_name)

    def load_parameter_data(self):

        self.synapse_parameter_data = ParameterBookkeeper(old_book_file=self.parameter_data_file_name)
        self.synapse_parameter_data.check_integrity()

        best_dataset = self.synapse_parameter_data.get_best_dataset()

        if best_dataset is not None:
            self.synapse_section_id = best_dataset["section_id"]
            self.synapse_section_x = best_dataset["section_x"]

    ############################################################################

    # parDict is the parameters that are associated with cellID

    def plot_data(self,
                  params=None,
                  show=True,
                  skip_time=0.050,
                  pretty_plot=None):

        if params is None:
            params = self.synapse_parameters

        if pretty_plot is None:
            pretty_plot = self.pretty_plot

        if pretty_plot:
            matplotlib.rcParams.update({'font.size': 24})
        else:
            matplotlib.rcParams.update({'font.size': 5})

        if self.volt is None:
            self.write_log("Nothing to plot (volt)")
            return

        best_dataset = self.synapse_parameter_data.get_best_dataset()
        best_params = best_dataset["parameters"]

        #synapse_position_override = (self.get_parameter_cache("sectionID"),
        #                             self.get_parameter_cache("sectionX"))
        dt = best_dataset["dt"]
        v_plot = np.array(best_dataset["volt"])
        t_plot = np.arange(0, len(v_plot)*dt, dt)
        min_error = best_dataset["error"]

        t_idx = np.where(skip_time <= self.time)[0]

        plt.figure()

        plt.plot(self.time[t_idx] * 1e3, self.volt[t_idx] * 1e3, 'r-')
        if v_plot is not None:
            t2_idx = np.where(skip_time <= t_plot)[0]
            plt.plot(t_plot[t2_idx] * 1e3, v_plot[t2_idx] * 1e3, 'k-')

        if not pretty_plot:
            title_str = self.cell_type

            u, tau_r, tau_f, tau_ratio, cond = best_params

            title_str += f"\nU={u:.3g}, tauR={tau_r:.3g}, tauF={tau_f:.3g}, tau={tau_r*tau_ratio:.3g},cond={cond:.3g}"

            plt.title(title_str)

        if pretty_plot:
            # Draw scalebars
            v_scale_x = 1200
            # vMax = np.max(vPlot[np.where(tPlot > 0.050)[0]])
            v_base = v_plot[-1]
            y_scale_bar = v_base * 1e3 + float(np.diff(plt.ylim())) / 4
            v_scale_y1 = y_scale_bar + 1
            v_scale_y2 = y_scale_bar
            t_scale_y = y_scale_bar
            t_scale_x1 = v_scale_x
            t_scale_x2 = v_scale_x + 100

            plt.plot([v_scale_x, v_scale_x], [v_scale_y1, v_scale_y2], color="black")
            plt.plot([t_scale_x1, t_scale_x2], [t_scale_y, t_scale_y], color="black")

            plt.text(v_scale_x - 100, v_scale_y2 + 0.20 * float(np.diff(plt.ylim())),
                     ("%.0f" % (v_scale_y1 - v_scale_y2)) + " mV",
                     rotation=90)
            plt.text(v_scale_x, v_scale_y2 - float(np.diff(plt.ylim())) / 10,
                     ("%.0f" % (t_scale_x2 - t_scale_x1) + " ms"))

            # Mark optogenetical stimulation
            y_height = float(np.diff(plt.ylim())) / 13

            t_stim = self.stim_time
            y_stim_marker1 = v_plot[-1] * 1e3 - 1.5 * y_height
            y_stim_marker2 = v_plot[-1] * 1e3 - 2.5 * y_height
            for ts in t_stim:
                plt.plot([ts, ts], [y_stim_marker1, y_stim_marker2], color="cyan")

            plt.axis("off")

        plt.xlabel("Time (ms)")
        plt.ylabel("Volt (mV)")

        if not os.path.exists("figures/"):
            os.makedirs("figures/")

        base_name = os.path.splitext(os.path.basename(self.data_file))[0]
        fig_name = f"figures/{base_name}.pdf"
        plt.savefig(fig_name, dpi=self.fig_resolution)

        if show:
            plt.ion()
            plt.show()
        else:
            plt.ioff()
            plt.close()

    ############################################################################

    def get_cell_properties(self):

        if self.cell_properties is None:
            with open(self.neuron_set_file, 'r') as f:
                self.cell_properties = json.load(f)

        cell_type = self.data["metadata"]["cell_type"]

        return self.cell_properties[cell_type].copy()

    ############################################################################

    def extract_input_res_tau(self, t, v, cur_amp, cur_start, cur_end, base_start, base_end):

        # We assume SI units
        t_idx_base = np.where(np.logical_and(base_start < t, t < base_end))[0]
        v_base = np.mean([v[x] for x in t_idx_base])

        t_idx_peak = np.where(np.logical_and(cur_start < t, t < cur_end))[0]
        v_peak = np.min(v[t_idx_peak])
        # vPeak = np.max([v[x] for x in tIdxPeak])

        assert np.abs(v_peak - v_base) > np.abs(np.max(v[t_idx_peak]) - v_base), \
            "The code assumes a hyperpolarising pulse, not a peak maximum"

        rm = (v_peak - v_base) / cur_amp

        idx_post_pulse = np.where(cur_end < t)[0]
        idx_max_post_pulse = idx_post_pulse[np.argmax(v[idx_post_pulse])]
        t_max_post_pulse = t[idx_max_post_pulse]

        t_idx_decay = np.where(np.logical_and(cur_end < t, t < t_max_post_pulse))[0]

        decay_func = lambda x, a, b, c: a * np.exp(-x / b) + c

        t_ab_fit = t[t_idx_decay] - t[t_idx_decay[0]]
        v_ab_fit = v[t_idx_decay]

        p0d = [-0.005, 0.01, -0.06]

        if np.isnan(v_ab_fit).any() or np.isinf(v_ab_fit).any():
            self.write_log("We have inifinite or nan values in the voltage")
            import pdb
            pdb.set_trace()

        try:
            fit_params, pcov = scipy.optimize.curve_fit(decay_func, t_ab_fit, v_ab_fit, p0=p0d)
            tau = fit_params[1]

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            plt.figure()
            plt.plot(t_ab_fit, v_ab_fit)
            plt.ion()
            plt.show()

            import pdb
            pdb.set_trace()

        if False:
            plt.figure()
            plt.plot(t_ab_fit, v_ab_fit, '-r')
            plt.plot(t_ab_fit, decay_func(t_ab_fit,
                                          fit_params[0],
                                          fit_params[1],
                                          fit_params[2]))
            plt.xlabel("t")
            plt.ylabel("v")
            plt.title("Tau extraction")
            plt.ion()
            plt.show()

            # self.writeLog("RM = " + str(RM) + " tau = " + str(tau))

        # Return membrane resistance and tau
        return rm, tau

    ############################################################################

    def get_peak_idx(self):

        p_time = np.array(self.data["metadata"]["stim_time"]) * 1e-3
        freq = self.data["metadata"]["freq"]

        assert np.abs(1.0 - freq / (p_time[1] - p_time[0])) < 0.01, "frequency mismatch"

        p_window = 1.0 / (2 * freq) * np.ones(p_time.shape)
        p_window[-1] *= 5

        peak_info = self.find_peaks_helper(p_time=p_time,
                                           p_window=p_window,
                                           time=self.time,
                                           volt=self.volt)

        return peak_info["peakIdx"]

    ############################################################################

    def get_peak_idx2(self, stim_time, time, volt):

        freq = 1.0 / (stim_time[1] - stim_time[0])

        p_window = 1.0 / (2 * freq) * np.ones(stim_time.shape)
        p_window[-1] *= 5

        peak_info = self.find_peaks_helper(p_time=stim_time,
                                           p_window=p_window,
                                           time=time,
                                           volt=volt)

        return peak_info["peakIdx"]

    ############################################################################

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
            assert len(t_idx) > 0, f"No time points within {t_start} and {t_end}"

            if self.synapse_type == "glut":
                p_idx = t_idx[np.argmax(volt[t_idx])]
            elif self.synapse_type == "gaba":
                # We assume that neuron is more depolarised than -65, ie gaba is
                # also depolarising
                p_idx = t_idx[np.argmax(volt[t_idx])]
            else:
                self.write_log("Unknown synapse type : " + str(self.synapse_type))
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

    # (peakIdx,peakTime,peakVolt)

    ############################################################################

    def find_trace_heights(self, time, volt, peak_idx):

        decay_func = lambda x, a, b, c: a * np.exp(-x / b) + c

        v_base = np.mean(volt[int(0.3 * peak_idx[0]):int(0.8 * peak_idx[0])])

        peak_height = np.zeros((len(peak_idx, )))
        peak_height[0] = volt[peak_idx[0]] - v_base

        decay_fits = []

        for idx_b in range(1, len(peak_idx)):

            if peak_height[0] > 0:
                if idx_b < len(peak_idx) - 1:
                    p0d = [0.06, 0.05, -0.074]
                else:
                    p0d = [1e-5, 100, -0.074]

                    if self.synapse_type == "gaba":
                        p0d = [1e-8, 10000, -0.0798]
            else:
                # In some cases for GABA we had really fast decay back
                if idx_b < len(peak_idx) - 1:
                    p0d = [-0.06, 0.05, -0.0798]
                else:
                    p0d = [-1e-5, 1e5, -0.0798]

            idx_a = idx_b - 1

            peak_idx_a = peak_idx[idx_b - 1]  # Prior peak
            peak_idx_b = peak_idx[idx_b]  # Next peak

            if idx_b < len(peak_idx) - 1:
                # Not the last spike
                idx_start = int(peak_idx_a * 0.9 + peak_idx_b * 0.1)
                idx_end = int(peak_idx_a * 0.1 + peak_idx_b * 0.9)
            else:
                # Last spike, use only last half of decay trace
                idx_start = int(peak_idx_a * 0.5 + peak_idx_b * 0.5)
                idx_end = int(peak_idx_a * 0.05 + peak_idx_b * 0.85)  # might need 0.85 as last

            try:
                assert idx_start < idx_end
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)

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
                self.write_log(tstr)

                if True:
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

        # import pdb
        # pdb.set_trace()

        return peak_height.copy(), decay_fits, v_base

    ############################################################################

    def setup_model(self, params=None,
                    synapse_density_override=None,
                    n_synapses_override=None,
                    synapse_position_override=None):

        if params is None:
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
            synapse_density = c_prop["synapseDensity"]

        if n_synapses_override is not None:
            n_synapses = n_synapses_override
        else:
            n_synapses = c_prop["nSynapses"]

        # !!! We need to get the baseline depolarisation in another way

        self.rsr_synapse_model = \
            RunSynapseRun(neuron_morphology=c_prop["neuronMorphology"],
                          neuron_mechanisms=c_prop["neuronMechanisms"],
                          neuron_parameters=c_prop["neuronParameters"],
                          neuron_modulation=c_prop["neuronModulation"],
                          stim_times=t_stim,
                          num_synapses=n_synapses,
                          synapse_density=synapse_density,
                          holding_voltage=c_prop["baselineVoltage"],
                          synapse_type=self.synapse_type,
                          params=params,
                          time=self.sim_time,
                          log_file=self.log_file,
                          synapse_section_id=synapse_section_id,
                          synapse_section_x=synapse_section_x)

        return self.rsr_synapse_model

    ############################################################################

    def neuron_synapse_swarm_helper(self, pars, t_spikes, peak_height, smooth_exp_trace):

        assert False, "Is this still used?"

        if self.debug_pars_flag:
            self.debug_pars.append(pars)

        res = np.zeros((pars.shape[0]))

        for idx, p in enumerate(pars):
            peak_h, t_sim, v_sim = self._neuron_synapse_swarm_helper(p, t_spikes)

            # Calculating error in peak height
            h_diff = np.abs(peak_h - peak_height)
            h_diff[0] *= 3
            h_diff[-1] *= 3
            h_error = np.sum(h_diff) / len(h_diff)

            # Calculate error in decay fit
            sim_trace, sim_time = self.smoothing_trace(v_sim, self.num_smoothing,
                                                       time=t_sim,
                                                       start_time=self.decayStartFit,
                                                       end_time=self.decayEndFit)

            # We only want to use the bit of the trace after max
            idx_max = np.argmax(smooth_exp_trace)

            # We divide by number of points in vector, to get the average deviation
            # then we multiply by 10000 to get an error comparable to the others
            decay_error = np.sum((smooth_exp_trace[idx_max:] - sim_trace[idx_max:]) ** 2) \
                          / (self.num_smoothing - idx_max + 1) * 2000

            if False:
                plt.figure()
                plt.plot(smooth_exp_trace[idx_max:])
                plt.plot(sim_trace[idx_max:])
                plt.ion()
                plt.show()
                import pdb
                pdb.set_trace()

            res[idx] = h_error + decay_error

        return res

    ############################################################################

    def smoothing_trace(self, original_trace, num_parts, time=None, start_time=None, end_time=None):

        if time is not None:
            t_flag = np.ones((len(original_trace),), dtype=bool)

            if end_time is not None:
                t_flag[np.where(time > end_time)[0]] = False

            if start_time is not None:
                t_flag[np.where(time < start_time)[0]] = False

            # tIdx = np.where(tFlag)[0]
            trace = original_trace[t_flag]
            t = time[t_flag]
        else:
            trace = original_trace
            t = time

        N = int(np.round(len(trace) / num_parts))

        smooth_trace = np.convolve(trace, np.ones((N,)) / N, mode='valid')

        idx = np.linspace(0, len(smooth_trace) - 1, num=num_parts, dtype=int)

        return smooth_trace[idx], t[idx]

    ############################################################################

    def _neuron_synapse_swarm_helper(self,
                                     pars,
                                     t_spikes):

        assert False, "Remove this function?"

        u, tau_r, tau_f, tau_ratio, cond = pars
        tau = tau_r * tau_ratio

        # TODO: Check where parms come from, is it a dictionary?
        peak_heights, t_sim, v_sim = self.run_model(t_spikes, u,
                                                    tau_r, tau_f, cond, tau,
                                                    params=params,
                                                    return_trace=True)

        return peak_heights, t_sim, v_sim

    ############################################################################

    def neuron_synapse_helper(self,
                              t_spike, u, tau_r, tau_f,
                              tau_ratio=None,
                              cond=1e-7, tau=None):

        assert tau_ratio is None or tau is None, \
            "Only one of tau and tauRatio should be defined"

        if tau_ratio is not None:
            tau = tau_r * tau_ratio
        elif tau is None:
            assert False, "tau or tauRatio must be specified"

        peak_heights = self.run_model(t_spike, u, tau_r, tau_f, cond, tau)

        return peak_heights

    ############################################################################

    # !!! OBS tauRatio is inparameter

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
        h_diff = np.abs(peak_h - exp_peak_height)
        h_diff[0] *= 3
        h_diff[-2] *= 2
        h_diff[-1] *= 3

        # This is to prevent the model spiking
        spike_penalty = np.sum(peak_h > 0.03) * 1

        h_error = np.sum(h_diff) / len(h_diff)

        decay_error8 = np.sum((smooth_exp_trace8[idx_max8:] - sim_trace8[idx_max8:]) ** 2) \
                        / (self.num_smoothing - idx_max8 + 1) * 10000

        decay_error9 = np.sum((smooth_exp_trace9[idx_max9:] - sim_trace9[idx_max9:]) ** 2) \
                        / (self.num_smoothing - idx_max9 + 1) * 10000

        fit_error = h_error + decay_error8 + decay_error9 + spike_penalty

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

    ############################################################################

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

        (t_sim, v_sim, i_sim) = self.rsr_synapse_model.run2(pars=params)

        if t_sim.shape != v_sim.shape:
            self.write_log("Shape are different, why?!")
            import pdb
            pdb.set_trace()

        peak_idx = self.get_peak_idx2(time=t_sim, volt=v_sim, stim_time=t_spike)
        peak_height, decay_fits, v_base = self.find_trace_heights(t_sim, v_sim, peak_idx)

        if return_trace:
            return peak_height, t_sim, v_sim
        else:
            return peak_height

    ############################################################################

    # This should read from a JSON file instead

    def get_model_bounds(self):

        mb = self.data["modeldata"]

        param_list = ["U", "tauR", "tauF", "tauRatio", "cond"]
        lower_bound = [mb[x][0] for x in param_list]
        upper_bound = [mb[x][1] for x in param_list]

        return lower_bound, upper_bound

    ############################################################################

    def sobol_scan(self, synapse_model,
                   t_stim, h_peak,
                   model_bounds,
                   smooth_exp_trace8, smooth_exp_trace9,
                   n_trials=5, load_params_flag=False,
                   parameter_sets=None,
                   return_min_error=False):

        assert self.synapse_type == "glut", \
            "GABA synapse not supported yet in new version"

        print(f"sobol_scan n_trials = {n_trials}")

        if parameter_sets is None:
            parameter_sets = self.setup_parameter_set(model_bounds, n_trials)

        # zip(*xxx) unzips xxx -- cool.
        u_sobol, tau_r_sobol, tau_f_sobol, tau_ratio_sobol, cond_sobol = zip(*parameter_sets)

        # tauSobol = np.multiply(tauRatioSobol,tauRSobol)

        min_pars = None
        min_error = np.inf

        if load_params_flag:
            # If we should load params then do so first
            min_pars = self.synapse_parameter_data.get_best_parameterset()

            # What was error of the cached parameterset
            if min_pars is not None:
                min_error = self.neuron_synapse_helper_glut(t_stim,
                                                            u=min_pars[0],
                                                            tau_r=min_pars[1],
                                                            tau_f=min_pars[2],
                                                            tau_ratio=min_pars[3] / min_pars[1],
                                                            cond=min_pars[4],
                                                            smooth_exp_trace8=smooth_exp_trace8,
                                                            smooth_exp_trace9=smooth_exp_trace9,
                                                            exp_peak_height=h_peak,
                                                            return_type="error")

        idx = 0

        for u, tau_r, tau_f, tau_ratio, cond \
                in zip(u_sobol, tau_r_sobol, tau_f_sobol, tau_ratio_sobol, cond_sobol):

            idx += 1
            if idx % 50 == 0:
                self.write_log("%d / %d : minError = %g" % (idx, len(u_sobol), min_error))

            error, peaks, t, v = self.neuron_synapse_helper_glut(t_stim, u, tau_r, tau_f, tau_ratio, cond,
                                                                 smooth_exp_trace8=smooth_exp_trace8,
                                                                 smooth_exp_trace9=smooth_exp_trace9,
                                                                 exp_peak_height=h_peak,
                                                                 return_type="full")

            min_error = np.minimum(error, min_error)

            param_this_run = np.array([u, tau_r, tau_f, tau_ratio, cond])

            self.synapse_parameter_data.add_parameters(parameter_set=param_this_run,
                                                       section_id=self.rsr_synapse_model.synapse_section_id,
                                                       section_x=self.rsr_synapse_model.synapse_section_x,
                                                       error=error,
                                                       dt=t[1] - t[0],
                                                       volt=v)

        return self.synapse_parameter_data.book

    ############################################################################

    def best_random(self, synapse_model,
                    t_stim, h_peak,
                    model_bounds,
                    n_trials=5000, load_params_flag=False):

        assert n_trials >= 1, "nTrials should be a positive integer"

        min_error = np.inf
        min_par = None

        for idx in range(0, n_trials):
            if idx % 100 == 0:
                self.write_log(f"Pre-trial : {idx}/{n_trials}")

            if idx == 0 and load_params_flag:
                # If we should load params then do so first
                pars = self.get_parameter_cache("synapse")
            else:
                pars = None

            if pars is not None:
                u = pars["U"]
                tau_r = pars["tauR"]
                tau_f = pars["tauF"]
                tau = pars["tau"]
                cond = pars["cond"]
            else:
                u = np.random.uniform(model_bounds[0][0], model_bounds[1][0])
                tau_r = np.random.uniform(model_bounds[0][1], model_bounds[1][1])
                tau_f = np.random.uniform(model_bounds[0][2], model_bounds[1][2])
                tau = tau_r * np.random.uniform(model_bounds[0][3], model_bounds[1][3])
                cond = np.random.uniform(model_bounds[0][4], model_bounds[1][4])

            try:
                peak_heights = self.run_model(t_stim, u, tau_r, tau_f, tau, cond)

                error = np.abs(peak_heights - h_peak)
                error[0] *= 3
                error[1] *= 2
                error[-1] *= 3
                error = np.sum(error)

                if error < min_error:
                    min_error = error
                    min_par = np.array([u, tau_r, tau_f, tau / tau_r, cond])

            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)

                # For big runs we do no want to give up. Let's try again...
                continue

        return min_par

    ############################################################################

    def parallel_optimise_single_cell(self, n_trials=10000, post_opt=False):

        # !!! Future improvement. Allow continuation of old optimisation by
        # reading synapse location and old parameter set, so that is not thrown away

        if self.role == "master":

            # 1. Setup workers
            params = self.synapse_parameters

            if self.synapse_section_id is not None:
                syn_override = self.synapse_section_id, self.synapse_section_x

            # 2. Setup one cell to optimise, randomise synapse positions
            synapse_model = self.setup_model(params=params,
                                             synapse_position_override=syn_override)

            # (volt,time) = self.getData(dataType,cellID)
            peak_idx = self.get_peak_idx2(stim_time=self.stim_time,
                                          time=self.time,
                                          volt=self.volt)
            t_spikes = self.time[peak_idx]

            sigma = np.ones(len(peak_idx))
            sigma[-1] = 1. / 3

            peak_height, decay_fits, v_base = self.find_trace_heights(self.time, self.volt, peak_idx)

            # 2b. Create list of all parameter points to investigate
            model_bounds = self.get_model_bounds()
            parameter_points = self.setup_parameter_set(model_bounds, n_trials)

            # 3. Send synapse positions to all workers, and split parameter points
            #    between workers

            if self.d_view is not None:
                self.setup_parallel(self.d_view)

                self.d_view.scatter("parameterPoints", parameter_points, block=True)

                self.d_view.push({"params": params,
                                 "synapseSectionID": synapse_model.synapse_section_id,
                                 "synapseSectionX": synapse_model.synapse_section_x,
                                 "modelBounds": model_bounds,
                                 "stimTime": self.stim_time,
                                 "peakHeight": peak_height},
                                 block=True)

                cmd_str_setup = \
                    "ly.sobol_worker_setup(params=params," \
                    + "synapsePositionOverride=(synapseSectionID,synapseSectionX))"

                self.d_view.execute(cmd_str_setup, block=True)

                cmd_str = "res = ly.sobolScan(synapseModel=ly.synapseModel, \
                                     tStim = stimTime, \
                                     hPeak = peakHeight, \
                                     parameterSets=parameterPoints, \
                                     modelBounds=modelBounds, \
                                     smoothExpTrace8=ly.smoothExpVolt8, \
                                     smoothExpTrace9=ly.smoothExpVolt9, \
                                     returnMinError=True)"

                self.write_log("Executing workers, bang bang")
                self.d_view.execute(cmd_str, block=True)

                # 5. Gather worker data
                self.write_log("Gathering results from workers")
                res = self.d_view["res"]

                for r in res:
                    self.synapse_parameter_data.merge(r)

                self.save_parameter_data()

                import pdb
                pdb.set_trace()

            else:

                # No dView, run in serial mode...
                self.sobol_worker_setup(params=params,
                                        synapse_position_override=(synapse_model.synapse_section_id,
                                                                   synapse_model.synapse_section_x))

                self.sobol_scan(synapse_model=synapse_model,
                                t_stim=self.stim_time,
                                h_peak=peak_height,
                                model_bounds=model_bounds,
                                smooth_exp_trace8=ly.smooth_exp_volt8,
                                smooth_exp_trace9=ly.smooth_exp_volt9,
                                return_min_error=True)

            self.write_log(f"Sobol search done. Best parameter {self.synapse_parameter_data.get_best_parameterset()}")

            if post_opt:
                # This updates parameters and saves new parameter cache
                self.get_refined_parameters()
                self.save_parameter_data()

    ############################################################################

    # This sets up the model also, so can be run in a self-contained loop
    # We might later want to let the workers do this, but then they cant
    # write to cache --- THAT WILL LEAD TO DATA CORRUPTION!!

    def get_refined_parameters(self):

        assert self.role == "master", \
            "You do not want to run this on workers in parallel, " \
            + " since it writes directly to parameter cache. " \
            + " That could lead to corrupted data."

        # Load parameters from disk
        self.load_parameter_data()
        model_bounds = self.get_model_bounds()

        data_set = self.synapse_parameter_data.get_best_dataset()

        start_par = data_set["parameters"]
        section_x = data_set["section_x"]
        section_id = data_set["section_id"]
        start_par_error_val = data_set["error"]

        synapse_position_override = (section_id, section_x)

        # Make sure we have correct taus etc for synapse
        params = self.synapse_parameters

        peak_idx = self.get_peak_idx2(stim_time=self.stim_time,
                                      time=self.time,
                                      volt=self.volt)
        t_spikes = self.time[peak_idx]

        peak_height, decay_fits, v_base = self.find_trace_heights(self.time, self.volt, peak_idx)

        self.sobol_worker_setup(params, synapse_position_override=synapse_position_override)

        func = lambda x: \
            self.neuron_synapse_helper_glut(t_spike=self.stim_time,
                                            u=x[0],
                                            tau_r=x[1],
                                            tau_f=x[2],
                                            tau_ratio=x[3],
                                            cond=x[4],
                                            smooth_exp_trace8=self.smooth_exp_volt8,
                                            smooth_exp_trace9=self.smooth_exp_volt9,
                                            exp_peak_height=peak_height,
                                            return_type="error")

        m_bounds = [x for x in zip(model_bounds[0], model_bounds[1])]

        res = scipy.optimize.minimize(func,
                                      x0=start_par,
                                      bounds=m_bounds)

        fit_params = res.x
        min_error = res.fun

        if min_error >= start_par_error_val:
            print("Refinement failed. Sobol parameters are better match than new fitting")
            # Dont overwrite the old parameters

        else:

            fit_error, peak_h, t_sim, v_sim = self.neuron_synapse_helper_glut(t_spike=self.stim_time,
                                                                              u=fit_params[0],
                                                                              tau_r=fit_params[1],
                                                                              tau_f=fit_params[2],
                                                                              tau_ratio=fit_params[3],
                                                                              cond=fit_params[4],
                                                                              smooth_exp_trace8=self.smooth_exp_volt8,
                                                                              smooth_exp_trace9=self.smooth_exp_volt9,
                                                                              exp_peak_height=peak_height,
                                                                              return_type="full")

            self.synapse_parameter_data.add_parameters(parameter_set=fit_params,
                                                       section_x=section_x,
                                                       section_id=section_id,
                                                       error=min_error,
                                                       dt=t_sim[1]-t_sim[0],
                                                       volt=v_sim)

            print(f"Old error: {start_par_error_val}, New error: {min_error}")

    ############################################################################

    def sobol_worker_setup(self, params, synapse_position_override=None):

        # TODO: These variables should be defined as None in init
        self.synapse_model = self.setup_model(params=params,
                                              synapse_position_override=synapse_position_override)
        self.smooth_exp_volt8, self.smooth_exp_time8 \
            = self.smoothing_trace(self.volt, self.num_smoothing,
                                   time=self.time,
                                   start_time=self.decay_start_fit8,
                                   end_time=self.decay_end_fit8)

        self.smooth_exp_volt9, self.smooth_exp_time9 \
            = self.smoothing_trace(self.volt, self.num_smoothing,
                                   time=self.time,
                                   start_time=self.decay_start_fit9,
                                   end_time=self.decay_end_fit9)

    ############################################################################

    def setup_parameter_set(self, model_bounds, n_sets):

        import chaospy
        distribution = chaospy.J(chaospy.Uniform(model_bounds[0][0],
                                                 model_bounds[1][0]),
                                 chaospy.Uniform(model_bounds[0][1],
                                                 model_bounds[1][1]),
                                 chaospy.Uniform(model_bounds[0][2],
                                                 model_bounds[1][2]),
                                 chaospy.Uniform(model_bounds[0][3],
                                                 model_bounds[1][3]),
                                 chaospy.Uniform(model_bounds[0][4],
                                                 model_bounds[1][4]))

        u_sobol, tau_r_sobol, tau_f_sobol, tau_ratio_sobol, cond_sobol \
            = distribution.sample(n_sets, rule="sobol")

        parameter_sets = [x for x in zip(u_sobol,
                                         tau_r_sobol, tau_f_sobol, tau_ratio_sobol,
                                         cond_sobol)]

        return parameter_sets

    def setup_parallel(self, d_view=None):

        assert self.role == "master", "Only master should call setupParallel"

        if d_view is None:
            self.write_log("No dView, no parallel")
            return
        else:
            self.d_view = d_view

        if self.parallel_setup_flag:
            # Already setup servants
            return

        with self.d_view.sync_imports():
            from run_synapse_run import RunSynapseRun
            from optimise_synapses_full import NumpyEncoder
            from optimise_synapses_full import OptimiseSynapsesFull

        self.write_log("Setting up workers: " \
                       + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        # Create unique log file names for the workers
        if self.log_file_name is not None:
            engine_log_file = [self.log_file_name + "-" \
                             + str(x) for x in range(0, len(self.d_view))]
        else:
            engine_log_file = [[] for x in range(0, len(self.d_view))]

        n_workers = len(self.d_view)
        self.d_view.scatter("engineLogFile", engine_log_file)

        self.d_view.push({"datafile": self.data_file,
                         "synapseType": self.synapse_type,
                         "synapseparameters": self.synapse_parameter_file,
                         "loadCache": self.load_cache,
                         "role": "servant"})

        cmd_str = "ly = OptimiseSynapsesFull(datafile=datafile, synapseParameterFile=synapseparameters, synapseType=synapseType,loadCache=loadCache,role=role,logFileName=engineLogFile[0])"
        self.d_view.execute(cmd_str, block=True)
        self.parallel_setup_flag = True

    ############################################################################

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

                ############################################################################

    def plot_debug_pars(self):

        try:
            n_iter = len(self.debug_pars)
            n_points = self.debug_pars[0].shape[0]
            n_pars = self.debug_pars[0].shape[1]

            assert n_pars == 6, "Should be six parameters"

            u_all = np.zeros((n_points, n_iter))
            tau_r = np.zeros((n_points, n_iter))
            tau_f = np.zeros((n_points, n_iter))
            tau_ratio = np.zeros((n_points, n_iter))
            cond = np.zeros((n_points, n_iter))

            for ctr, par in enumerate(self.debug_pars):
                u_all[:, ctr] = par[:, 0]
                tau_r[:, ctr] = par[:, 1]
                tau_f[:, ctr] = par[:, 2]
                tau_ratio[:, ctr] = par[:, 3]
                cond[:, ctr] = par[:, 4]

            plt.figure()
            plt.plot(u_all, cond, '-')
            plt.xlabel("U")
            plt.ylabel("cond")

            plt.figure()
            plt.plot(tau_r, tau_f, '-')
            plt.xlabel("tauR")
            plt.ylabel("tauF")

            plt.figure()
            plt.plot(u_all)
            plt.xlabel("Uall")

            plt.ion()
            plt.show()

            import pdb
            pdb.set_trace()

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            import pdb
            pdb.set_trace()

    ############################################################################


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Extract synaptic parameters from electrophysiological recordings")
    parser.add_argument("datafile", help="JSON DATA file")
    parser.add_argument("--synapseParameters", help="Static synapse parameters (JSON)",
                        default=None)
    parser.add_argument("--st", help="Synapse type (glut or gaba)",
                        choices=["glut", "gaba"], default="glut")
    parser.add_argument("--optMethod",
                        help="Optimisation method",
                        choices=["sobol", "stupid", "swarm"],
                        default="sobol")
    parser.add_argument("--plot", action="store_true",
                        help="plotting previous optimised model")
    parser.add_argument("--prettyplot", action="store_true",
                        help="plotting traces for article")

    args = parser.parse_args()

    optMethod = args.optMethod

    print(f"Reading file : {args.datafile}")
    print(f"Synapse type : {args.st}")
    print(f"Synapse params : {args.synapseParameters}")
    print(f"Optimisation method : {optMethod}")

    print("IPYTHON_PROFILE = " + str(os.getenv('IPYTHON_PROFILE')))

    if (os.getenv('IPYTHON_PROFILE') is not None or os.getenv('SLURMID') is not None):
        from ipyparallel import Client

        rc = Client(profile=os.getenv('IPYTHON_PROFILE'),
                    debug=False)

        # http://davidmasad.com/blog/simulation-with-ipyparallel/
        # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
        d_view = rc.direct_view(targets='all')  # rc[:] # Direct view into clients
        lb_view = rc.load_balanced_view(targets='all')
    else:
        d_view = None

    log_file_name = os.path.join("logs", f"{os.path.basename(args.datafile)}-log.txt")
    if not os.path.exists("logs/"):
        os.makedirs("logs/")

    ly = OptimiseSynapsesFull(data_file=args.datafile,
                              synapse_parameter_file=args.synapseParameters,
                              synapse_type=args.st, d_view=d_view,
                              role="master",
                              log_file_name=log_file_name, opt_method=optMethod)

    if args.plot or args.prettyplot:

        if args.prettyplot:
            pretty_plot_flag = True
        else:
            pretty_plot_flag = False

        ly.plot_data(show=True, pretty_plot=pretty_plot_flag)

        exit(0)

    ly.parallel_optimise_single_cell(n_trials=1)
