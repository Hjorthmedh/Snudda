import numpy as np
import json
import os
import scipy
import scipy.optimize

import pyswarms as ps
from run_little_synapse_run import RunLittleSynapseRun
import matplotlib.pyplot as plt
import matplotlib

from snudda.utils.numpy_encoder import NumpyEncoder

############################################################################


class Planert2010part2(object):

    def __init__(self, data_type, cell_id=None, pretty_plot=True, fit_data=True):

        self.data_legend = {"II": "iSPN to iSPN",
                            "ID": "iSPN to dSPN",
                            "DD": "dSPN to dSPN",
                            "DI": "dSPN to iSPN",
                            "FD": "FSN to dSPN",
                            "FI": "FSN to iSPN"}

        # Since we reuse code that was parallel, need to say that this is master
        self.role = "master"
        self.fig_resolution = 300
        self.pretty_plot = pretty_plot
        self.parameter_cache = dict([])
        self.rsr_synapse_model = None

        file_name = os.path.join("DATA", "Planert2010", f"PlanertFitting-{data_type}-cache.json")
        print(f"Loading data {data_type}: {file_name}")

        with open(file_name, "r") as f:
            self.data = json.load(f)

        if cell_id is None:
            cell_id = np.arange(0, len(self.data["simAmp"]))

        # Also load old parameters if they exist?
        self.cache_file_name = os.path.join("DATA", "Planert2010", f"PlanertFitting-{data_type}-tmgaba-fit.json")

        self.load_parameter_cache()

        if fit_data:
            # Save parameters after fitting

            if type(cell_id) in [list, np.ndarray]:
                for cID in cell_id:
                    self.fit_peaks(data_type, cID)
            else:
                self.fit_peaks(data_type, int(cell_id))

            self.save_parameter_cache()
        else:
            self.setup_model_synapse_fitting(data_type)

            if type(cell_id) in [list, np.ndarray]:
                for cID in cell_id:
                    self.plot_data(data_type, cID, run_sim=True, show=False)
            else:
                self.plot_data(data_type, int(cell_id), run_sim=True, show=False)

    ############################################################################

    def __delete__(self):
        # Save cache file before exiting
        self.save_parameter_cache()

    def add_parameter_cache(self, cell_id, name, value):

        if cell_id not in self.parameter_cache:
            self.parameter_cache[int(cell_id)] = dict([])

        self.parameter_cache[int(cell_id)][name] = value

    ############################################################################

    def get_parameter_cache(self, cell_id, name):

        if cell_id in self.parameter_cache and name in self.parameter_cache[cell_id]:
            return self.parameter_cache[cell_id][name]
        else:
            return None

    ############################################################################

    def load_parameter_cache(self):

        if os.path.exists(self.cache_file_name):
            try:
                print(f"Loading cache file {self.cache_file_name}")
                with open(self.cache_file_name, "r") as f:
                    tmp_dict = json.load(f)

                    self.parameter_cache = dict([])
                    for k in tmp_dict:
                        self.parameter_cache[int(k)] = tmp_dict[k]

                    f.close()
            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)

                print(f"Unable to open {self.cache_file_name}")
                self.parameter_cache = dict([])
        else:
            # No cache file to load, create empty dictionary
            self.parameter_cache = dict([])

            ############################################################################

    def save_parameter_cache(self):

        if self.role != "master":
            print("No servants are allowed to write output to json, ignoring call.")
            return

        print("Saving parameters to cache file: " + str(self.cache_file_name))

        try:
            with open(self.cache_file_name, "w") as f:
                json.dump(self.parameter_cache, f, indent=2, cls=NumpyEncoder)
                f.close()
        except:
            import traceback
            tstr = traceback.format_exc()
            print(tstr)

            print(f"Failed to save cache file ... {self.cache_file_name}")

    ############################################################################

    def fit_peaks(self, data_type, cell_id):

        print(f"Optimising for ID = {cell_id}")

        # Get peaks
        peak_height = np.array(self.data["simAmp"][cell_id])
        t_stim = np.array(self.data["tStim"])

        assert len(t_stim) == len(peak_height), "Inconsistent data lengths for fitPeaks"
        sigma = np.ones(len(peak_height))

        # Setup neuron model
        self.setup_model_synapse_fitting(data_type)

        # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR)
        model_bounds = ([1e-3, 1e-4, 1e-4, 0, 1e-5],
                        [1.0, 2, 2, 0.9999999, 1e-1])

        # Pyswarm options
        options = {"c1": 0.5, "c2": 0.3, "w": 0.9}  # default

        n_particles = 200
        n_iterations = 10

        optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,
                                            dimensions=len(model_bounds[0]),
                                            options=options,
                                            bounds=model_bounds)

        cost, fit_params = optimizer.optimize(self.neuron_synapse_swarm_helper_planert,
                                              iters=n_iterations,
                                              t_stim=t_stim,
                                              peak_height=peak_height)

        model_heights, t_sim, v_sim = self._neuron_synapse_swarm_helper(fit_params, t_stim)

        # tau < tauR, so we use tauRatio for optimisation
        fit_params[3] *= fit_params[1]  # tau = tauR * tauRatio

        print("Parameters: U = %.3g, tauR = %.3g, tauF = %.3g, tau = %.3g, cond = %3.g" % tuple(fit_params))

        print(f"peakHeight = {peak_height}")
        print(f"modelHeights = {model_heights}")

        self.add_parameter_cache(cell_id, "synapse",
                                 {"U": fit_params[0],
                                  "tauR": fit_params[1],
                                  "tauF": fit_params[2],
                                  "tau": fit_params[3],
                                  "cond": fit_params[4]})

        self.add_parameter_cache(cell_id, "surrogatePeaks", peak_height)
        self.add_parameter_cache(cell_id, "fittedPeaks", model_heights)

        self.plot_data(data_type=data_type, cell_id=cell_id,
                       t_stim=t_stim,
                       surrogate_peaks=peak_height,
                       model_peaks=model_heights,
                       t_sim=t_sim,
                       v_sim=v_sim,
                       model_params=fit_params)

    ############################################################################

    def plot_data(self, data_type, cell_id,
                  t_stim=None,
                  surrogate_peaks=None,
                  model_peaks=None,
                  t_sim=None, v_sim=None,
                  model_params=None,
                  show=True,
                  t_skip=0.01,
                  pretty_plot=None,
                  run_sim=False):

        matplotlib.rcParams.update({'font.size': 24})

        if surrogate_peaks is None:
            surrogate_peaks = np.array(self.data["simAmp"][cell_id])

        if t_stim is None:
            t_stim = np.array(self.data["tStim"])

        if model_params is None:
            p_dict = self.get_parameter_cache(cell_id, "synapse")
            model_params = [p_dict[x] for x in ["U", "tauR", "tauF", "tau", "cond"]]

        if run_sim:

            assert model_params is not None, \
                "plot_data: model_params must be given if run_sim = True"
            if t_sim is not None or v_sim is not None:
                print("Ignoring t_sim and v_sim when run_sim = True")

            print(f"Run simulation for {data_type} {cell_id}")
            U, tau_r, tau_f, tau, cond = model_params

            model_peaks, t_sim, v_sim = self.synapse_model_neuron(t_stim, U,
                                                                  tau_r, tau_f, cond, tau,
                                                                  params={},
                                                                  return_trace=True)

        if pretty_plot is None:
            pretty_plot = self.pretty_plot

        plt.figure()
        if t_sim is not None:

            t_idx = np.where(t_sim > t_skip)[0][0]

            v_base = np.mean(v_sim[-50:-1] * 1e3)

            for t, sp, mp in zip(t_stim * 1e3, surrogate_peaks * 1e3, model_peaks * 1e3):
                idx = np.argmin(np.abs(t_sim - t))
                plt.plot([t, t], [v_base, v_base + sp], color=(1, 0.31, 0), linewidth=3)
                # plt.plot([t,t],[vBase,vBase+sp],color="red",linewidth=3)
                plt.plot([t, t], [v_base, v_base + mp], color="blue")

            plt.plot(1e3 * t_sim[t_idx:], 1e3 * v_sim[t_idx:], color="black")

        else:
            for t, sp, mp in zip(t_stim * 1e3, surrogate_peaks * 1e3, model_peaks * 1e3):
                # plt.plot([t,t],[0,sp],color="red",linewidth=3)
                plt.plot([t, t], [0, sp], color=(1, 0.31, 0), linewidth=3)
                plt.plot([t, t], [0, mp], color="blue")

            v_base = 0  # No sim data, base it at 0

        if pretty_plot:

            v_bar_length = 0.05
            if data_type[0] == "F":
                v_bar_length = 0.5

            # Draw scalebars
            v_scale_x = 1200
            # vMax = np.max(vPlot[np.where(tPlot > 0.050)[0]])
            y_scale_bar = v_base + float(np.diff(plt.ylim())) / 4
            v_scale_y1 = y_scale_bar + v_bar_length
            v_scale_y2 = y_scale_bar
            t_scale_y = y_scale_bar
            t_scale_x1 = v_scale_x
            t_scale_x2 = v_scale_x + 100

            plt.plot([v_scale_x, v_scale_x], [v_scale_y1, v_scale_y2], color="black")
            plt.plot([t_scale_x1, t_scale_x2], [t_scale_y, t_scale_y], color="black")

            plt.text(v_scale_x - 100, v_scale_y2 + 0.30 * float(np.diff(plt.ylim())),
                     ("%.2f" % (v_scale_y1 - v_scale_y2)) + " mV",
                     rotation=90)
            plt.text(v_scale_x, v_scale_y2 - float(np.diff(plt.ylim())) / 10,
                     ("%.0f" % (t_scale_x2 - t_scale_x1) + " ms"))

            plt.axis("off")

            if False:
                plt.ion()
                plt.show()
                import pdb
                pdb.set_trace()

        if not pretty_plot:
            if model_params is not None:
                title_str = "\nU=%.3g, tauR=%.3g, tauF=%.3g, tau=%.3g,\ncond=%.3g" \
                           % (model_params[0],
                              model_params[1],
                              model_params[2],
                              model_params[3],
                              model_params[4])
                plt.title(title_str)

            plt.xlabel("Time (ms)")
            plt.ylabel("Volt (mV)")

            # Remove part of the frame
            plt.gca().spines["right"].set_visible(False)
            plt.gca().spines["top"].set_visible(False)

        else:
            if data_type in self.data_legend:
                plt.title(self.data_legend[data_type])

            plt.axis("off")

        if not os.path.exists("figures/"):
            os.makedirs("figures/")

        if pretty_plot:
            fig_name = os.path.join("figures", f"PlanertFitting-{data_type}-{cell_id}-noaxis.pdf")
        else:
            fig_name = os.path.join("figures", f"PlanertFitting-{data_type}-{cell_id}.pdf")

        plt.savefig(fig_name, dpi=self.fig_resolution)

        if show:
            plt.ion()
            plt.show()
        else:
            plt.ioff()
            plt.close()

    ############################################################################

    def neuron_synapse_swarm_helper_planert(self, pars, t_stim, peak_height):

        res = np.zeros((pars.shape[0]))

        for idx, p in enumerate(pars):
            peak_h, t_sim, v_sim = self._neuron_synapse_swarm_helper(p, t_stim)

            # Calculating error in peak height
            h_diff = np.abs(peak_h - peak_height)
            h_diff[0] *= 3
            h_diff[-1] *= 3
            h_error = np.sum(h_diff) / len(h_diff)

            res[idx] = h_error

        return res

    ############################################################################

    def _neuron_synapse_swarm_helper(self,
                                     pars,
                                     t_spikes):

        U, tau_r, tau_f, tau_ratio, cond = pars
        tau = tau_r * tau_ratio
        params = {}

        peak_heights, t_sim, v_sim = self.synapse_model_neuron(t_spikes, U,
                                                               tau_r, tau_f, cond, tau,
                                                               params=params,
                                                               return_trace=True)

        return peak_heights, t_sim, v_sim

    ############################################################################

    def setup_model_synapse_fitting(self, data_type, params=None):

        # If the delta model existed, clear it
        # self.rsrDeltaModel = None

        soma_diameter = 20e-6
        soma_gleak = 3

        if not params:
            params = {"somaDiameter": soma_diameter, "somaGleak": soma_gleak}

        t_stim = np.array(self.data["tStim"])
        max_time = np.max(t_stim) + 0.5

        baseline_depol = -80e-3

        self.rsr_synapse_model = RunLittleSynapseRun(stim_times=t_stim,
                                                     holding_voltage=baseline_depol,
                                                     synapse_type="GABA",
                                                     params=params,
                                                     time=max_time)

    ############################################################################

    def synapse_model_neuron(self, t_spike, U, tau_r, tau_f, cond, tau,
                             params=None,
                             return_trace=False):

        # print("Running neuron model")

        assert self.rsr_synapse_model is not None, \
            "!!! Need to call setupModelSynapseFitting first"

        if not params:
            params = {}

        # Should we make a copy of params, to not destroy it? ;)
        params["U"] = U
        params["tauR"] = tau_r
        params["tauF"] = tau_f
        params["cond"] = cond
        params["tau"] = tau

        # print("params=" + str(params))

        (tSim, vSim, iSim) = \
            self.rsr_synapse_model.run2(pars=params)

        if tSim.shape != vSim.shape:
            print("Shape are different, why?!")
            import pdb
            pdb.set_trace()

        peak_idx = self.get_peak_idx2(time=tSim, volt=vSim, stim_time=t_spike)
        peakHeight, decayFits, vBase = self.find_trace_heights(tSim, vSim, peak_idx)

        if return_trace:
            return peakHeight, tSim, vSim
        else:
            return peakHeight

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

    @staticmethod
    def find_peaks_helper(p_time, p_window, time=None, volt=None):

        peak_idx = []
        peak_time = []
        peak_volt = []

        for pt, pw in zip(p_time, p_window):
            t_start = pt
            t_end = pt + pw

            t_idx = np.where(np.logical_and(t_start <= time, time <= t_end))[0]

            # We assume that neuron is more depolarised than -65, ie gaba is
            # also depolarising
            p_idx = t_idx[np.argmax(volt[t_idx])]

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

    # We are using a simplified function that skips decay fits

    @staticmethod
    def find_trace_heights(time, volt, peak_idx):

        decay_func = lambda x, a, b, c: a * np.exp(-x / b) + c

        v_base = np.mean(volt[int(0.3 * peak_idx[0]):int(0.8 * peak_idx[0])])

        peak_height = np.zeros((len(peak_idx, )))
        peak_height[0] = volt[peak_idx[0]] - v_base

        decay_fits = []

        for idx_b in range(1, len(peak_idx)):

            if peak_height[0] > 0:
                if idx_b < len(peak_idx) - 1:
                    p0d = [0.06, -0.05, -0.074]
                else:
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

                idx_start = int(peak_idx_a * 0.2 + peak_idx_b * 0.8)
                idx_end = int(peak_idx_a * 0.08 + peak_idx_b * 0.92)
            else:
                # Last spike, use only last half of decay trace
                idx_start = int(peak_idx_a * 0.6 + peak_idx_b * 0.4)
                idx_end = int(peak_idx_a * 0.08 + peak_idx_b * 0.92)

            try:
                assert idx_start < idx_end
            except:
                import traceback
                t_str = traceback.format_exc()
                print(t_str)

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
                    fit_params, p_cov = scipy.optimize.curve_fit(decay_func, t_ab_fit, v_ab_fit, p0=p0d)
                except:
                    import traceback
                    t_str = traceback.format_exc()
                    print(t_str)

                    print("!!! Failed to converge, trying with smaller decay constant")
                    p0d[1] *= 0.01
                    fit_params, p_cov = scipy.optimize.curve_fit(decay_func, t_ab_fit, v_ab_fit,
                                                                 p0=p0d)

                t_b = time[peak_idx_b] - t_ab[0]
                v_base_b = decay_func(t_b, fit_params[0], fit_params[1], fit_params[2])

                peak_height[idx_b] = volt[peak_idx_b] - v_base_b

                v_fit = decay_func(t_ab - t_ab[0], fit_params[0], fit_params[1], fit_params[2])
                decay_fits.append((t_ab, v_fit))

                ################################################################

                if False:
                    plt.figure()
                    plt.plot(t_ab, v_ab, 'r')
                    plt.title("Error in findTraceHeights")
                    plt.xlabel("time")
                    plt.ylabel("volt")
                    plt.plot(t_ab, v_fit, 'k-')
                    plt.ion()
                    plt.show()

                    import pdb
                    pdb.set_trace()

                ########################################

            except:

                print("Check that the threshold in the peak detection before is OK")
                # self.plot(name)
                import traceback
                t_str = traceback.format_exc()
                print(t_str)

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

    ############################################################################


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Fit GABA model to peaks")
    parser.add_argument("dataType", choices=["DD", "ID", "DI", "II", "FI", "FD"])
    parser.add_argument("--plotOnly", help="Only plot, no new fitting",
                        action="store_true")
    args = parser.parse_args()

    if args.plotOnly:
        fitData = False
    else:
        fitData = True

    pp = Planert2010part2(data_type=args.dataType, cell_id=None,
                          pretty_plot=True, fit_data=fitData)
    # pp = Planert2010part2(dataType=args.dataType,cellID=0)
