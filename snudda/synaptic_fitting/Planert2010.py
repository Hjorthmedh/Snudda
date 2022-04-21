# Plan:
#
# Randomise dtau, ftau, U -- check that it is within the ranges specified
# for pair pulse ratio, and recovery test ratio
#
# Fit the neuron synapse model to the modelled traces
#
#

import json
import os

import numpy as np
import scipy.stats as stats

import matplotlib.pyplot as plt

from snudda.utils.numpy_encoder import NumpyEncoder

################################################################################


class Planert2010(object):

    def __init__(self, par_type, pre_num=1000000, max_num=100):

        self.con_info = None
        self.t_stim = None  # Set in def_exp_data

        self.def_exp_data()
        pars = self.pick_random_param(par_type=par_type, num=pre_num)
        amps = None

        if True:
            # Just checking that original starting distributions are ok
            self.plot_model_params(par_type, pars)

        amps = self.check_data_match(par_type=par_type, pars=pars, amps=amps)

        for ctr in range(0, 3):
            if False:
                # This tries to prune the distribution of the variables back to
                # the range
                pars, amps = self.prune_parameters(par_type=par_type, pars=pars, amps=amps)
                self.check_data_match(par_type=par_type, pars=pars, amps=amps)

            pars, amps = self.prune_data(par_type=par_type, pars=pars, amps=amps)
            self.check_data_match(par_type=par_type, pars=pars, amps=amps)
            self.plot_model_params(par_type, pars)

        pars, amps = self.max_num_params(pars, amps, max_num=max_num)

        self.save_data(par_type, pars, amps)

        self.check_data_match(par_type=par_type, pars=pars, amps=amps)
        self.plot_model_params(par_type, pars)

    ############################################################################

    def prune_data(self, par_type, pars, amps=None):

        if amps.shape[0] == 0:
            print("No data left?")
            import pdb
            pdb.set_trace()

        ppr_mean, ppr_std = self.con_info[par_type]["ppr"]
        rtr_mean, rtr_std = self.con_info[par_type]["rtr"]

        if amps is None:
            amps = self.synapse_model_wrapper(self.t_stim, pars)

        ppr = np.divide(amps[:, 1], amps[:, 0])
        rtr = np.divide(amps[:, -1], amps[:, 0])

        n_bins = 20
        ppr_max = np.max(ppr)
        rtr_max = np.max(rtr)

        # Which bin does each ppr belong to
        ppr_idx = (n_bins * ppr / (ppr_max + 1e-6)).astype(int)
        rtr_idx = (n_bins * rtr / (rtr_max + 1e-6)).astype(int)

        ppr_count = np.zeros((n_bins,))
        rtr_count = np.zeros((n_bins,))

        for p, r in zip(ppr_idx, rtr_idx):
            ppr_count[p] += 1
            rtr_count[r] += 1

        ppr_bin_width = ppr_max / (n_bins - 1)
        rtr_bin_width = rtr_max / (n_bins - 1)

        ppr_centre = ppr_bin_width / 2 + np.arange(0, n_bins) * ppr_bin_width
        rtr_centre = rtr_bin_width / 2 + np.arange(0, n_bins) * rtr_bin_width

        ppr_density = stats.norm.pdf(ppr_centre, ppr_mean, ppr_std)
        rtr_density = stats.norm.pdf(rtr_centre, rtr_mean, rtr_std)

        # We pick random param sets, and see if that PPR and RTR are over
        # represented, if so we remove it.

        n_left = len(pars["u"])
        keep_flag = np.ones((n_left,), dtype=bool)

        idx_rand = np.random.permutation(n_left)

        for idx in idx_rand:
            ppr_bin = ppr_idx[idx]
            rtr_bin = rtr_idx[idx]
            ppr_p = ppr_count[ppr_bin] / n_left / ppr_bin_width
            rtr_p = rtr_count[rtr_bin] / n_left / rtr_bin_width

            # if(pprP > pprDensity[pprBin] or rtrP > rtrDensity[rtrBin] \
            #   or pprBin == 0 or rtrBin == 0):
            if (ppr_p / ppr_density[ppr_bin] + rtr_p / rtr_density[rtr_bin] > 2
                    or ppr_bin == 0 or rtr_bin == 0
                    or ppr_density[ppr_bin] < 1e-4 or rtr_density[rtr_bin] < 1e-4):
                # This point is over represented, lets remove it
                keep_flag[idx] = False
                n_left -= 1
                ppr_count[ppr_idx[idx]] -= 1
                rtr_count[rtr_idx[idx]] -= 1

        pars["u"] = pars["u"][keep_flag]
        pars["ftau"] = pars["ftau"][keep_flag]
        pars["dtau"] = pars["dtau"][keep_flag]
        pars["amp"] = pars["amp"][keep_flag]

        print(f"(Peaks) Reduced {len(keep_flag)} down to {sum(keep_flag)}")

        return pars, amps[keep_flag, :]

    ############################################################################

    @staticmethod
    def max_num_params(pars, amps, max_num=100):

        n_par_sets = len(pars["u"])
        if n_par_sets < max_num:
            return pars, amps

        print(f"Reducing from {n_par_sets} down to {max_num} parameter sets")

        pars["u"] = pars["u"][:max_num]

        pars["ftau"] = pars["ftau"][:max_num]
        pars["dtau"] = pars["dtau"][:max_num]
        pars["amp"] = pars["amp"][:max_num]

        return pars, amps[:max_num, :]

        ############################################################################

    # This prunes the parameters to get back their distribution to match

    def prune_parameters(self, par_type, pars, amps):

        (d_tau_mean, d_tau_std) = self.con_info[par_type]["dtau"]
        (f_tau_mean, f_tau_std) = self.con_info[par_type]["ftau"]
        (u_mean, u_std) = self.con_info[par_type]["U"]

        n_bins = 20

        u_max = np.max(pars["u"])
        d_tau_max = np.max(pars["dtau"])
        f_tau_max = np.max(pars["ftau"])

        u_idx = (n_bins * pars["u"] / (u_max + 1e-6)).astype(int)
        d_tau_idx = (n_bins * pars["dtau"] / (d_tau_max + 1e-6)).astype(int)
        f_tau_idx = (n_bins * pars["ftau"] / (f_tau_max + 1e-6)).astype(int)

        u_count = np.zeros((n_bins,))
        d_tau_count = np.zeros((n_bins,))
        f_tau_count = np.zeros((n_bins,))

        for u, d, f in zip(u_idx, d_tau_idx, f_tau_idx):
            u_count[u] += 1
            d_tau_count[d] += 1
            f_tau_count[f] += 1

        u_bin_width = u_max / (n_bins - 1)
        d_tau_bin_width = d_tau_max / (n_bins - 1)
        f_tau_bin_width = f_tau_max / (n_bins - 1)

        u_centre = u_bin_width / 2 + np.arange(0, n_bins) * u_bin_width
        d_tau_centre = d_tau_bin_width / 2 + np.arange(0, n_bins) * d_tau_bin_width
        f_tau_centre = f_tau_bin_width / 2 + np.arange(0, n_bins) * f_tau_bin_width

        u_density = stats.norm.pdf(u_centre, u_mean, u_std)
        d_tau_density = stats.norm.pdf(d_tau_centre, d_tau_mean * 1e-3, d_tau_std * 1e-3)
        f_tau_density = stats.norm.pdf(f_tau_centre, f_tau_mean * 1e-3, f_tau_std * 1e-3)

        n_left = len(pars["u"])
        keep_flag = np.ones((n_left,), dtype=bool)

        idx_rand = np.random.permutation(n_left)

        for idx in idx_rand:
            u_bin = u_idx[idx]
            d_tau_bin = d_tau_idx[idx]
            f_tau_bin = f_tau_idx[idx]

            u_p = u_count[u_bin] / n_left / u_bin_width
            d_tau_p = d_tau_count[d_tau_bin] / n_left / d_tau_bin_width
            f_tau_p = f_tau_count[f_tau_bin] / n_left / f_tau_bin_width

            if (u_p / u_density[u_bin] + d_tau_p / d_tau_density[d_tau_bin] + f_tau_p / f_tau_density[f_tau_bin] > 3
                    or u_density[u_bin] < 1e-4 or d_tau_density[d_tau_bin] < 1e-4 or f_tau_density[f_tau_bin] < 1e-4):
                keep_flag[idx] = False
                n_left -= 1
                u_count[u_idx[idx]] -= 1
                d_tau_count[d_tau_idx[idx]] -= 1
                f_tau_count[f_tau_idx[idx]] -= 1

                # import pdb
                # pdb.set_trace()

        pars["u"] = pars["u"][keep_flag]
        pars["ftau"] = pars["ftau"][keep_flag]
        pars["dtau"] = pars["dtau"][keep_flag]
        pars["amp"] = pars["amp"][keep_flag]

        print(f"(Pars) Reduced {len(keep_flag)} down to {sum(keep_flag)}")

        return pars, amps[keep_flag, :]

    ############################################################################

    def check_data_match(self, par_type, pars, amps=None):

        ppr_mean, ppr_std = self.con_info[par_type]["ppr"]
        rtr_mean, rtr_std = self.con_info[par_type]["rtr"]

        if amps is None:
            amps = self.synapse_model_wrapper(self.t_stim, pars)

        ppr = np.divide(amps[:, 1], amps[:, 0])
        rtr = np.divide(amps[:, -1], amps[:, 0])

        # self.plotSpikes(amps,self.tSpike,pars,0)

        if True:
            self.plot_data_match(values=ppr,
                                 mean=ppr_mean,
                                 std=ppr_std,
                                 title=par_type,
                                 xlabel="PPR")

            self.plot_data_match(values=rtr,
                                 mean=rtr_mean,
                                 std=rtr_std,
                                 title=par_type,
                                 xlabel="RTR")

        return amps

    ############################################################################

    @staticmethod
    def plot_data_match(values, mean, std, title, xlabel, n_bins=20):

        print(f"{xlabel} Mean : {mean}, std : {std}")
        print(f"Data : {np.mean(values)}, std : {np.std(values)}")
        n_points = 100
        x = np.linspace(mean - 3 * std, mean + 3 * std, n_points)

        plt.figure()
        plt.plot(x, stats.norm.pdf(x, mean, std))
        plt.hist(values, density=True, bins=n_bins)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ion()
        plt.show()
        plt.tight_layout()

        fig_name = os.path.join("DATA", "Planert2010", "figs",
                                f"Planert2010-surrogate-data-pre-fit-{title}-{xlabel}.pdf")
        plt.savefig(fig_name)
        plt.close()

    ############################################################################

    @staticmethod
    def plot_spikes(amps, t_stim, pars=None, idx=0):

        plt.figure()
        for a, t in zip(amps[idx, :], t_stim):
            plt.plot([t, t], [0, a], color="black")
        plt.xlabel("Time")
        plt.ylabel("Amplitude")

        if pars is not None:
            title_str = "U = %.3f, ftau = %.3f, dtau = %.3f, amp = %.3f" \
                        % (pars["u"][idx], pars["ftau"][idx], pars["dtau"][idx], pars["amp"][idx])
            plt.title(title_str)

        plt.ion()
        plt.show()

    ############################################################################

    # parType = "DD", "DI", "ID", "II", "FD" or "FI"

    def pick_random_param(self, par_type, num):

        print(f"Picking {num} {par_type} parameter sets...")

        # Need to randomise amp, dtau, ftau and U

        amp_mean, amp_std = self.con_info[par_type]["amp"]
        d_tau_mean, d_tau_std = self.con_info[par_type]["dtau"]
        f_tau_mean, f_tau_std = self.con_info[par_type]["ftau"]
        u_mean, u_std = self.con_info[par_type]["U"]

        amp = (np.random.normal(size=num, scale=amp_std) + amp_mean) * 1e-3
        d_tau = (np.random.normal(size=num, scale=d_tau_std) + d_tau_mean) * 1e-3
        f_tau = (np.random.normal(size=num, scale=f_tau_std) + f_tau_mean) * 1e-3
        u = np.random.normal(size=num, scale=u_std) + u_mean

        # If any of the conditions are true, we are outside OK range,
        # so only OK if all conditions are false, and the result sums to False
        ok_idx = False == (np.array(amp < 0) + np.array(d_tau < 0)
                           + np.array(f_tau < 0) + np.array(u < 0) + np.array(u > 1))

        # Then we need to verify that the ppr and rtr are within range

        print(f"Npoints = {np.sum(ok_idx)}")

        pars = {"amp": amp[ok_idx],
                "dtau": d_tau[ok_idx],
                "ftau": f_tau[ok_idx],
                "u": u[ok_idx]}

        return pars

    ############################################################################

    def synapse_model_wrapper(self, t_stim, pars):

        print(f"Running {len(pars['u'])} simulations")

        amps = np.zeros((len(pars["u"]), len(t_stim)))

        for i, (u, d_tau, f_tau, Asc) in enumerate(zip(pars["u"],
                                                       pars["dtau"],
                                                       pars["ftau"],
                                                       pars["amp"])):
            if i > 0 and i % 100000 == 0:
                # Print progress...
                print(str(i) + "/" + str(amps.shape[0]))

            amps[i, :] = self.synapse_model(t_stim, U=u, d_tau=d_tau, f_tau=f_tau, asc=Asc)

        return amps

    ############################################################################

    @staticmethod
    def synapse_model(t_stim, U, d_tau, f_tau, asc):

        n_spikes = len(t_stim)
        u = np.zeros((n_spikes,))
        r = np.zeros((n_spikes,))
        isi = np.diff(t_stim)

        # Init
        u[0] = U
        r[0] = 1

        for idx, ii in enumerate(isi):
            u[idx + 1] = u[idx] * np.exp(-isi[idx] / f_tau) + U * (1 - u[idx] * np.exp(-isi[idx] / f_tau))

            r[idx + 1] = r[idx] * (1 - u[idx]) * np.exp(-isi[idx] / d_tau) + 1 - np.exp(-isi[idx] / d_tau)

        amp = np.multiply(u, r) * asc

        return amp

    ############################################################################

    def def_exp_data(self):

        self.con_info = dict()

        self.con_info["DD"] = {"P": (3, 43),
                               "amp": (0.24, 0.15),
                               "ppr": (0.91, 0.63),
                               "rtr": (1.23, 0.50),
                               "dtau": (192, 114),
                               "ftau": (1266, 1427),
                               "U": (0.39, 0.22),
                               "FDratio": 4.76}

        self.con_info["DI"] = {"P": (3, 66),
                               "amp": (0.33, 0.15),
                               "ppr": (0.84, 0.3),
                               "rtr": (1.16, 0.29),
                               "dtau": (96, 9),
                               "ftau": (313, 363),
                               "U": (0.46, 0.24),
                               "FDratio": 4.76}

        self.con_info["ID"] = {"P": (10, 80),
                               "amp": (0.27, 0.09),
                               "ppr": (1.1, 0.6),
                               "rtr": (1.51, 0.64),
                               "dtau": (365, 471),
                               "ftau": (570, 783),
                               "U": (0.36, 0.18),
                               "FDratio": 4.76}

        self.con_info["II"] = {"P": (7, 31),
                               "amp": (0.45, 0.44),
                               "ppr": (0.95, 0.48),
                               "rtr": (1.39, 0.69),
                               "dtau": (149, 90),
                               "ftau": (1462, 1800),
                               "U": (0.34, 0.19),
                               "FDratio": 4.76}

        self.con_info["FD"] = {"P": (8, 9),
                               "amp": (4.8, 4.9),
                               "ppr": (0.62, 0.12),
                               "rtr": (0.72, 0.08),
                               "dtau": (740, 350),
                               "ftau": (3.1, 2.4),
                               "U": (0.24, 0.07),
                               "FDratio": 0.16}

        self.con_info["FI"] = {"P": (6, 9),
                               "amp": (3.1, 4.1),
                               "ppr": (0.66, 0.14),
                               "rtr": (0.63, 0.19),
                               "dtau": (850, 500),
                               "ftau": (4.5, 2.7),
                               "U": (0.23, 0.07),
                               "FDratio": 0.16}

        self.t_stim = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.9]

    ############################################################################

    def save_data(self, par_type, pars, amp):
        file_name = os.path.join("DATA", "Planert2010", f"PlanertFitting-{par_type}-cache.json")

        data = dict()

        for p in pars:
            data[p] = pars[p]

        data["simAmp"] = amp
        data["tStim"] = self.t_stim

        with open(file_name, "w") as f:
            json.dump(data, f, indent=2, cls=NumpyEncoder)

    ############################################################################

    @staticmethod
    def load_data(par_type):

        file_name = os.path.join("Data", "Planert2010", f"PlanertFitting-{par_type}-cache.json")

        if not os.path.exists(file_name):
            return None, None

        with open(file_name, "r") as f:
            data = json.load(f)

        pars = dict([])
        for p in data:
            if p in ["u", "dtau", "ftau", "amp"]:
                pars[p] = data[p]
            elif p != "simAmp":
                print(f"Unknown parameter : {p}")

        amps = data["simAmp"]

        return pars, amps

    ############################################################################

    def plot_model_params(self, data_type, pars):

        (d_tau_mean, d_tau_std) = self.con_info[data_type]["dtau"]
        (f_tau_mean, f_tau_std) = self.con_info[data_type]["ftau"]
        (u_mean, u_std) = self.con_info[data_type]["U"]

        self._plot_model_params_helper(data_type, d_tau_mean * 1e-3, d_tau_std * 1e-3, pars["dtau"], "dtau")
        self._plot_model_params_helper(data_type, f_tau_mean * 1e-3, f_tau_std * 1e-3,
                                       pars["ftau"], "ftau")
        self._plot_model_params_helper(data_type, u_mean, u_std, pars["u"], "U")

        ############################################################################

    @staticmethod
    def _plot_model_params_helper(data_type, data_mean, data_std, data_points, x_label):

        x_min = np.min(data_points)
        x_max = np.max(data_points)

        x = np.linspace(x_min, x_max, 100)
        x_density = stats.norm.pdf(x, data_mean, data_std)  # / (x[1]-x[0])

        plt.figure()
        plt.hist(data_points, bins=20, density=True)
        plt.plot(x, x_density)
        plt.xlabel(x_label)
        plt.ylabel("Density")

        plt.tight_layout()

        fig_name = os.path.join("DATA", "Planert2010", "figs",
                                f"Surrogate-variables-distribution-{data_type}-{x_label}.pdf")
        plt.savefig(fig_name)
        plt.close()

    ############################################################################


if __name__ == "__main__":
    n_runs = 1000000
    max_num = 50

    pp = Planert2010(par_type="FI", pre_num=n_runs, max_num=max_num)
    pp = Planert2010(par_type="FD", pre_num=n_runs, max_num=max_num)

    pp = Planert2010(par_type="DD", pre_num=n_runs, max_num=max_num)
    pp = Planert2010(par_type="DI", pre_num=n_runs, max_num=max_num)
    pp = Planert2010(par_type="ID", pre_num=n_runs, max_num=max_num)
    pp = Planert2010(par_type="II", pre_num=n_runs, max_num=max_num)


