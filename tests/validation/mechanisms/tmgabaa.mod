TITLE GABA_A synapse with short-term plasticity


NEURON {
    THREADSAFE
    POINT_PROCESS tmGabaA
    RANGE tau1, tau2, e, i, q
    RANGE tau, tauR, tauF, U, u0
    RANGE failRate
    NONSPECIFIC_CURRENT i

    USEION PKAc READ PKAci VALENCE 0
    RANGE mod_pka_g_min, mod_pka_g_max, mod_pka_g_half, mod_pka_g_slope
    RANGE mod_pka_fail_min, mod_pka_fail_max, mod_pka_fail_half, mod_pka_fail_slope
    RANGE modulation_factor, modulation_factor_fail

    RANGE g
    RANDOM release_probability
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    : q = 2, now included in tau1,tau2 parameters.
    tau1= 0.25 (ms) : ORIG: 0.5ms
    tau2 = 3.75 (ms)  : ORIG: 7.5ms, tau2 > tau1
    e = -60 (mV)
    tau = 3 (ms)
    tauR = 500 (ms)  : tauR > tau
    tauF = 0 (ms)    : tauF >= 0
    U = 0.1 (1) <0, 1>
    u0 = 0 (1) <0, 1>
    mod_pka_g_min = 1 (1)
    mod_pka_g_max = 1 (1)
    mod_pka_g_half = 0.000100 (mM)
    mod_pka_g_slope = 0.01 (mM)
    mod_pka_fail_min = 0 (1)
    mod_pka_fail_max = 0 (1)
    mod_pka_fail_half = 0.000100 (mM)
    mod_pka_fail_slope = 0.01 (mM)
    failRate = 0
    failRateModulationScaling = 0
}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
    factor
    x
    PKAci (mM)

    modulation_factor (1)
    modulation_factor_fail (1)

}

STATE {
    A (uS)
    B (uS)
}

INITIAL {
    LOCAL tp
    A = 0
    B = 0
    tp = (tau1*tau2)/(tau2-tau1) * log(tau2/tau1)
    factor = -exp(-tp/tau1) + exp(-tp/tau2)
    factor = 1/factor
}

BREAKPOINT {
     SOLVE state METHOD cnexp
     modulation_factor=modulation(PKAci, mod_pka_g_min, mod_pka_g_max, mod_pka_g_half, mod_pka_g_slope)
     modulation_factor_fail=modulation(PKAci, mod_pka_fail_min, mod_pka_fail_max, mod_pka_fail_half, mod_pka_fail_slope)
    g = (B - A)*modulation_factor
    i = g*(v - e)
}

DERIVATIVE state {
    A' = -A/tau1
    B' = -B/tau2
}

NET_RECEIVE(weight (uS), y, z, u, tsyn (ms)) {
    LOCAL result
    INITIAL {
        y = 0
        z = 0
        u = u0
        tsyn = t
    }
    if ( weight <= 0 ) {
        VERBATIM
        return;
        ENDVERBATIM
    }
    if( urand() > failRate*(1 + modulation_factor_fail)) {
      z = z*exp(-(t-tsyn)/tauR)
      z = z + (y*(exp(-(t-tsyn)/tau) - exp(-(t-tsyn)/tauR)) / (tau/tauR - 1) )
      y = y*exp(-(t-tsyn)/tau)
      x = 1-y-z
      if (tauF > 0) {
          u = u*exp(-(t-tsyn)/tauF)
          u = u + U*(1-u)
      } else {
          u = U
      }
    A = A + weight*factor*x*u / U
    B = B + weight*factor*x*u / U
    y = y + x*u
    tsyn = t
    }
}

FUNCTION urand() {
    urand = random_uniform(release_probability)
}

FUNCTION modulation(conc (mM), mod_min (1), mod_max (1), mod_half (mM), mod_slope (mM)) (1) {
    : returns modulation factor
    modulation = mod_min + (mod_max-mod_min) / (1 + exp(-(conc - mod_half)/mod_slope))
}

COMMENT
(2025-10-08) NEURON 9.0+ compatibility. Replaced scop_random with the
new RANDOM keyword.
See: https://nrn.readthedocs.io/en/latest/nmodl/language/nmodl_neuron_extension.html#random

(2019-11-25) Synaptic failure rate (failRate) added. Random factor, no
reproducibility guaranteed in parallel sim.

(2019-09-12) Set GABA reversal potential to -65mV

(2019-08-21) Normalise activation by U, to make sure first activation has
amplitude set by g

(2019-06-05) Q-factor was calculated in INITAL block, which meant if
the synapse was reinitalised then the time constants changed with each
initalise. Updated: Johannes Hjorth, hjorth@kth.se

(2025-04-02) Set GABA reversal potential to -60mV as per [4]

Implementation of GABA_A synapse model with short-term facilitation
and depression based on modified tmgsyn.mod [1] by Tsodyks et al [2].
Choice of time constants follows [3].  NEURON implementation by Alexander
Kozlov <akozlov@kth.se>.

[1] tmgsyn.mod, ModelDB (https://senselab.med.yale.edu/ModelDB/),
accession number 3815.

[2] Tsodyks M, Uziel A, Markram H (2000) Synchrony generation in recurrent
networks with frequency-dependent synapses. J Neurosci. 20(1):RC50.

[3] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.

[4] Day M, Belal M, Surmeier WC, Melendez A, Wokosin D, Tkatch T,
Clarke VRJ, Surmeier DJ (2024) GABAergic regulation of striatal spiny
projection neurons depends upon their activity state. PLOS Biology.
22(7):e3002752.
ENDCOMMENT
