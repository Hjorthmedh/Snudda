TITLE AMPA synapse with short-term plasticity

NEURON {
    POINT_PROCESS tmAmpa
    RANGE tau1, tau2, e, g, i, q
    RANGE tau, tauR, tauF, U, u0
    RANGE ca_ratio, ical
    NONSPECIFIC_CURRENT i
    USEION cal WRITE ical VALENCE 2
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    : q = 2 correction applied
    tau1= 1.1 (ms) : ORIG: 2.2ms
    tau2 = 5.75 (ms)  : ORIG: 11.5 ms, tau2 > tau1
    e = 0 (mV)
    tau = 3 (ms)
    tauR = 100 (ms)   : tauR > tau
    tauF = 800 (ms)   : tauF >= 0
    U = 0.3 (1) <0, 1>
    u0 = 0 (1) <0, 1>
    ca_ratio = 0.005

}

ASSIGNED {
    v (mV)
    i (nA)
    ical (nA)
    g (uS)
    factor
    x
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
    LOCAL itotal
    SOLVE state METHOD cnexp
    g = B - A
    itotal = g*(v - e)
    ical = ca_ratio*itotal
    i = itotal - ical
}

DERIVATIVE state {
    A' = -A/tau1
    B' = -B/tau2
}

NET_RECEIVE(weight (uS), y, z, u, tsyn (ms)) {
    INITIAL {
        y = 0
        z = 0
        u = u0
        tsyn = t
    }
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
    A = A + weight*factor*x*u
    B = B + weight*factor*x*u
    y = y + x*u
    tsyn = t
}

COMMENT

Implementation of AMPA synapse model with short-term facilitation
and depression based on modified tmgsyn.mod [1] by Tsodyks et al [2].
Choice of time constants and calcium current model follows [3].
NEURON implementation by Alexander Kozlov <akozlov@kth.se>.

[1] tmgsyn.mod, ModelDB (https://senselab.med.yale.edu/ModelDB/),
accession number 3815.

[2] Tsodyks M, Uziel A, Markram H (2000) Synchrony generation in recurrent
networks with frequency-dependent synapses. J Neurosci. 20(1):RC50.

[3] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.

ENDCOMMENT
