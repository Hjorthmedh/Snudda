TITLE Fast delayed rectifier potassium current (Kv3.1/3.2)

NEURON {
    SUFFIX kdr_ms
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 (S/cm2) 
    q = 1
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    gk (S/cm2)
    minf
    mtau (ms)
}

STATE { m }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gbar*m
    ik = gk*(v-ek)
}

DERIVATIVE states {
    rates()
    m' = (minf-m)/mtau*q
}

INITIAL {
    rates()
    m = minf
}

PROCEDURE rates() {
    UNITSOFF
    minf = 1/(1+exp((v-(-13))/(-6)))
    mtau = 11.1
    UNITSON
}

COMMENT

Original data and model by Baranauskas et al (1999) [1] for globus
pallidus neurons from young adult rat.  The temperature was not
specified. Potentials were not corrected for the liquid junction
potential, which was estimated to be 1-2 mV.

Kinetics of m^1 type is used as in [2,3]. Room temperature 20-23 C
is assumed.

NEURON implementation by Alexander Kozlov <akozlov@kth.se>.

[1] Baranauskas G, Tkatch T, Surmeier DJ (1999) Delayed rectifier currents
in rat globus pallidus neurons are attributable to Kv2.1 and Kv3.1/3.2
K(+) channels. J Neurosci 19(15):6394-404.

[2] Migliore M, Hoffman DA, Magee JC, Johnston D (1999) Role of an
A-type K+ conductance in the back-propagation of action potentials in the
dendrites of hippocampal pyramidal neurons. J Comput Neurosci 7(1):5-15.

[3] Evans RC, Morera-Herreras T, Cui Y, Du K, Sheehan T, Kotaleski JH,
Venance L, Blackwell KT (2012) The effects of NMDA subunit composition on
calcium influx and spike timing-dependent plasticity in striatal medium
spiny neurons. PLoS Comput Biol 8(4):e1002493.

ENDCOMMENT
