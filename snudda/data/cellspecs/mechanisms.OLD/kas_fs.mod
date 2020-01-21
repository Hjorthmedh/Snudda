TITLE Slowly inactivating A-type potassium current (Kv1.2)

NEURON {
    SUFFIX kas_fs
    USEION k READ ek WRITE ik
    RANGE gbar, gk, ik, shift, q
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
}

PARAMETER {
    gbar = 0.0 	(S/cm2) 
    a = 0.8
    :q = 1	: room temperature 22-24 C
    q = 3	: body temperature 33 C
    shift = 0
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    gk (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk = gbar*m*m*(h*a+1-a)
    ik = gk*(v-ek)
}

DERIVATIVE states {
    rates()
    m' = (minf-m)/mtau*q
    h' = (hinf-h)/htau*q
}

INITIAL {
    rates()
    m = minf
    h = hinf
}

PROCEDURE rates() {
    UNITSOFF
    minf = 1/(1+exp((v-(-27)-shift)/(-16)))
    mtau = 3.4+89.2*exp(-((v-(-34.3))/30.1)^2)
    hinf = 1/(1+exp((v-(-33.5)-shift)/21.5))
    htau = 548.7*6/(exp((v-(-96))/(-29.01))+exp((v-(-96))/100))
    UNITSON
}

COMMENT

Experimental data by Shen et al (2004) [1]. Medium spiny neurons were
acutely dissociated from from young adult (P21-P28) Sprague-Dawley rat
brain. All recordings were conducted at 22-24 C. No correction for the
liquid junction potential was reported.

Conductance kinetics of m2h type is used [1] with partial inactivation,
m2 (a h + (1-a)). Fraction a is set to 0.8, as in [1, Fig.6B]; other
values for a are possible [2] (see also kas.mod in companion code).
Equation for htau [1] is corrected to match the authors' data [1, Fig.6B]
by Alexander Kozlov <akozlov@kth.se>.  Time constants were corrected to
body temperature with factor q=3 [1,3].

[1] Shen W, Hernandez-Lopez S, Tkatch T, Held JE, Surmeier DJ (2004)
Kv1.2-containing K+ channels regulate subthreshold excitability of
striatal medium spiny neurons. J Neurophysiol 91(3):1337-49.

[2] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.

[3] Evans RC, Maniar YM, Blackwell KT (2013) Dynamic modulation of
spike timing-dependent calcium influx during corticostriatal upstates. J
Neurophysiol 110(7):1631-45.

ENDCOMMENT
