TITLE BK-type calcium activated K channel (KCa1.1)

UNITS {
    (molar) = (1/liter)
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (millimolar)
    FARADAY = (faraday) (kilocoulombs)
    R = (k-mole) (joule/degC)
}

NEURON {
    SUFFIX bk_fs
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gbar, ik
}

PARAMETER {
    gbar = 0.0 (mho/cm2)
    k1 = 0.003 (mM)
    k4 = 0.009 (mM)
    d1 = 0.84
    d4 = 1.0
    q = 1	: body temperature 35 C
}

ASSIGNED {
    v (mV)
    ik (mA/cm2)
    celsius (degC)
    cai (mM) 
    ek (mV)
    oinf
    otau (ms)
}

STATE { o }

BREAKPOINT {
    SOLVE state METHOD cnexp
    ik = gbar*o*(v-ek)
}

DERIVATIVE state {
    rate(v, cai)
    o' = (oinf-o)/otau*q
}

INITIAL {
    rate(v, cai)
    o = oinf
}

PROCEDURE rate(v (mV), ca (mM)) {
    LOCAL a, b, sum, z
    UNITSOFF
    z = 1e-3*2*FARADAY/(R*(celsius+273.15))
    a = 0.48*ca/(ca+k1*exp(-z*d1*v))
    b = 0.28/(1+ca/(k4*exp(-z*d4*v)))
    sum = a+b
    oinf = a/sum
    otau = 1/sum
    UNITSON
}

COMMENT

Experimental data was obtained from BKCa channels from rat brain injected
as cRNAs into Xenopus oocytes [1].  Electrophysiological recordings were
performed at room temperature 22-24 C [1, supporting online material].

Original model [2, model 3 in Tab.1] was implemented by De Schutter
and adapted by Kai Du.  In the model revisions [3,4] parameters k1 and k4
[2, channel A in Tab.2] were adjusted to fit rat/Xenopus data
[1, Fig.3C and Fig.4A, 10 uM Ca] at body temperature 35 C.

NEURON implementation by Alexander Kozlov <akozlov@kth.se>.

[1] Berkefeld H, Sailer CA, Bildl W, Rohde V, Thumfart JO, Eble S,
Klugbauer N, Reisinger E, Bischofberger J, Oliver D, Knaus HG, Schulte U,
Fakler B (2006) BKCa-Cav channel complexes mediate rapid and localized
Ca2+-activated K+ signaling. Science 314(5799):615-20.

[2] Moczydlowski E, Latorre R (1983) Gating kinetics of Ca2+-activated K+
channels from rat muscle incorporated into planar lipid bilayers. Evidence
for two voltage-dependent Ca2+ binding reactions. J Gen Physiol
82(4):511-42.

[3] Evans RC, Morera-Herreras T, Cui Y, Du K, Sheehan T, Kotaleski JH,
Venance L, Blackwell KT (2012) The effects of NMDA subunit composition on
calcium influx and spike timing-dependent plasticity in striatal medium
spiny neurons. PLoS Comput Biol 8(4):e1002493.

[4] Evans RC, Maniar YM, Blackwell KT (2013) Dynamic modulation of
spike timing-dependent calcium influx during corticostriatal upstates. J
Neurophysiol 110(7):1631-45.

ENDCOMMENT
