TITLE SK-type calcium activated K channel (KCa2.2)

UNITS {
    (molar) = (1/liter)
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (millimolar)
}

NEURON {
    SUFFIX sk_ms
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gbar, ik
}

PARAMETER {
    gbar = 0.0 	(mho/cm2)
    :q = 1	: room temperature
    q = 1	: body temperature
}

ASSIGNED {
    v (mV)
    ik (mA/cm2)
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
    LOCAL a
    a = (ca/0.57e-3)^5.2
    oinf = a/(1+a)
    otau = 4.9
}

COMMENT

Experimental data was obtained for the apamin-sensitive clone rSK2 from
rat brain cDNA expressed in Xenopus oocytes [1,2].  All experiments were
performed at room tempretaure.

Original model [3] used calcium dependence from [2, Fig.2] and calcium
activation time constant from [1,  Fig.13]. NEURON implementation by
Alexander Kozlov <akozlov@kth.se> follows the revised model [4].

[1] Hirschberg B, Maylie J, Adelman JP, Marrion NV (1998) Gating of
recombinant small-conductance Ca-activated K+ channels by calcium. J
Gen Physiol 111(4):565-81.

[2] Maylie J, Bond CT, Herson PS, Lee WS, Adelman JP (2004) Small
conductance Ca2+-activated K+ channels and calmodulin. J Physiol 554(Pt
2):255-61.

[3] Evans RC, Morera-Herreras T, Cui Y, Du K, Sheehan T, Kotaleski JH,
Venance L, Blackwell KT (2012) The effects of NMDA subunit composition on
calcium influx and spike timing-dependent plasticity in striatal medium
spiny neurons. PLoS Comput Biol 8(4):e1002493.

[4] Evans RC, Maniar YM, Blackwell KT (2013) Dynamic modulation of
spike timing-dependent calcium influx during corticostriatal upstates. J
Neurophysiol 110(7):1631-45.

ENDCOMMENT
