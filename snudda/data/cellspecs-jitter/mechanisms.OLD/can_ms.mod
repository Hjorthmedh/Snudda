TITLE N-type calcium current (Cav2.2)

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

NEURON {
    SUFFIX can_ms
    USEION ca READ cai, cao WRITE ica VALENCE 2
    RANGE pbar, ica
}

PARAMETER {
    pbar = 0.0 	(cm/s)
    a = 0.21
    :q = 1	: room temperature 22-25 C
    q = 3	: body temperature 35 C
} 

ASSIGNED { 
    v (mV)
    ica (mA/cm2)
    eca (mV)
    celsius (degC)
    cai (mM)
    cao (mM)
    minf
    mtau (ms)
    hinf
    htau (ms)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica = pbar*m*m*(h*a+1-a)*ghk(v, cai, cao)
}

INITIAL {
    rates()
    m = minf
    h = hinf
}

DERIVATIVE states { 
    rates()
    m' = (minf-m)/mtau*q
    h' = (hinf-h)/htau*q
}

PROCEDURE rates() {
    UNITSOFF
    minf = 1/(1+exp((v-(-3))/(-8)))
    mtau = 0.06+1/(exp((v-25)/18)+exp((v-(-31))/(-44)))
    hinf = 1/(1+exp((v-(-74.8))/6.5))
    htau = 70
    UNITSON
}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco
    z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
    eco = co*efun(z)
    eci = ci*efun(-z)
    ghk = (1e-3)*2*FARADAY*(eci-eco)
}

FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
        efun = 1-z/2
    }else{
        efun = z/(exp(z)-1)
    }
}

COMMENT

Model is based on mixed data. Activation curve is from neostriatal
medium spiny neurons of adult P28+ rats [1, Fig.12F], unspecified
recording temperature. Potentials were not corrected for the liquid
junction potential, which was estimated to be 7 mV.  Activation time
constant is from the rodent neuron culture (both rat and mouse cells),
room temperature 22-25 C [2, Fig.15B].  Inactivation data is from human
(HEK) cells [3, Tab.1, Tab.2], supposedly at room temperature.

Kinetics of m2h type is used [2, Fig.5]. Activation of m^2 is fitted
to the experimental data [1,4], activation time constant corresponds
to m^2 already [2].  Original model [5,4] was modified by Alexander
Kozlov <akozlov@kth.se>. Activation time constant was refitted to avoid
singularity in the expression. Temperature correction factor 3 is used
for body temperature [4,5].

[1] Bargas J, Howe A, Eberwine J, Cao Y, Surmeier DJ (1994) Cellular
and molecular characterization of Ca2+ currents in acutely isolated,
adult rat neostriatal neurons. J Neurosci 14(11 Pt 1):6667-86.

[2] Kasai H, Neher E (1992) Dihydropyridine-sensitive and
omega-conotoxin-sensitive calcium channels in a mammalian
neuroblastoma-glioma cell line. J Physiol 448:161-88.

[3] McNaughton NC, Randall AD (1997) Electrophysiological properties of
the human N-type Ca2+ channel: I. Channel gating in Ca2+, Ba2+ and Sr2+
containing solutions. Neuropharmacology 36(7):895-915.

[4] Evans RC, Maniar YM, Blackwell KT (2013) Dynamic modulation of
spike timing-dependent calcium influx during corticostriatal upstates. J
Neurophysiol 110(7):1631-45.

[5] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.

ENDCOMMENT
