TITLE HVA L-type calcium current (Cav1.2)

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
    SUFFIX cal12_ms
    USEION cal READ cali, calo WRITE ical VALENCE 2
    RANGE pbar, ical
}

PARAMETER {
    pbar = 0.0 (cm/s)
    a = 0.17
    :q = 1	          : room temperature 22-25 C
    q = 2	          : body temperature 35 C
} 

ASSIGNED { 
    v (mV)
    ical (mA/cm2)
    ecal (mV)
    celsius (degC)
    cali (mM)
    calo (mM)
    minf
    mtau (ms)
    hinf
    htau (ms)
}

STATE { m h }

BREAKPOINT {
    SOLVE states METHOD cnexp
    ical = pbar*m*(h*a+1-a)*ghk(v, cali, calo)
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
    minf = 1/(1+exp((v-(-8.9))/(-6.7)))
    mtau = 0.06+1/(exp((v-10)/20)+exp((v-(-17))/-48))
    hinf = 1/(1+exp((v-(-13.4))/11.9))
    htau = 44.3
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

Activation curve was reconstructed for cultured NAc neurons from P5-P32
Charles River rat pups [1].   Activation time constant is from the
rodent neuron culture (both rat and mouse cells), room temperature 22-25
C [2, Fig.15A]. Inactivation curve of CaL v1.3 current was taken from HEK
cells [3, Fig.2 and p.819] at room temperature.

Original NEURON model by Wolf (2005) [4] was modified by Alexander Kozlov
<akozlov@csc.kth.se>. Kinetics of m1h type was used [5,6]. Activation
time constant was refitted to avoid singularity.

[1] Churchill D, Macvicar BA (1998) Biophysical and pharmacological
characterization of voltage-dependent Ca2+ channels in neurons isolated
from rat nucleus accumbens. J Neurophysiol 79(2):635-47.

[2] Kasai H, Neher E (1992) Dihydropyridine-sensitive and
omega-conotoxin-sensitive calcium channels in a mammalian
neuroblastoma-glioma cell line. J Physiol 448:161-88.

[3] Bell DC, Butcher AJ, Berrow NS, Page KM, Brust PF, Nesterova A,
Stauderman KA, Seabrook GR, Nurnberg B, Dolphin AC (2001) Biophysical
properties, pharmacology, and modulation of human, neuronal L-type
(alpha(1D), Ca(V)1.3) voltage-dependent calcium currents. J Neurophysiol
85:816-827.

[4] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.

[5] Evans RC, Morera-Herreras T, Cui Y, Du K, Sheehan T, Kotaleski JH,
Venance L, Blackwell KT (2012) The effects of NMDA subunit composition on
calcium influx and spike timing-dependent plasticity in striatal medium
spiny neurons. PLoS Comput Biol 8(4):e1002493.

[6] Tuckwell HC (2012) Quantitative aspects of L-type Ca2+ currents. Prog
Neurobiol 96(1):1-31.

ENDCOMMENT
