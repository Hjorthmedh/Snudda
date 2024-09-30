TITLE Non-inactivating inwardly rectifying potassium current (Kir2.3)


NEURON {
    SUFFIX kirrxd
    USEION k READ ek WRITE ik
    USEION PKA READ PKAi VALENCE 0

    RANGE gbar, gk, ik, shift
    RANGE  mod_pka_g_min, mod_pka_g_max, mod_pka_g_half, mod_pka_g_slope		     
    RANGE modulation_factor							   
}

UNITS {
    (S) = (siemens)
    (mV) = (millivolt)
    (mA) = (milliamp)
    (molar) =  (1/liter)
    (mM) = (millimolar)
}

PARAMETER {
    gbar = 0.0 	(S/cm2) 
    shift = 0.0 (mV)
    q = 1 	: body temperature 35 C
    mod_pka_g_min = 1 (1)
    mod_pka_g_max = 1 (1)
    mod_pka_g_half = 0.000100 (mM)
    mod_pka_g_slope = 0.01 (mM)
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    gk (S/cm2)
    minf
    mtau (ms)
    PKAi (mM)
    modulation_factor (1)
}

STATE { m }

BREAKPOINT {
    SOLVE states METHOD cnexp
    modulation_factor=modulation(PKAi, mod_pka_g_min, mod_pka_g_max, mod_pka_g_half, mod_pka_g_slope)	   
    gk = gbar*m*modulation_factor
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

: TODO: These parameters should NOT be hardcoded, and use function instead of procedure?
PROCEDURE rates() {
    UNITSOFF
    minf = 1/(1+exp((v-(-82)-shift)/13))
    mtau = 1/(exp((v-(-103))/(-14.5))+0.125/(1+exp((v-(-35))/(-19))))
    UNITSON
}

FUNCTION modulation(conc (mM), mod_min (1), mod_max (1), mod_half (mM), mod_slope (mM)) (1) {
    : returns modulation factor
    modulation = mod_min + (mod_max-mod_min) / (1 + exp(-(conc - mod_half)/mod_slope))
}


COMMENT

2024-05-28 : Wilhelm Thunberg, Johannes Hjorth (KTH, Stockholm)
Adding neuromodulation using RxD

Original model by Wolf et al (2005) [1] for the rat MSN cells from the
nucleus accumbens.  The activation curve was fitted to a mouse Kir2.1
channel expressed in HEK cells [2] and shifted to match extracellular
concentration of K in rat. Measured half-activation values are -109.3
mV (striatonigral MSN) and -113.2 mV (striatopallidal MSN) [6, Supp
Tab.1]. Time constants were derived from Aplysia data [3] and adjusted
to match the rat experiments [1]. Time constant was further tuned [4]
to fit the rat data below -80 mV [5].  Kinetics is corrected to the body
temperature 35 C [4].

Non-inactivating Kir current was observed in cells expressing Kir2.2
and/or Kir2.3 [5]. Activation variable with m^1 kinetics is used [1,4].
Smooth fit of the time constants by Alexander Kozlov <akozlov@kth.se>.

[1] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.

[2] Kubo Y, Murata Y (2001) Control of rectification and permeation by
two distinct sites after the second transmembrane region in Kir2.1 K+
channel. J Physiol 531, 645-660.

[3] Hayashi H, Fishman HM (1988) Inward rectifier K+ channel kinetics
from analysis of the complex conductance of aplysia neuronal membrane.
Biophys J 53, 747-757.

[4] Steephen JE, Manchanda R (2009) Differences in biophysical properties
of nucleus accumbens medium spiny neurons emerging from inactivation of
inward rectifying potassium currents. J Comput Neurosci 27(3):453-70

[5] Mermelstein PG, Song WJ, Tkatch T, Yan Z, Surmeier DJ (1998)
Inwardly rectifying potassium (IRK) currents are correlated with IRK
subunit expression in rat nucleus accumbens medium spiny neurons. J
Neurosci 18(17):6650-61.

[6] Shen W, Tian X, Day M, Ulrich S, Tkatch T, Nathanson NM, Surmeier DJ
(2007) Cholinergic modulation of Kir2 channels selectively elevates
dendritic excitability in striatopallidal neurons. Nat Neurosci
10(11):1458-66.

ENDCOMMENT
