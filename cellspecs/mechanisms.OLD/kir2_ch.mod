:Kir2_ch.MOD
: Kir2, inwardly rectifying channel

NEURON {
	SUFFIX kir2_ch
	USEION k READ ek WRITE ik
	RANGE g, ninf, tn, ik, gbar
	GLOBAL C_tn, vh, vc
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 3	(S/cm2)
	ek		(mV)
	vh = -80	(mV)
	vc = 5		(mV)
	C_tn = 1	(ms)
}

ASSIGNED {
	v	(mV)
	ninf
	tn	(ms)
	ik	(mA/cm2)
	g	(S/cm2)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*n
	ik = g*(v-ek)
}

DERIVATIVE states{
	values()
	n' = (ninf - n)/tn
}

INITIAL {
	values()
	n = ninf
}

PROCEDURE values() {
	ninf = 1/(1 + exp((v - vh)/vc))
	tn = C_tn
}