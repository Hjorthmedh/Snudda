: HH P-type Calcium current
: Created 8/13/02 - nwg

: copy by josh for cholinergic interneuron

NEURON {
	SUFFIX cap_ch
	USEION ca READ cai, cao WRITE ica
	RANGE gbar, ica ,g
	GLOBAL minf,mtau
	GLOBAL monovalConc, monovalPerm
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(mM) = (milli/liter)
	F = 9.6485e4   (coul)
	R = 8.3145 (joule/degC)
}

PARAMETER {
	v (mV)

	gbar = .00005	(cm/s)
	monovalConc = 140     (mM)
	monovalPerm = 0
	celsius = 35
	cai             (milli/liter)
	cao             (milli/liter)
}

ASSIGNED {
	ica            (mA/cm2)
        minf
	mtau           (ms)
	T              (degC)
	E              (volts)
	g	(S/cm2)
}

STATE {
	m
}

INITIAL {
	rates(v)
	m = minf
}

BREAKPOINT {
     SOLVE states METHOD cnexp
	g = (1e3) * gbar * m 
	ica = g * ghk(v, cai, cao, 2)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (coul/cm3) { LOCAL Ci
	T = celsius + 273.19  : Kelvin
        E = (1e-3) * v
        Ci = ci + (monovalPerm) * (monovalConc)        : Monovalent permeability
	if (fabs(1-exp(-z*(F*E)/(R*T))) < 1e-6) { : denominator is small -> Taylor series
		ghk = (1e-6) * z * F * (Ci-co*exp(-z*(F*E)/(R*T)))*(1-(z*(F*E)/(R*T)))
	} else {
		ghk = (1e-6) * z^2*(E*F^2)/(R*T)*(Ci-co*exp(-z*(F*E)/(R*T)))/(1-exp(-z*(F*E)/(R*T)))
	}
}

PROCEDURE rates (v (mV)) {
        UNITSOFF
	minf = 1/(1+exp(-(v - (-19)) / 5.5))
	mtau = (mtau_func(v)) * 1e3
        UNITSON
}

FUNCTION mtau_func( v (mV) ) (ms) {
        UNITSOFF
        if (v > -50) {
            mtau_func = .000191 + .00376*exp(-((v-(-41.9))/27.8)^2)
        } else {
            mtau_func = .00026367 + .1278 * exp(.10327*v)
        }
        UNITSON
}
