NEURON {
	POINT_PROCESS tmNmda
	RANGE tau1, tau2, e, i, q, mg
	RANGE tau, tauR, tauF, U, u0
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
        : q = 2, now included in parameters
	tau1= 2.815 (ms) <1e-9, 1e9> : ORIG : 5.63 ms
	tau2 = 160 (ms) <1e-9, 1e9> : ORIG : 320 ms
	e = 0	(mV)
	tau = 3 (ms) <1e-9, 1e9>
	tauR = 100 (ms) <1e-9, 1e9>
	tauF = 800 (ms) <0, 1e9>
	U = 0.3 (1) <0, 1>
	u0 = 0 (1) <0, 1>
        mg = 1 (mM)

}

ASSIGNED {
	v (mV)
	i (nA)
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
	LOCAL mggate
	SOLVE state METHOD cnexp
	mggate = 1 / (1 + exp(0.062 (/mV) * -(v)) * (mg / 3.57 (mM)))
	g = B - A
	i = g*(v - e)*mggate
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
	z = z + (y*(exp(-(t-tsyn)/tau) - exp(-(t-tsyn)/tauR)) / ((tau/tauR)-1) )
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
