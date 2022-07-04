COMMENT
Two state kinetic scheme synapse which comes from Expsyn modfile which was fitted with new time constants to describe the effect of nitric oxide coupled to a spike integrating mechanisms (Infire1)
ENDCOMMENT

NEURON {
	POINT_PROCESS NO
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	RANGE g, m
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1 = 2000 (ms) <1e-9,1e9>
	tau2 = 3000 (ms) <1e-9,1e9>
	e=0	(mV)
	tau = 1000 (ms)
	refrac = 10 (ms)
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	m
	t0(ms)
	refractory
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	if (tau1/tau2 < 1e-9) {
		tau1 = tau2*1e-9
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
				     factor = 1/factor
						     m = 0
	t0 = t
	refractory = 0 : 0-integrates input, 1-refractory
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}


NET_RECEIVE(weight (uS)) {
     if (refractory == 0) { : inputs integrated only when excitable
		m = m*exp(-(t - t0)/tau)
		t0 = t
		m = m + 0.075
		if (m > 1) {
			refractory = 1
			m = 2
		
			A = A + weight*factor
			B = B + weight*factor
			
		}
	}else if (m==2) { : ready to integrate again
		t0 = t
		refractory = 0
		m = 0
	}

}
