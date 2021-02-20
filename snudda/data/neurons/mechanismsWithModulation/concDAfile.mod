COMMENT
Two state kinetic scheme  which reproduces the concentration change of acetylcholine in response to spiking in cholinergic interneuron. The time constants, tau1 and tau2 were fitted to reproduce the same concentration as in Blackwell et al. 2019 with estimated input.

ENDCOMMENT

NEURON {
	POINT_PROCESS concDAfile
	RANGE tau1, tau2, e
	RANGE concentration, amplitude
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	
}

PARAMETER {
	tau1 = 400 (ms) <1e-9,1e9>
	tau2 = 300 (ms) <1e-9,1e9>
	amplitude = 1     
        concentration
				      
}

ASSIGNED {
     factor
}

STATE {
	A 
	B 
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
}

BREAKPOINT {
	concentration = amplitude*concentration		   
	
}
