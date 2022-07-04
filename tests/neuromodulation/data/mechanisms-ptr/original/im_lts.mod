TITLE Cortical M current
:
:   M-current, responsible for the adaptation of firing rate and the 
:   afterhyperpolarization (AHP) of cortical pyramidal cells
:
:   First-order model described by hodgkin-Hyxley like equations.
:   K+ current, activated by depolarization, noninactivating.
:
:   Model taken from Yamada, W.M., Koch, C. and Adams, P.R.  Multiple 
:   channels and calcium dynamics.  In: Methods in Neuronal Modeling, 
:   edited by C. Koch and I. Segev, MIT press, 1989, p 97-134.
:
:   See also: McCormick, D.A., Wang, Z. and Huguenard, J. Neurotransmitter 
:   control of neocortical neuronal activity and excitability. 
:   Cerebral Cortex 3: 387-398, 1993.
:
:   Written by Alain Destexhe, Laval University, 1995
:


COMMENT

Neuromodulation is added as functions:
    
    modulationDA = 1 + modDA*(maxModDA-1)*levelDA

where:
    
    modDA  [0]: is a switch for turning modulation on or off {1/0}
    maxModDA [1]: is the maximum modulation for this specific channel (read from the param file)
                    e.g. 10% increase would correspond to a factor of 1.1 (100% +10%) {0-inf}
    levelDA  [0]: is an additional parameter for scaling modulation. 
                Can be used simulate non static modulation by gradually changing the value from 0 to 1 {0-1}
									
	  Further neuromodulators can be added by for example:
          modulationDA = 1 + modDA*(maxModDA-1)
	  modulationACh = 1 + modACh*(maxModACh-1)
	  ....

	  etc. for other neuromodulators
	  
	   
								     
[] == default values
{} == ranges
    
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX im_lts
	USEION k READ ek WRITE ik
        RANGE gkbar, m_inf, tau_m, ik
	GLOBAL taumax
        RANGE modDA, maxModDA, levelDA
	RANGE modACh, maxModACh, levelACh
	

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {
	v		(mV)
	celsius = 36    (degC)
	ek		(mV)
	gkbar	= 1e-6	(mho/cm2)
	taumax	= 1000	(ms)		: peak value of tau
        modDA = 0
        maxModDA = 1
        levelDA = 0
        modACh = 0
        maxModACh = 1 
        levelACh = 0
}



STATE {
	m
}

ASSIGNED {
	ik	(mA/cm2)
	m_inf
	tau_m	(ms)
	tau_peak	(ms)
	tadj
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar * m * (v - ek)*modulationDA()*modulationACh()
}

DERIVATIVE states { 
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {
	evaluate_fct(v)
	m = 0
:
:  The Q10 value is assumed to be 2.3
:
        tadj = 2.3 ^ ((celsius-36)/10)
	tau_peak = taumax / tadj
}

PROCEDURE evaluate_fct(v(mV)) {

	m_inf = 1 / ( 1 + exptable(-(v+35)/10) )
	tau_m = tau_peak / ( 3.3 * exptable((v+35)/20) + exptable(-(v+35)/20) )
}
UNITSON


FUNCTION exptable(x) { 
	TABLE  FROM -25 TO 25 WITH 10000

	if ((x > -25) && (x < 25)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}

FUNCTION modulationDA() {
    : returns modulation factor
    
    modulationDA = 1 + modDA*(maxModDA-1)*levelDA 
}

FUNCTION modulationACh() {
    : returns modulation factor
    
    modulationACh = 1 + modACh*(maxModACh-1)*levelACh 
}
