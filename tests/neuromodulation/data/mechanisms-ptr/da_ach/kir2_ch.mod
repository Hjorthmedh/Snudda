:Kir2_ch.MOD
: Kir2, inwardly rectifying channel


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

NEURON {
	SUFFIX kir2_ch
	USEION k READ ek WRITE ik
	RANGE g, ninf, tn, ik, gbar
	GLOBAL C_tn, vh, vc
        RANGE modDA, maxModDA, levelDA
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
        modDA = 0
        maxModDA = 1
        levelDA = 0
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
	g = gbar*n*modulationDA()
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

FUNCTION modulationDA() {
    : returns modulation factor
    
    modulationDA = 1 + modDA*(maxModDA-1)*levelDA 
}
