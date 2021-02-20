:Kir2_ch.MOD
: Kir2, inwardly rectifying channel


COMMENT

neuromodulation is added as functions:
    
    modulation = 1 + damod*(maxMod-1)

where:
    
    damod  [0]: is a switch for turning modulation on or off {1/0}
    maxMod [1]: is the maximum modulation for this specific channel (read from the param file)
                    e.g. 10% increase would correspond to a factor of 1.1 (100% +10%) {0-inf}

[] == default values
{} == ranges
    
ENDCOMMENT

NEURON {
	SUFFIX kir2_ch
	USEION k READ ek WRITE ik
	RANGE g, ninf, tn, ik, gbar
	GLOBAL C_tn, vh, vc
        RANGE maxMod
        POINTER damod
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
    damod = 0
    maxMod = 1
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
	g = gbar*n*modulation()
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

FUNCTION modulation() {
    : returns modulation factor
    
    modulation = 1 + damod*(maxMod-1)
}
