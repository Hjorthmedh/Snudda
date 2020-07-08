TITLE I-h channel from Magee 1998 for distal dendrites

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
    SUFFIX hd_lts
    NONSPECIFIC_CURRENT i
    RANGE ghdbar, vhalfl
    GLOBAL linf,taul
    RANGE maxMod
    POINTER damod
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
    v 		        (mV)
    ehd  = -30		(mV)        
    celsius 	    (degC)
    ghdbar=.0001 	(mho/cm2)
    vhalfl=-81   	(mV)
    kl=-8
    vhalft=-75   	(mV)
    a0t=0.011      	(/ms)
    zetat=2.2    	(1)
    gmt=.4   	    (1)
    q10=4.5
    qtl=1
    damod = 0
    maxMod = 1
}



STATE {
    l
}

ASSIGNED {
	i (mA/cm2)
    linf      
    taul
    ghd
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = ghdbar*l*modulation()
	i = ghd*(v-ehd)

}


FUNCTION alpt(v(mV)) {
    alpt = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
    bett = exp(0.0378*zetat*gmt*(v-vhalft)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
    rate(v)
    l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) { :callable from hoc
    LOCAL a,qt
    qt=q10^((celsius-33)/10)
    a = alpt(v)
    linf = 1/(1 + exp(-(v-vhalfl)/kl))
    :linf = 1/(1+ alpl(v))
    taul = bett(v)/(qtl*qt*a0t*(1+a))
}



FUNCTION modulation() {
    : returns modulation factor
    
    modulation = 1 + damod*(maxMod-1)
}


