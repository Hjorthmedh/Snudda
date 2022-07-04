TITLE I-h channel from Magee 1998 for distal dendrites

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
	SUFFIX hd_lts
	NONSPECIFIC_CURRENT i
        RANGE ghdbar, vhalfl
        GLOBAL linf,taul
        RANGE modDA, maxModDA, levelDA
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
		        modDA = 0
                        maxModDA = 1
                        levelDA = 0
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
	ghd = ghdbar*l*modulationDA()
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


FUNCTION modulationDA() {
    : returns modulation factor
    
    modulationDA = 1 + modDA*(maxModDA-1)*levelDA 
}


