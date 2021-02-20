TITLE KCNQ/M current from striatal medium spiny neurons, used for ch. interneuron.

:KCNQ_CH.MOD
:
: 11/3/2003

NEURON {
        SUFFIX kcnq_ch
        USEION k READ ek WRITE ik
        RANGE ik, ek, g, gbar
        GLOBAL a0, b0, ah, bh, ac, bc
}

UNITS {
        (mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
        gbar    = 1		(S/cm2)
        ek			(mV)
        a0      = .018		(/ms)	: parameters for alpha and beta
        b0      = .01		(/ms)
        ah      = -43.367	(mV)
        bh      = -43.367	(mV)
        ac      = 9.7054	(mV)
        bc      = -9.7054	(mV)
        q10v    = 3
        celsius			(degC)
}

ASSIGNED {
        v	(mV)
        g	(S/cm2)
        ik	(mA/cm2)
        alpha   (/ms)
        beta    (/ms)
}

STATE {
	c
	o
}

INITIAL {
	SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
        SOLVE kin METHOD sparse
        g = gbar*o
        ik = g*(v-ek) 
}

KINETIC kin {
        rates(v)
        ~ c <-> o       (alpha, beta)
        CONSERVE c + o = 1
}

PROCEDURE rates(v(mV)) {
        LOCAL qv
        qv = q10v^((celsius-22 (degC))/10 (degC))
        alpha = a0*qv / (1 + exp(-(v-ah)/ac))
        beta = b0*qv / (1 + exp(-(v-bh)/bc))
}
