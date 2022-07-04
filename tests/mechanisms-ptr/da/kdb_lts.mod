TITLE  K-D channel
: M.Migliore jun 2006

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
        ek (mV)		: must be explicitely def. in hoc
	celsius		(degC)
	gkdbar=.0 (mho/cm2)
        vhalfn=-33   (mV)
        a0n=0.005      (/ms)
        zetan=3    (1)
        gmn=0.7  (1)
	nmax=2  (1)
	q10=1
	sh = 0
}


NEURON {
	SUFFIX kdb_lts
	USEION k READ ek WRITE ik
        RANGE gkd,gkdbar, sh
	GLOBAL ninf,taun
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        gkd
        taun
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkd = gkdbar*n
	ik = gkd*(v-ek)

}

INITIAL {
	rates(v)
	n=ninf
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn-sh)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn-sh)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(qt*a0n*(1+a))
	if (taun<nmax) {taun=nmax/qt}
}














