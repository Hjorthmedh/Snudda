: KV2_CH.MOD
: KV2 for Cholinergic Interneuron

NEURON {
	SUFFIX kv2_ch
	USEION k READ ek WRITE ik
	RANGE g, ik, an, bn, gbar
}

UNITS {
	(mV) = (millivolt)
	(S) = (siemens)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 1	(S/cm2)
}

ASSIGNED {
	v	(mV)
	ek	(mV)
	g	(S/cm2)
	ik	(mA/cm2)

	an	(1/ms)
	bn	(1/ms)
	kf1	(1/ms)
	kb1	(1/ms)
	kf2	(1/ms)
	kb2	(1/ms)
	kf3	(1/ms)
	kb3	(1/ms)
	kf4	(1/ms)
	kb4	(1/ms)
}

STATE {
	c1
	c2
	c3
	c4
	o
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*o
	ik = g*(v-ek)
}

INITIAL {
	SOLVE kin STEADYSTATE sparse
}



KINETIC kin{
	rates(v)
	~ c4 <-> c3     (kf1,kb1)
	~ c3 <-> c2     (kf2,kb2)
	~ c2 <-> c1     (kf3,kb3)
	~ c1 <-> o      (kf4,kb4)
	CONSERVE c4+c3+c2+c1+o=1
}


PROCEDURE rates(v(millivolt)) {
		    : an = (51.743 (1/ms) - (0.612 (1/ms-mV))*v)/(exp((-84.54 (mV)+v)/-11.84 (mV)) - 1)
        an = ancalc(v)								      
	bn = (0.0051 (1/ms) )/exp(v/22.02 (mV))

	kf1 = 4*an
	kb1 = bn
	kf2 = 3*an
	kb2 = 2*bn
	kf3 = 2*an
	kb3 = 3*bn
	kf4 = an
	kb4 = 4*bn
}

FUNCTION ancalc(vx(millivolt)) {
     UNITSOFF
     if(80 < vx && vx < 90 ) {
          ancalc = 0.004186*vx^2  -0.40239*vx + 11.34  
       }
     else {
	  ancalc = (51.743 - 0.612*vx)/(exp((-84.54 +vx)/-11.84 ) - 1)     
     }
     UNITSON
}     
