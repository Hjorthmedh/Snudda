COMMENT

c1 - c2 - c3 - c4 - o
|    |    |    |    |
i1 - i2 - i3 - i4 - i5 - is

ENDCOMMENT




NEURON {
	SUFFIX kv4_ch
	USEION k READ ek WRITE ik
	RANGE g, ik, gbar
	GLOBAL alpha, beta
	GLOBAL ci, ic, oi, io, a, b, am, bm, vc, gamma, delta, vha, vhb
	GLOBAL i5is, isi5
	GLOBAL q10i, q10v
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	gbar = 1	(S/cm2)
	gamma = 200	(1/ms)
	delta = 4	(1/ms)
	a = 3
	b = 40
	ic = 500	(/ms)
	oi = 1e-9	(/ms)
	io = .01	(/ms)
	ci = .2		(/ms)
	am = 1		(1/ms)
	bm = 7		(1/ms)
	vc = 10		(mV)
	vha = -75	(mV)
	vhb = -30	(mV)
	i5is = .001	(/ms)
	isi5 = .001	(/ms)
	q10i = 3
	q10v = 3
	celsius		(degC)
}

ASSIGNED {
	v	(mV)
	ek	(mV)
	g	(S/cm2)
	ik	(mA/cm2)
	alpha	(/ms)
	beta	(/ms)
}

STATE {
	c1
	c2
	c3
	c4
	o
	i1
	i2
	i3
	i4
	i5
	is
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*o
	ik = g*(v-ek)
}

INITIAL {
	SOLVE kin STEADYSTATE sparse
}

KINETIC kin{LOCAL q10
	q10 = q10i^((celsius - 22 (degC))/10 (degC))
	rates(v)
	~ c1 <-> c2 (3*alpha,beta)
	~ c2 <-> c3 (2*alpha,2*beta)
	~ c3 <-> c4 (alpha,3*beta)
	~ c4 <-> o  (q10*gamma,q10*delta)

	~ i1 <-> i2 (3*alpha*a,beta/b)
	~ i2 <-> i3 (2*alpha*a,2*beta/b)
	~ i3 <-> i4 (alpha*a,3*beta/b)
	~ i4 <-> i5 (q10*gamma,q10*delta)
	~ i5 <-> is (q10*i5is,q10*isi5)

	~ i1 <-> c1 (q10*ic,q10*ci)
	~ i2 <-> c2 (q10*ic/b,q10*ci*a)
	~ i3 <-> c3 (q10*ic/b^2,q10*ci*a^2)
	~ i4 <-> c4 (q10*ic/b^3,q10*ci*a^3)
	~ i5 <-> o  (q10*io,q10*oi)

	CONSERVE c1+c2+c3+c4+i1+i2+i3+i4+i5+o+is=1
}

PROCEDURE rates(v(millivolt)) {LOCAL q10
	q10 = q10v^((celsius - 22 (degC))/10 (degC))
	alpha = q10*am*exp((v-vha)/vc)
	beta = q10*bm*exp((v-vhb)/-vc)
}