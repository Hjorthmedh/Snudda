COMMENT
Two state kinetic scheme synapse described by rise time tau1,
and decay time constant tau2. The normalized peak conductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 is very small compared to tau1, this is an alphasynapse with time constant tau2.
If tau1/tau2 is very small, this is single exponential decay with time constant tau2.

The factor is evaluated in the initial block
such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
    POINT_PROCESS nACh
    RANGE tau1, tau2, e, i, q
    RANGE tau, tauR, tauF, U, u0
    RANGE failRate
    NONSPECIFIC_CURRENT i

    RANDOM release_probability
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    : q = 2,
    tau1= 0.15 (ms)
    tau2 = 40 (ms)
    e = 0 (mV)
    tau = 100 (ms)
    tauR = 800 (ms)
    tauF = 0 (ms)
    U = 0.95 (1) <0, 1>
    u0 = 0 (1) <0, 1>
    failRate = 0
}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
    factor
    x
}

STATE {
    A (uS)
    B (uS)
}

INITIAL {
    LOCAL tp
    A = 0
    B = 0
    tp = (tau1*tau2)/(tau2-tau1) * log(tau2/tau1)
    factor = -exp(-tp/tau1) + exp(-tp/tau2)
    factor = 1/factor
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    g = B - A
    i = g*(v - e)
}

DERIVATIVE state {
    A' = -A/tau1
    B' = -B/tau2
}

NET_RECEIVE(weight (uS), y, z, u, tsyn (ms)) {
    LOCAL result
    INITIAL {
        y = 0
        z = 0
        u = u0
        tsyn = t
    }
    if ( weight <= 0 ) {
VERBATIM
        return;
ENDVERBATIM
    }
    if( urand() > failRate ){
      z = z*exp(-(t-tsyn)/tauR)
      z = z + (y*(exp(-(t-tsyn)/tau) - exp(-(t-tsyn)/tauR)) / (tau/tauR - 1) )
      y = y*exp(-(t-tsyn)/tau)
      x = 1-y-z
      if (tauF > 0) {
          u = u*exp(-(t-tsyn)/tauF)
          u = u + U*(1-u)
      } else {
          u = U
      }
    A = A + weight*factor*x*u / U
    B = B + weight*factor*x*u / U
    y = y + x*u
    tsyn = t
    }
}

FUNCTION urand() {
    urand = random_uniform(release_probability)
}

COMMENT
(2025-10-08) NEURON 9.0+ compatibility. Replaced scop_random with the
new RANDOM keyword.
See: https://nrn.readthedocs.io/en/latest/nmodl/language/nmodl_neuron_extension.html#random

ENDCOMMENT
