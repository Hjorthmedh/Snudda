NEURON {
    POINT_PROCESS DASyn
    RANGE quanta, tau, open
}
UNITS {
    (mM) = (milli / liter)
}

PARAMETER {
    quanta = 1e-4 (mM/ms)
    tau = 10 (ms)
}

INITIAL {
    open = 0
}

STATE { 
    open (1)
}

BREAKPOINT {SOLVE state METHOD cnexp}

DERIVATIVE state { 
    open' = -open/tau
}

NET_RECEIVE(weight) {
    open = open + weight
}
