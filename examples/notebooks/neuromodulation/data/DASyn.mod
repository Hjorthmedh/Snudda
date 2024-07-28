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


: 33000 dopamine molecules per vesicle : Omiatek, D., Bressler, A.,
: Cans, AS. et al. The real catecholamine content of secretory vesicles
: in the CNS revealed by electrochemical cytometry. Sci Rep 3, 1447
: (2013). https://doi.org/10.1038/srep01447



NET_RECEIVE(weight) {
    open = open + weight
}
