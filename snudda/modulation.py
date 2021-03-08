import numpy as np

def alpha(parameter=None):

    ht = parameter['ht']
    tstart = parameter['tstart']
    gmax = parameter['gmax']
    tau = parameter['tau']

    mag = list()
    
    '''
    calc and returns a "magnitude" using an alpha function -> used for modulation
        transients

    ht      = simulation time (h.t)
    tstart  = time when triggering the function
    gmax    = maximal amplitude of curve (default 1; transient must lie between 0-1)
    tau     = time constant of alpha function
    '''

    for t_step in ht:

        if t_step >= tstart:
            t = (t_step - tstart) / tau
            e = np.exp(1 - t)
            mag.append(gmax * t * e)

        else:
            mag.append(0)
    return mag


def step(parameter=None):

    ht = parameter['ht']
    tstart = parameter['tstart']
    tstop = parameter['tstop']
    gmax = parameter['gmax']

    mag = list()
    
    for t_step in ht:
        if ht > tstart and ht < tstop:
            mag.append(gmax)
        else:
            mag.append(0)
            

    return mag


def bath_application(parameter=None):

    ht = parameter['ht']
    gmax = parameter['gmax']

    mag = list()

    for t_step in ht:

        mag.append(gmax)

    return mag


def alpha_background(parameter=None):

    ht = parameter['ht']
    tstart = parameter['tstart']
    gmax_decrease = parameter['gmax_decrease']
    tau = parameter['tau']
    tonic = parameter['tonic']

    mag = list()

    for t_step in ht:

        for start in tstart:

            mag_intermediate = 0

            if t_step >= tstart:
                t = (t_step - tstart) / tau
                e = np.exp(1 - t)
                mag_intermediate = mag_intermediate + gmax_decrease * t * e

        if mag_intermediate > 0.99:
            mag_intermediate = 1
        mag.append(tonic - mag_intermediate)

    return mag


def time_series(parameter=None):

    mag = eval(parameter['array'])
    
    return mag
