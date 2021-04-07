
'''

Different functions for modulating the cells

The functions return array which contain the modulation level per time step of the simulation


'''


import numpy as np

def alpha(parameter=None):

    time_step_array = parameter['time_step_array']
    tstart = parameter['tstart']
    gmax = parameter['gmax']
    tau = parameter['tau']

    mag = np.zeros_like(time_step_array)

    mag = list()

    for t_step in time_step_array:

        if t_step >= tstart:
            t = (t_step - tstart) / tau
            e = np.exp(1 - t)
            mag.append(gmax * t * e)

        else:
            mag.append(0)
    return mag


def step(parameter=None):

    time_step_array = parameter['time_step_array']
    tstart = parameter['tstart']
    step_stop = parameter['duration'] + parameter['tstart']
    gmax = parameter['gmax']

    mag = list()
    
    for t_step in time_step_array:
        if t_step > tstart and t_step < step_stop:
            mag.append(gmax)
        else:
            mag.append(0)
            

    return mag


def bath_application(parameter=None):

    time_step_array = parameter['time_step_array']
    gmax = parameter['gmax']

    mag = list()

    for t_step in time_step_array:

        mag.append(gmax)

    return mag


def alpha_background(parameter=None):

    time_step_array = parameter['time_step_array']
    tstart = parameter['tstart']
    gmax_decrease = parameter['gmax_decrease']
    tau = parameter['tau']
    tonic = parameter['tonic']

    mag = list()

    for t_step in time_step_array:

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
