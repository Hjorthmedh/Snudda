import numpy as np
import numexpr


def alpha_sub_function(t_step, tau, tstart, gmax):

    return gmax * ((t_step - tstart) / tau) * np.exp(1 - ((t_step - tstart) / tau))


def alpha(parameter=None):

    time_step_array = parameter['time_step_array']
    tstart = parameter['tstart']
    gmax = parameter['gmax']
    tau = parameter['tau']

    magnitude = np.zeros_like(time_step_array)
    index = np.where(time_step_array > tstart)
    start_index = index[0][0]

    magnitude[start_index:] = alpha_sub_function(np.take(time_step_array, index), tau, tstart, gmax)

    return magnitude


def step(parameter=None):

    time_step_array = parameter['time_step_array']
    tstart = parameter['tstart']
    step_stop = parameter['duration'] + parameter['tstart']
    gmax = parameter['gmax']

    magnitude = np.zeros_like(time_step_array)

    start_index = np.where(np.logical_and(time_step_array > tstart, time_step_array < step_stop))[0][0]

    magnitude[start_index:] = gmax
    
    return magnitude


def bath_application(parameter=None):

    time_step_array = parameter['time_step_array']
    gmax = parameter['gmax']

    magnitude = np.ones_like(time_step_array) * gmax

    return magnitude


def alpha_background(parameter=None):

    time_step_array = parameter['time_step_array']
    tstart = parameter['tstart']
    if 'gmax_decrease' in parameter.keys():
        gmax_shift = parameter['gmax_decrease'] * (-1)
    elif 'gmax_increase' in parameter.keys():
        gmax_shift = parameter['gmax_increase']
    else:
        raise ValueError('Include a gmax increase or decrease')

    tau = parameter['tau']
    tonic = parameter['tonic']

    magnitude = np.ones_like(time_step_array) * tonic

    index = np.where(time_step_array > tstart)[0]
    start_index = np.where(time_step_array > tstart)[0][0]

    magnitude[start_index:] = tonic + alpha_sub_function(np.take(time_step_array, index), tau, tstart, gmax_shift)

    if min(magnitude) < 0 or max(magnitude) > 1:
        raise ValueError(' Modulation is outside the range (0,1). Modify parameters')

    return magnitude


def time_series(parameter=None):

    magnitude = numexpr.evalute(parameter['array'])
    
    return magnitude
