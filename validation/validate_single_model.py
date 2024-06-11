import os
import json
from sciunit import Model, Capability

from snudda.synaptic_fitting.run_little_synapse_run import RunLittleSynapseRun

try:
    from clb_nb_utils import oauth
    have_collab_token_handler = True
except ImportError:
    have_collab_token_handler = False


class SnuddaSingleModel(Model, Capability):
    """Snudda single model"""

    def __init__(self,
                 stim_times,
                 holding_voltage=-70e-3,
                 synapse_type='glut',
                 params={},
                 time=2.0,
                 tau=None,
                 tau_r=None,
                 tau_f=None,
                 u=None,
                 cond=None):

        self.rlsr = RunLittleSynapseRun(stim_times=stim_times, holding_voltage=holding_voltage,
                                        synapse_type=synapse_type, params=params, time=time)


class

        pass
