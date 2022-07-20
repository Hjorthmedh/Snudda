import json
from neuron import hoc
import numpy as np
from snudda.neurons.neuron_model_extended import NeuronModel
import types
import nrn


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif type(obj) in [bool, np.bool_]:
            return int(obj)
        elif type(obj) in [bytes, np.bytes_]:
            return obj.decode()
        elif isinstance(obj, hoc.HocObject):
            if "to_python" in dir(obj):
                return obj.to_python()
            else:
                return str(obj)
        elif isinstance(obj, NeuronModel):
            return str(obj)
        elif isinstance(obj, types.BuiltinMethodType):
            return str(obj)
        elif isinstance(obj, nrn.Segment):
            return str(obj)
        else:
            # return super(NumpyEncoder, self).default(obj)
            return json.JSONEncoder.default(self, obj)
