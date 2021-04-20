import json
import numpy as np


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
        else:
            # return super(NumpyEncoder, self).default(obj)
            return json.JSONEncoder.default(self, obj)

