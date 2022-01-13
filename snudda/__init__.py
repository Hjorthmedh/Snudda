from .core import Snudda

__version__ = "1.2.15"

from .init import SnuddaInit
from .place import SnuddaPlace
from .detect import SnuddaDetect
from .detect import SnuddaProject
from .detect import SnuddaPrune
from .utils import SnuddaLoad

# The user has to explicity import snudda.input and snudda.simulate
# We dont want to have NEURON dependencies for the standard imports

import os
if os.environ.get('READTHEDOCS') == 'True':
    from .input import SnuddaInput
    from .simulate import SnuddaSimulate
