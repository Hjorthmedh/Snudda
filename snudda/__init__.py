from .core import Snudda

__version__ = "1.0.11"

from .init import SnuddaInit
from .place import SnuddaPlace
from .detect import SnuddaDetect
from .detect import SnuddaProject
from .detect import SnuddaPrune
from .input import SnuddaInput

# The user has to explicity import SnuddaSimulate

# Network creation should work without having NEURON installed, you only need
# it to run the simulation

# from .simulate import SnuddaSimulate
