from .core import Snudda

<<<<<<< HEAD
__version__ = "1.0.9"
=======
__version__ = "1.0.11"

from .init import SnuddaInit
from .place import SnuddaPlace
from .detect import SnuddaDetect
from .detect import SnuddaProject
from .detect import SnuddaPrune

# The user has to explicity import snudda.input and snudda.simulate
# We dont want to have NEURON dependencies for the standard imports

# from .input import SnuddaInput
# from .simulate import SnuddaSimulate
>>>>>>> origin/dev
