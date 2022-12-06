import numpy as np

from snudda.neurons.morphology_data import MorphologyData, SectionMetaData
from snudda.utils.snudda_path import snudda_parse_path

class NeuronMorphologyExtended:

    def __init__(self, name, position, rotation, swc_filename, snudda_data,
                 parameter_key, morphology_key, modulation_key):

        self.name = name
        self.position = position
        self.rotation = rotation
        self.swc_filename = swc_filename
        self.snudda_data = snudda_data

        self.parameter_key = parameter_key
        self.morphology_key = morphology_key
        self.modulation_key = modulation_key

        pass

    def load_swc(self, swc_file, append_flag=False):

        with open(swc_file, "r") as f:
