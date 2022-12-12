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

        self.morphology_data = dict()

    def add_morphology(self, swc_file, name="neuron", position=None, rotation=None, parent=None):

        self.morphology_data[name] = MorphologyData(swc_file=swc_file, parent=parent)
        self.morphology_data[name].place(position=position, rotation=rotation)

    def section_iterator(self, section_type=None):
        for subtree in self.morphology_data.values():
            for section in subtree.section_iterator(section_type=section_type):
                yield section
