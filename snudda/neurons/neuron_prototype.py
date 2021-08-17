import glob
import json
import os

from snudda.neurons import NeuronMorphology


class NeuronPrototype:

    """ Helper class, returns a neuron prototype based on parameter_id, morph_id and modulation_id """

    def __init__(self, neuron_path, neuron_name,
                 morphology_path=None,
                 parameter_path=None,
                 mechanism_path=None,
                 modulation_path=None,
                 virtual_neuron=False):

        if neuron_path:
            self.neuron_path = neuron_path
        else:
            assert morphology_path is not None \
                and parameter_path is not None \
                and mechanism_path is not None, \
                ("If neuron_path is empty then morphology_path, parameter_path, " 
                 "mechanism_path, modulation_path must be set")

            self.neuron_path = None

        if morphology_path:
            self.morphology_path = morphology_path
        else:
            self.morphology_path = os.path.join(self.neuron_path, "morphology")

        if mechanism_path:
            self.mechanism_path = mechanism_path
        else:
            self.mechanism_path = os.path.join(self.neuron_path, "mechanisms.json")

        if parameter_path:
            self.parameter_path = parameter_path
        else:
            self.parameter_path = os.path.join(self.neuron_path, "parameters.json")

        if modulation_path:
            self.modulation_path = modulation_path
        else:
            self.modulation_path = os.path.join(self.neuron_path, "modulation.json")

        self.neuron_name = neuron_name
        self.parameter_info = None
        self.modulation_info = None
        self.virtual_neuron = virtual_neuron

        self.morphology_cache = dict()
        self.morphology_lookup = dict()

        self.load_info()

    def load_info(self):
        """ Reads information about the neuron prototype. """

        with open(self.parameter_path, "r") as f:
            self.parameter_info = json.load(f)

        # We expect a list of parameter sets, but if the user just provided one, convert it to a list
        if type(self.parameter_info[0]) == dict:
            self.parameter_info = [self.parameter_info]

        if os.path.exists(self.modulation_path):
            with open(self.modulation_path, "r") as f:
                self.modulation_info = json.load(f)
        else:
            self.modulation_info = None

        self.morphology_cache = dict()
        self.morphology_lookup = dict()

    def get_morphology(self, parameter_id, morphology_id):

        """
        Returns morphology for a given parameter_id, morphology_id combination.
        (Each parameter set has a set of morphologies that it is valid for)
        """

        par_set = self.parameter_info[parameter_id % len(self.parameter_info)]

        if "morphology" in par_set[0]:
            morph_list = par_set[0]["morphology"]
            morph_path = os.path.join(self.morphology_path, morph_list[morphology_id % len(morph_list)])
        elif os.path.isfile(self.morphology_path):
            morph_path = self.morphology_path
        else:
            # No morphology
            file_list = glob.glob(os.path.join(self.neuron_path, "*swc"))
            assert len(file_list) == 1, f"There should only be one SWC file in {self.neuron_path}"
            morph_path = file_list[0]

        return morph_path

    def get_parameters(self, parameter_id):
        """
        Returns neuron parameter set
        """
        par_set = self.parameter_info[parameter_id % len(self.parameter_info)]
        return par_set

    def get_modulation_parameters(self, modulation_id):
        if self.modulation_info:
            modulation_set = self.modulation_info[modulation_id % len(self.modulation_info)]
        else:
            modulation_set = []

        return modulation_set

    def clone(self, parameter_id, morphology_id, modulation_id,
              position, rotation, load_morphology=False):

        if (parameter_id, morphology_id) not in self.morphology_lookup:
            morph_path = self.get_morphology(parameter_id=parameter_id, morphology_id=morphology_id)
            morph_tag = os.path.basename(morph_path)
            self.morphology_lookup[parameter_id, morphology_id] = morph_tag

            if morph_tag not in self.morphology_cache:
                # TODO: hoc file will depend on both morphology_id and parameter_id, we ignore it for now
                self.morphology_cache[morph_tag] = NeuronMorphology(swc_filename=morph_path,
                                                                    param_data=self.parameter_path,
                                                                    mech_filename=self.mechanism_path,
                                                                    name=self.neuron_name,
                                                                    hoc=None,
                                                                    virtual_neuron=self.virtual_neuron)
        else:
            morph_tag = self.morphology_lookup[parameter_id, morphology_id]

        morph = self.morphology_cache[morph_tag].clone(position=position,
                                                       rotation=rotation,
                                                       parameter_id=parameter_id,
                                                       modulation_id=modulation_id,
                                                       load_morphology=load_morphology)

        return morph
