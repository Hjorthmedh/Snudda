import glob
import json
import os

from snudda.neurons import NeuronMorphology
from snudda.utils.snudda_path import snudda_parse_path


class NeuronPrototype:

    """ Helper class, returns a neuron prototype based on parameter_id, morph_id and modulation_id """

    def __init__(self, neuron_path, neuron_name,
                 morphology_path=None,
                 parameter_path=None,
                 mechanism_path=None,
                 modulation_path=None,
                 virtual_neuron=False,
                 load_morphology=True,
                 axon_stump_id_flag=False):

        if neuron_path:
            self.neuron_path = snudda_parse_path(neuron_path)
        elif parameter_path:
            self.neuron_path = os.path.dirname(snudda_parse_path(parameter_path))
        else:
            self.neuron_path = None

        if morphology_path:
            self.morphology_path = snudda_parse_path(morphology_path)
        elif self.neuron_path:
            self.morphology_path = snudda_parse_path(os.path.join(self.neuron_path, "morphology"))
        else:
            self.morphology_path = None

        if mechanism_path:
            self.mechanism_path = snudda_parse_path(mechanism_path)
        elif self.neuron_path:
            self.mechanism_path = snudda_parse_path(os.path.join(self.neuron_path, "mechanisms.json"))
        else:
            self.mechanism_path = None

        if parameter_path:
            self.parameter_path = snudda_parse_path(parameter_path)
        elif self.neuron_path:
            self.parameter_path = snudda_parse_path(os.path.join(self.neuron_path, "parameters.json"))
        else:
            self.parameter_path = None

        if modulation_path:
            self.modulation_path = snudda_parse_path(modulation_path)
        elif self.neuron_path:
            self.modulation_path = snudda_parse_path(os.path.join(self.neuron_path, "modulation.json"))

            if not os.path.exists(self.modulation_path):
                self.modulation_path = None
        else:
            self.modulation_path = None

        self.neuron_name = neuron_name
        self.parameter_info = None
        self.modulation_info = None
        self.virtual_neuron = virtual_neuron
        self.axon_stump_id_flag = axon_stump_id_flag
        self.load_morphology = load_morphology

        self.morphology_cache = dict()
        self.morphology_lookup = dict()

        self.load_info()

    def load_info(self):
        """ Reads information about the neuron prototype. """

        self.morphology_cache = dict()
        self.morphology_lookup = dict()

        par_path = self.parameter_path

        if par_path is None or not os.path.exists(par_path):
            print(f"Missing parameter.json : {par_path}")
            self.parameter_info = None
            return

        with open(par_path, "r") as f:
            self.parameter_info = json.load(f)

        # We expect a list of parameter sets, but if the user just provided one, convert it to a list
        if type(self.parameter_info[0]) == dict:
            self.parameter_info = [self.parameter_info]

        mod_path = self.modulation_path

        if mod_path is not None and os.path.exists(mod_path):
            with open(mod_path, "r") as f:
                self.modulation_info = json.load(f)
        else:
            self.modulation_info = None

    def get_morphology(self, parameter_id, morphology_id):

        """
        Returns morphology for a given parameter_id, morphology_id combination.
        (Each parameter set has a set of morphologies that it is valid for)
        """

        if self.parameter_info:
            par_set = self.parameter_info[parameter_id % len(self.parameter_info)]
        else:
            par_set = None

        if par_set is not None and len(par_set) > 0 and "morphology" in par_set[0]:
            morph_list = par_set[0]["morphology"]
            morph_path = os.path.join(self.morphology_path, morph_list[morphology_id % len(morph_list)])
        elif os.path.isfile(self.morphology_path):
            morph_path = self.morphology_path
        else:
            # No morphology
            file_list = glob.glob(os.path.join(self.neuron_path, "*swc"))

            if len(file_list) > 0:
                morph_path = file_list[0]

                assert len(file_list) == 1, f"More than one morphology available in {self.neuron_path}"

            else:
                print(f"No morphology in {self.neuron_path}")
                morph_path = None

        if morph_path is None:
            print("morph_path is None. Is SNUDDA_DATA set correctly?")

        return morph_path

    def get_parameters(self, parameter_id):
        """
        Returns neuron parameter set

        Args:
            parameter_id (int) : parameter ID, which of the parameter sets to use

        Returns:
            One parameter set as a list of dictionaries
        """

        if self.parameter_info:
            par_set = self.parameter_info[parameter_id % len(self.parameter_info)]
        else:
            par_set = []

        return par_set

    def get_modulation_parameters(self, modulation_id):

        """
        Returns modulation parameter set

        Args:
            modulation_id (int) : modulation ID, which of the modulation parameter sets to use

        Returns:
            One parameter set as a list of dictionaries
        """

        if self.modulation_info:
            modulation_set = self.modulation_info[modulation_id % len(self.modulation_info)]
        else:
            modulation_set = []

        return modulation_set

    def instantiate(self):
        """
        Instantiates all morphologies at once, instead of on demand.
        """
        for par_id in range(0, len(self.parameter_info)):
            print(f"Instantiates par_id = {par_id}")

            par_data = self.get_parameters(parameter_id=par_id)
            if par_data is not None and len(par_data) > 0 and "morphology" in par_data[0]:
                n_morph = len(par_data[0]["morphology"])
            else:
                n_morph = 1

            for morph_id in range(0, n_morph):
                morph_path = self.get_morphology(parameter_id=par_id, morphology_id=morph_id)
                morph_tag = os.path.basename(morph_path)

                if morph_tag not in self.morphology_cache:
                    print(f"morph_tag = {morph_tag}")
                    self.morphology_cache[morph_tag] = NeuronMorphology(swc_filename=morph_path,
                                                                        param_data=self.parameter_path,
                                                                        mech_filename=self.mechanism_path,
                                                                        name=self.neuron_name,
                                                                        hoc=None,
                                                                        load_morphology=self.load_morphology,
                                                                        virtual_neuron=self.virtual_neuron,
                                                                        axon_stump_id_flag=self.axon_stump_id_flag)

    def apply(self, function_name, arguments):
        """
        Applies function to all morphology prototypes
        """

        res = []

        for m in self.morphology_cache.values():
            function = getattr(m, function_name)
            res.append(function(*arguments))

        return res

    def all_have_axon(self):
        return all([len(m.axon) > 0 for m in self.morphology_cache.values()])

    def all_have_dend(self):
        return all([len(m.dend) > 0 for m in self.morphology_cache.values()])

    def clone(self, parameter_id, morphology_id, modulation_id,
              position=None, rotation=None, get_cache_original=False):
        """
        Creates a clone of the neuron prototype, with given position and rotation.

        Args:
            parameter_id (int) : parameter set to use
            morphology_id (int) : morphology set to use, a parameter set can have multiple morphology variations
            modulation_id (int) : neuro-modulation parameter set to use
            position (float,float,float) : position of neuron clone
            rotation : rotation (3x3 numpy matrix)
            get_cache_original (bool) : return the original rather than a clone

        """

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
                                                                    load_morphology=self.load_morphology,
                                                                    virtual_neuron=self.virtual_neuron,
                                                                    axon_stump_id_flag=self.axon_stump_id_flag)
        else:
            morph_tag = self.morphology_lookup[parameter_id, morphology_id]

        if get_cache_original:
            assert position is None and rotation is None and modulation_id is None, \
                    "If get_cache_original is passed position, rotation and modulation_id must be None"
            morph = self.morphology_cache[morph_tag]
        else:
            morph = self.morphology_cache[morph_tag].clone(position=position,
                                                           rotation=rotation,
                                                           parameter_id=parameter_id,
                                                           morphology_id=morphology_id,
                                                           modulation_id=modulation_id)
        return morph

