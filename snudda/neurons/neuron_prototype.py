import glob
import json
import os
from collections import OrderedDict

from snudda.neurons import NeuronMorphology
from snudda.utils.snudda_path import snudda_parse_path


class NeuronPrototype:

    """ Helper class, returns a neuron prototype based on parameter_id, morph_id and modulation_id """

    def __init__(self,
                 neuron_path,
                 neuron_name,
                 morphology_path=None,
                 parameter_path=None,
                 mechanism_path=None,
                 modulation_path=None,
                 meta_path=None,
                 virtual_neuron=False,
                 load_morphology=True,
                 axon_stump_id_flag=False,
                 verbose=False):

        self.verbose = verbose

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

        if meta_path:
            self.meta_path = snudda_parse_path(meta_path)
        elif self.neuron_path:
            self.meta_path = snudda_parse_path(os.path.join(self.neuron_path, "meta.json"))
        else:
            self.meta_path = None

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
        self.meta_info = None
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

        if self.meta_path and os.path.exists(self.meta_path):
            with open(self.meta_path, "r") as fm:
                self.meta_info = json.load(fm, object_pairs_hook=OrderedDict)

        with open(par_path, "r") as f:
            self.parameter_info = json.load(f, object_pairs_hook=OrderedDict)

        # We now expect a dictionary of parameter sets. If it is a list, we convert it to a dictionary
        if type(self.parameter_info) == list:
            self.parameter_info = {"default": self.parameter_info}

        mod_path = self.modulation_path

        if mod_path is not None and os.path.exists(mod_path):
            with open(mod_path, "r") as f:
                self.modulation_info = json.load(f, object_pairs_hook=OrderedDict)
        else:
            self.modulation_info = None

    def get_num_morphologies(self, parameter_id):

        if self.parameter_info:
            par_key_list = list(self.parameter_info.keys())
            par_key = par_key_list[parameter_id % len(par_key_list)]
            par_set = self.parameter_info[par_key]
        else:
            par_key = None

        if self.meta_info and par_key:
            morph_key_list = self.meta_info[par_key].keys()
            return len(morph_key_list)
        else:
            return 1

    def get_parameter_key(self, parameter_id):

        if self.parameter_info:
            par_key_list = list(self.parameter_info.keys())
            par_key = par_key_list[parameter_id % len(par_key_list)]
        else:
            par_key = None

        return par_key

    def get_morph_key(self, parameter_id, morphology_id):

        par_key = self.get_parameter_key(parameter_id=parameter_id)
        if self.meta_info:
            assert par_key in self.meta_info, f"Parameter key {par_key} missing in {self.meta_path}"
            morph_key_list = list(self.meta_info[par_key].keys())
            assert len(morph_key_list) > 0, f"No morphologies available for parameter key {par_key} in {self.meta_path}"
            morph_key = morph_key_list[morphology_id % len(morph_key_list)]
        else:
            morph_key = None

        return morph_key

    def get_input_parameters(self, parameter_id, morphology_id):

        par_key = self.get_parameter_key(parameter_id=parameter_id)
        morph_key = self.get_morph_key(parameter_id=parameter_id, morphology_id=morphology_id)

        if self.meta_info and "input" in self.meta_info[par_key][morph_key]:
            input_info = self.meta_info[par_key][morph_key]["input"]
        else:
            input_info = dict()

        return input_info

    def get_morphology(self, parameter_id, morphology_id):

        """
        Returns morphology for a given parameter_id, morphology_id combination.
        (Each parameter set has a set of morphologies that it is valid for)
        """

        par_key = self.get_parameter_key(parameter_id=parameter_id)

        if self.meta_info and par_key:
            morph_key = self.get_morph_key(parameter_id=parameter_id, morphology_id=morphology_id)
            morph_path = os.path.join(self.morphology_path, self.meta_info[par_key][morph_key]["morphology"])
            assert os.path.isfile(morph_path), f"Morphology file {morph_path} is missing (listed in {self.meta_path})"

        elif self.morphology_path and os.path.isfile(self.morphology_path):
            # Fallback if morphology file is specified in the path
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
            par_key_list = list(self.parameter_info.keys())
            par_key = par_key_list[parameter_id % len(par_key_list)]
            par_set = self.parameter_info[par_key]
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
            if self.verbose:
                print(f"Instantiates par_id = {par_id}")

            par_data = self.get_parameters(parameter_id=par_id)
            n_morph = self.get_num_morphologies(parameter_id=par_id)

            for morph_id in range(0, n_morph):
                morph_path = self.get_morphology(parameter_id=par_id, morphology_id=morph_id)
                morph_tag = os.path.basename(morph_path)

                if morph_tag not in self.morphology_cache:
                    if self.verbose:
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

