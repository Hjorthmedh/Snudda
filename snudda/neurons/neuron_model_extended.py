""" Based on BluePyOpt exampel code by Werner van Geit, 
modified by Johannes Hjorth """

import os
import json
from collections import OrderedDict

import numpy as np

import bluepyopt.ephys as ephys

from snudda.neurons.neuron_prototype import NeuronPrototype


class NeuronModel(ephys.models.CellModel):

    """ Extended NeuronModel for simulation. """

    def __init__(self,
                 cell_name="Unknown",
                 morph_path=None,
                 mech_file=None,
                 param_file=None,
                 modulation_file=None,
                 parameter_id=None,
                 morphology_id=None,
                 modulation_id=None,
                 parameter_key=None,
                 morphology_key=None,
                 modulation_key=None):

        """
        Constructor

        Args:
            cell_name: Name of neuron, e.g. "dSPN"
            morph_path: Path to morphology directory, there may be multiple morphologies available
            mech_file: Path to mechanism file
            param_file: Path to parameter file
            modulation_file: Path to neuromodulation parameter file
            parameter_id: ID of parameter set
            morphology_id: ID of morphology set
            modulation_id: ID of neuromodulation parameter set
            parameter_key (str): parameter key for lookup in parameter.json
            morphology_key (str): morphology key, together with parameter_key lookup in meta.json
            modulation_key (str): modulation key, lookup in modulation.json

        """

        # OBS, modulation_id is not used. morphology_id, parameter_id and modulation_id will be deprecated
        # at some point. Key usage is preferred (less risk of mixing models)

        self.name = cell_name
        self.type = cell_name.split('_')[0]
        self.parameters = []

        self.script_dir = os.path.dirname(__file__)
        self.config_dir = os.path.join(self.script_dir, 'config')

        morph_file = None
        if parameter_key is not None and morphology_key is not None and param_file is not None:
            meta_file = os.path.join(os.path.dirname(param_file), "meta.json")
            if os.path.isfile(meta_file):
                with open(meta_file, "r") as f:
                    meta_data = json.load(f, object_pairs_hook=OrderedDict)
                    morph_file = os.path.join(morph_path, meta_data[parameter_key][morphology_key]["morphology"])

        if not morph_file:
            # We now allow multiple variations of morphologies for a given neuron name, so here NeuronPrototype
            # is used to acquire the actual morphology file we will use for this particular neuron
            # based on parameter_id and morphology_id.
            neuron_prototype = NeuronPrototype(neuron_name=cell_name,
                                               neuron_path=None,
                                               morphology_path=morph_path,
                                               parameter_path=param_file,
                                               mechanism_path=mech_file,
                                               modulation_path=modulation_file)

            morph_file = neuron_prototype.get_morphology(parameter_id=parameter_id, morphology_id=morphology_id,
                                                         parameter_key=parameter_key, morphology_key=morphology_key)

        assert morph_file, (f"Neuron {cell_name} with morph_path = {morph_path} ({morphology_id}, "
                            f"parameter_path = {param_file} ({parameter_id}) "
                            f"has morph_file = {morph_file} (Should not be None)")

        morph = self.define_morphology(replace_axon=True, morph_file=morph_file)
        mechs = self.define_mechanisms(mechanism_config=mech_file)
        params = self.define_parameters(param_file, parameter_id, parameter_key)
        
        if modulation_key and modulation_file:
            
            mod_params = self.define_parameters(parameter_config=modulation_file, parameter_key=modulation_key)
            params = params + mod_params

        super(NeuronModel, self).__init__(name=cell_name, morph=morph,
                                          mechs=mechs, params=params)
        self.syn_list = []
        self.section_lookup = None

    ##############################################################################

    # Helper function

    def define_mechanisms(self, mechanism_config=None):
        """Define mechanisms based on mechanism_config """

        assert (mechanism_config is not None)
        # print("Using mechanmism config: " + mechanism_config)

        with open(mechanism_config, "r") as f:
            mech_definitions = json.load(f, object_pairs_hook=OrderedDict)

        if "modpath" in mech_definitions:
            mod_path = os.path.join(self.script_dir, mech_definitions["modpath"])
            print("mod_path set to " + mod_path + " (not yet implemented)")
        else:
            mod_path = None

        # import pdb
        # pdb.set_trace()

        mechanisms = []
        for section_list in mech_definitions:
            channels = mech_definitions[section_list]

            # This allows us to specify a modpath in the file
            if section_list == "modpath":
                continue

            seclist_loc = \
                ephys.locations.NrnSeclistLocation(section_list,
                                                   seclist_name=section_list)
            for channel in channels:
                mechanisms.append(ephys.mechanisms.NrnMODMechanism(
                    name='%s.%s' % (channel, section_list),
                    mod_path=mod_path,
                    suffix=channel,
                    locations=[seclist_loc],
                    preloaded=True))

        return mechanisms

    ############################################################################

    # Helper function

    def define_parameters(self, parameter_config=None, parameter_id=None, parameter_key=None):
        """
        Define parameters based on parameter_config and parameter_id.
        If there are n parameter sets, and parameter_id is k, then the parameter set is n % k."""

        assert (parameter_config is not None)

        # print("Using parameter config: " + parameter_config)

        with open(parameter_config, "r") as f:
            param_configs = json.load(f, object_pairs_hook=OrderedDict)

        parameters = []

        if type(param_configs) == OrderedDict:
            # Multiple parameters, pick one

            if parameter_key is not None:
                assert parameter_key in param_configs, f"Missing parameter key {parameter_key}"
                p_config = param_configs[parameter_key]

                if parameter_id is not None:
                    par_key_list = list(param_configs.keys())

                    par_key = par_key_list[parameter_id % len(par_key_list)]
                    assert par_key == parameter_key, \
                        (f"parameter_key mismatch, got {parameter_key}, but parameter_id={parameter_id} "
                         f"implies param_key={par_key}")
            else:
                assert parameter_id is not None, "Multiple parameter sets require parameter_id set"

                par_key_list = list(param_configs.keys())
                par_key = par_key_list[parameter_id % len(par_key_list)]
                p_config = param_configs[par_key]

        elif type(param_configs) == list and type(param_configs[0]) == list:
            # This is old fallback code, for old version format of parameters.json, remove in the future.
            num_params = len(param_configs)
            p_config = param_configs[parameter_id % num_params]
        else:
            p_config = param_configs

        # Save this to be accessible in the future
        self.parameters += p_config

        for param_config in p_config:
            if 'value' in param_config:
                frozen = True
                value = param_config['value']
                bounds = None
            elif 'bounds' in param_config:
                frozen = False
                bounds = param_config['bounds']
                value = None
            else:
                raise Exception(
                    'Parameter config has to have bounds or value: %s'
                    % param_config)

            if param_config['type'] == 'global':
                parameters.append(
                    ephys.parameters.NrnGlobalParameter(
                        name=param_config['param_name'],
                        param_name=param_config['param_name'],
                        frozen=frozen,
                        bounds=bounds,
                        value=value))
            elif param_config['type'] in ['section', 'range']:
                if param_config['dist_type'] == 'uniform':
                    scaler = ephys.parameterscalers.NrnSegmentLinearScaler()
                elif param_config['dist_type'] in ['exp', 'distance']:
                    scaler = ephys.parameterscalers.NrnSegmentSomaDistanceScaler(
                        distribution=param_config['dist'])
                seclist_loc = ephys.locations.NrnSeclistLocation(
                    param_config['sectionlist'],
                    seclist_name=param_config['sectionlist'])

                name = '%s.%s' % (param_config['param_name'],
                                  param_config['sectionlist'])

                if param_config['type'] == 'section':
                    parameters.append(
                        ephys.parameters.NrnSectionParameter(
                            name=name,
                            param_name=param_config['param_name'],
                            value_scaler=scaler,
                            value=value,
                            frozen=frozen,
                            bounds=bounds,
                            locations=[seclist_loc]))
                elif param_config['type'] == 'range':
                    parameters.append(
                        ephys.parameters.NrnRangeParameter(
                            name=name,
                            param_name=param_config['param_name'],
                            value_scaler=scaler,
                            value=value,
                            frozen=frozen,
                            bounds=bounds,
                            locations=[seclist_loc]))
            else:
                raise Exception(
                    'Param config type has to be global, section or range: %s' %
                    param_config)

            # import pdb
            # pdb.set_trace()

        return parameters

    ##############################################################################

    # Helper function

    def define_morphology(self, replace_axon=True, morph_file=None):
        """
        Define morphology. Handles SWC and ASC.

        Args:
            replace_axon (bool): Is axon replaced with a stump for simulation, default True
            morph_file (str): Path to morphology
        """

        assert (morph_file is not None)

        # print("Using morphology: " + morph_file)

        return ephys.morphologies.NrnFileMorphology(morph_file, do_replace_axon=replace_axon)

    # OLD BUGFIX FOR segment pop
    # ,replace_axon_hoc=self.getReplacementAxon())

    ##############################################################################

    # Neuron_morphology defines sectionID, these must match what this returns
    # so that they point to the same compartment.
    #
    # Soma is 0
    # axons are negative values (currently all set to -1) in Neuron_morphology
    # dendrites are 1,2,3,4,5... ie one higher than what Neuron internally
    # uses to index the dendrites (due to us wanting to include soma)

    def map_id_to_compartment(self, section_id):

        """
        Map section_id to compartment.

        Neuron_morphology defines sectionID, these must match what this returns
        so that they point to the same compartment.

        Soma is 0
        axons are negative values (currently all set to -1) in Neuron_morphology
        dendrites are 1,2,3,4,5... ie one higher than what Neuron internally
        uses to index the dendrites (due to us wanting to include soma)

        """

        if self.section_lookup is None:

            self.section_lookup = dict([])

            # Soma is zero
            self.section_lookup[0] = self.icell.soma[0]

            # Dendrites are consequtive numbers starting from 1
            # Ie neurons dend(0) is in pos 1, dend(99) is in pos 100
            # This so we dont need to special treat soma (pos 0)

            for ic, c in enumerate(self.icell.dend):
                self.section_lookup[ic + 1] = c

            # Negative numbers for axon
            for ic, c in enumerate(self.icell.axon):
                self.section_lookup[-ic - 1] = c

        sec = [self.section_lookup[x] for x in section_id]

        return sec

    # OVERRIDE the create_empty_template from CellModel
    '''create an hoc template named template_name for an empty cell'''

    @staticmethod
    def create_empty_template(
            template_name,
            seclist_names=None,
            secarray_names=None):

        objref_str = 'objref this, CellRef, synlist'
        new_seclist_str = ''

        if seclist_names:
            for seclist_name in seclist_names:
                objref_str += ', %s' % seclist_name
                new_seclist_str += \
                    '       %s = new SectionList()\n' % seclist_name

        create_str = ''
        if secarray_names:
            create_str = 'create '
            create_str += ', '.join(
                '%s[1]' % secarray_name
                for secarray_name in secarray_names)
            create_str += '\n'

        template = '''\
    begintemplate %(template_name)s
      %(objref_str)s
      proc init() {\n%(newseclist_str)s
      forall delete_section()
      CellRef = this
      synlist = new List()
      }
      gid = 0
      proc destroy() {localobj nil
      CellRef = nil
      }
      %(create_str)s
    endtemplate %(template_name)s
         ''' % dict(template_name=template_name, objref_str=objref_str,
                    newseclist_str=new_seclist_str,
                    create_str=create_str)

        # print(">>>> begin template")
        # print(str(template))
        # print(">>>> end template")

        return template

    ############################################################################

    # OVERRIDE the create_empty_template from CellModel
    """Create an empty cell in Neuron"""

    @staticmethod
    def create_empty_cell(
            name,
            sim,
            seclist_names=None,
            secarray_names=None):

        # TODO minize hardcoded definition
        # E.g. sectionlist can be procedurally generated
        # hoc_template = ephys.models.CellModel.create_empty_template(name,
        #                              seclist_names,
        #                              secarray_names)

        hoc_template = NeuronModel.create_empty_template(name,
                                                         seclist_names,
                                                         secarray_names)

        sim.neuron.h(hoc_template)

        template_function = getattr(sim.neuron.h, name)

        return template_function()

        ############################################################################

    # OVERRIDE the create_empty_template from CellModel

    def instantiate(self, sim=None):
        """Instantiate model in simulator"""

        # TODO replace this with the real template name
        if not hasattr(sim.neuron.h, self.name):
            self.icell = self.create_empty_cell(
                self.name,
                sim=sim,
                seclist_names=self.seclist_names,
                secarray_names=self.secarray_names)
        else:
            self.icell = getattr(sim.neuron.h, self.name)()

        self.icell.gid = self.gid

        self.morphology.instantiate(sim=sim, icell=self.icell)

        for mechanism in self.mechanisms:
            mechanism.instantiate(sim=sim, icell=self.icell)
        for param in self.params.values():
            param.instantiate(sim=sim, icell=self.icell)

    ############################################################################

    def get_replacement_axon(self):

        assert False, "Old bugfix for segment stack pop, should not be needed anymore"

        new_axon_hoc = \
            '''
proc replace_axon(){ local nSec, D1, D2
  // preserve the number of original axonal sections
  nSec = sec_count(axonal)
  // Try to grab info from original axon
  if(nSec == 0) { //No axon section present
    D1 = D2 = 1
  } else if(nSec == 1) {
    access axon[0]
    D1 = D2 = diam
  } else {
    access axon[0]
    D1 = D2 = diam

    //access soma distance() //to calculate distance from soma
    soma distance() //to calculate distance from soma
    forsec axonal{
      D2 = diam
      //if section is longer than 60um then store diam and exit from loop
      if(distance(0.5) > 60){
        break
      }
    }
  }
  // get rid of the old axon
  forsec axonal{
    delete_section()
  }
  create axon[2]

  //access axon[0] {
  axon[0] {
    L = 30
    diam = D1
    nseg = 1 + 2*int(L/40)
    all.append()
    axonal.append()
  }
  // access axon[1] {
  axon[1] {
    L = 30
    diam = D2
    nseg = 1 + 2*int(L/40)
    all.append()
    axonal.append()
  }
  nSecAxonal = 2
  soma[0] connect axon[0](0), 1
  axon[0] connect axon[1](0), 1

  if(nSec > 0) {
    pop_section()
  }


}
        '''

        return new_axon_hoc
