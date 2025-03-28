""" Based on BluePyOpt exampel code by Werner van Geit, 
modified by Johannes Hjorth """

import json
import os
from collections import OrderedDict
import numpy as np

import bluepyopt.ephys as ephys

from snudda.neurons.neuron_prototype import NeuronPrototype

# Fix for axon diameters, need to refactor in the future, or integrate directly in bluepyopt
from snudda.neurons.NrnFileMorphology_axon_fix import NrnFileMorphology_axon_fix


class NeuronModel(ephys.models.CellModel):
    """ Extended NeuronModel for simulation.

        Note that the paths are here assumed to be absolute, ie SNUDDA_DATA is not replaced, that should be done
        before calling this class.

    """

    def __init__(self,
                 cell_name="Unknown",
                 morph_path=None,
                 mech_file=None,
                 param_file=None,
                 modulation_file=None,
                 reaction_diffusion_file=None,
                 parameter_id=None,
                 morphology_id=None,
                 modulation_id=None,
                 parameter_key=None,
                 morphology_key=None,
                 modulation_key=None,
                 use_rxd_neuromodulation=True,
                 replace_axon_length=60e-6,
                 replace_axon_nseg_frequency=40e-6,
                 replace_axon_diameter=None,  # If this is not None, it uses special code to get diameter also
                 replace_axon_myelin_length=None,
                 replace_axon_myelin_diameter=None,
                 position=None,
                 rotation=None,
                 volume_id=None):

        """
        Constructor

        Args:
            cell_name: Name of neuron, e.g. "dSPN"
            morph_path: Path to morphology directory, there may be multiple morphologies available
            mech_file: Path to mechanism file
            param_file: Path to parameter file
            modulation_file: Path to neuromodulation parameter file
            reaction_diffusion_file: Path to the RxD reaction diffusion file
            parameter_id: ID of parameter set
            morphology_id: ID of morphology set -- DEPRECATED
            modulation_id: ID of neuromodulation parameter set
            parameter_key (str): parameter key for lookup in parameter.json
            morphology_key (str): morphology key, together with parameter_key lookup in meta.json
            modulation_key (str): modulation key, lookup in modulation.json
            volume_id (str): name of volume (region) that neuron is located in

        """

        # OBS, modulation_id is not used. morphology_id, parameter_id and modulation_id will be deprecated
        # at some point. Key usage is preferred (less risk of mixing models)

        self.name = cell_name
        self.type = cell_name.split('_')[0]
        self.parameters = []

        self.script_dir = os.path.dirname(__file__)
        self.config_dir = os.path.join(self.script_dir, 'config')
        self.use_rxd_neuromodulation = use_rxd_neuromodulation

        self.is_relocated = False  # NEURON coordinates match Snudda simulation coordinates

        self.position = position
        self.rotation = rotation
        self.volume_id = volume_id

        if os.path.isfile(morph_path):
            # If morph_path is a swc file, use it directly
            morph_file = morph_path
        else:
            # Figure out what file in the morphology path we should use

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
                                                   modulation_path=modulation_file,
                                                   reaction_diffusion_path=reaction_diffusion_file)

                morph_file, _ = neuron_prototype.get_morphology(parameter_id=parameter_id,
                                                                morphology_id=morphology_id,
                                                                parameter_key=parameter_key,
                                                                morphology_key=morphology_key)

            assert morph_file, (f"Neuron {cell_name} with morph_path = {morph_path} ({morphology_id}, "
                                f"parameter_path = {param_file} ({parameter_id}) "
                                f"has morph_file = {morph_file} (Should not be None)")

        self.morph_file = morph_file

        morph = self.define_morphology(replace_axon=True, morph_file=morph_file,
                                       replace_axon_length=replace_axon_length,
                                       replace_axon_nseg_frequency=replace_axon_nseg_frequency,
                                       replace_axon_diameter=replace_axon_diameter,
                                       replace_axon_myelin_length=replace_axon_myelin_length,
                                       replace_axon_myelin_diameter=replace_axon_myelin_diameter)

        mechs = self.define_mechanisms(mechanism_config=mech_file)
        params = self.define_parameters(param_file, parameter_id, parameter_key)

        if modulation_file:
            if modulation_key:
                mod_params = self.define_parameters(parameter_config=modulation_file,
                                                    parameter_key=modulation_key)
                params = params + mod_params
            else:
                print(f"Warning! No modulation key specified, ignoring {modulation_file}")

        super(NeuronModel, self).__init__(name=cell_name, morph=morph,
                                          mechs=mechs, params=params)

        if reaction_diffusion_file and self.use_rxd_neuromodulation:
            # Only load the module if actually used, this avoids weird RxD stuff happening in NEURON when not needed
            from snudda.neurons.neuron_modulation import NeuronModulation

            self.modulation = NeuronModulation(neuron=self)
            self.modulation.config_file = reaction_diffusion_file
        else:
            self.modulation = None

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
            print(f"mod_path set to {mod_path} (not yet implemented)")
        else:
            mod_path = None

        mechanisms = []
        for section_list in mech_definitions:
            channels = mech_definitions[section_list]

            # This allows us to specify a modpath in the file
            if section_list == "modpath":
                continue

            seclist_loc = \
                ephys.locations.NrnSeclistLocation(section_list, seclist_name=section_list)
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

    def define_parameters(self, parameter_config, parameter_id=None, parameter_key=None):
        """
        Define parameters based on parameter_config and parameter_id.
        If there are n parameter sets, and parameter_id is k, then the parameter set is n % k."""

        assert parameter_config is not None

        # print("Using parameter config: " + parameter_config)

        with open(parameter_config, "r") as f:
            param_configs = json.load(f, object_pairs_hook=OrderedDict)

        parameters = []

        if isinstance(param_configs, (OrderedDict, dict)):
            # Multiple parameters, pick one

            if parameter_key is not None:
                assert parameter_key in param_configs, f"Missing parameter key {parameter_key} in {parameter_config}"
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
            # Parameter ID no longer exists, we default to parameter_id = 0
            # This is old fallback code, for old version format of parameters.json, remove in the future.
            print("Warning: Old format of parameter config, using parameter_id = 0.")
            parameter_id = 0
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
                raise Exception(f"Parameter config has to have bounds or value: {param_config}")

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
                    scaler = ephys.parameterscalers.NrnSegmentSomaDistanceScaler(distribution=param_config['dist'])
                else:
                    raise ValueError(f"Unknown dist_type = {param_config['dist_type']}, "
                                     f"expected 'uniform', 'exp' or 'distance'")
                # 2024-07-23: Updated format, so that "sectionlist" is allowed to be either a string (of one section type)
                #             or a list of strings with section types.

                section_list = param_config['sectionlist']
                if not isinstance(section_list, list):
                    section_list = [section_list]

                for seclist_item in section_list:

                    seclist_loc = ephys.locations.NrnSeclistLocation(seclist_item,
                                                                     seclist_name=seclist_item)

                    name = '%s.%s' % (param_config['param_name'],
                                      seclist_item)

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
                raise Exception(f"Param config type has to be global, section or range: {param_config}")

        return parameters

    ##############################################################################

    # Helper function

    def define_morphology(self, replace_axon=True, morph_file=None,
                          replace_axon_length=60e-6,
                          replace_axon_nseg_frequency=40e-6,
                          replace_axon_diameter=None,
                          replace_axon_myelin_length=None,
                          replace_axon_myelin_diameter=None):
        """
        Define morphology. Handles SWC and ASC.

        Args:
            replace_axon (bool): Is axon replaced with a stump for simulation, default True
            morph_file (str): Path to morphology
        """

        assert (morph_file is not None)

        if replace_axon_diameter is None and replace_axon_myelin_length is None:
            nrn_morph = ephys.morphologies.NrnFileMorphology(morph_file, do_replace_axon=replace_axon,
                                                             axon_stub_length=replace_axon_length*1e6,
                                                             axon_nseg_frequency=replace_axon_nseg_frequency*1e6)
        else:
            # Special treatment, if we have both length and diameter requirements
            print(f"Axon diameter given, assuming both axon diameter and length are lists or arrays")
            nrn_morph = NrnFileMorphology_axon_fix(morphology_path=morph_file,
                                                   axon_length=replace_axon_length,
                                                   axon_diameter=replace_axon_diameter,
                                                   axon_nseg_frequency=replace_axon_nseg_frequency,
                                                   replace_axon_myelin_length=replace_axon_myelin_length,
                                                   replace_axon_myelin_diameter=replace_axon_myelin_diameter)

        return nrn_morph

    ##############################################################################

    # Neuron_morphology defines sectionID, these must match what this returns
    # so that they point to the same compartment.
    #
    # Soma is -1
    # axons are negative values (-2, -3, ..) in Neuron_morphology
    # dendrites are 0,1,2,3,4,5... ie one higher than what Neuron internally
    # uses to index the dendrites (due to us wanting to include soma)

    def build_section_lookup(self):

        self.section_lookup = dict([])

        # Soma is -1
        self.section_lookup[-1] = self.icell.soma[0]

        for ic, c in enumerate(self.icell.dend):
            self.section_lookup[ic] = c

        # Negative numbers for axon
        for ic, c in enumerate(self.icell.axon):
            self.section_lookup[-ic - 2] = c

    def map_id_to_compartment(self, section_id):

        """
        Map section_id to compartment.

        Neuron_morphology defines sectionID, these must match what this returns
        so that they point to the same compartment.

        Soma is -1
        axons are negative values (currently all set to -2) in Neuron_morphology
        dendrites are 0, 1,2,3,4,5...

        """

        if self.section_lookup is None:
            self.build_section_lookup()

        if isinstance(section_id, int):
            sec = self.section_lookup[section_id]
        else:
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

    def instantiate(self, sim=None, extracellular_regions=None):
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
        self.relocate_NEURON(sim=sim)

        for mechanism in self.mechanisms:
            mechanism.instantiate(sim=sim, icell=self.icell)
        for param in self.params.values():
            param.instantiate(sim=sim, icell=self.icell)

        if self.modulation:
            # TODO: THIS IS NOT THE Simulation object but the NEURON sim object... cant use it!!
            self.modulation.load_json(extracellular_regions=extracellular_regions, neuron_region=self.volume_id)

    ############################################################################

    def get_replacement_axon(self, axon_length=60e-6, axon_diameter=None, axon_nseg=None):

        raise DeprecationWarning("This function is now orphaned. Might need to re-update it if we switch to using hoc")

        """
            If axon_length is given as a scalar, then the code is similar to the BluePyOpt hoc, with the modification
            that the total axon_length is specified by the user.

            If axon_length is a vector, then axon_diameter and axon_nseg must be given and be vectors of same length
            The returned hoc will then create the user specified axon stump.

        """

        if isinstance(axon_length, (int, float, np.integer, np.floating)):
            assert axon_diameter is None and axon_nseg is None
            # We have only the length specified, use two sections

            assert axon_length < 10000e-6, "Please specify replacement axon_length in SI units (meters)."

            axon_length_specified_hoc = \
                f'''
            proc replace_axon(){{ local nSec, D1, D2
              // preserve the number of original axonal sections
              nSec = sec_count(axonal)

              // Try to grab info from original axon
              if(nSec == 0) {{ //No axon section present
                D1 = D2 = 1
              }} else if(nSec == 1) {{
                axon[0] D1 = D2 = diam
              }} else {{
                axon[0] D1 = D2 = diam
                soma distance() //to calculate distance from soma
                forsec axonal{{
                  //if section is longer than 60um then store diam and exit from loop
                  if(distance(0.5) > {axon_length * 1e6}){{
                    D2 = diam
                    break
                  }}
                }}
              }}

              // get rid of the old axon
              forsec axonal{{
                delete_section()
              }}

              create axon[2]

              axon[0] {{
                L = {axon_length * 1e6 / 2}
                diam = D1
                nseg = 1 + 2*int(L/40)
                all.append()
                axonal.append()
              }}
              axon[1] {{
                L = {axon_length * 1e6 / 2}
                diam = D2
                nseg = 1 + 2*int(L/40)
                all.append()
                axonal.append()
              }}
              nSecAxonal = 2
              soma[0] connect axon[0](0), 1
              axon[0] connect axon[1](0), 1
            }}
                    '''

            return axon_length_specified_hoc

        assert axon_length is not None and axon_diameter is not None and axon_nseg is not None

        # In all remaining cases the user has to specify vectors for axon_length, axon_diameter, axon_nseg
        assert len(axon_length) == len(axon_diameter) == len(axon_nseg), \
            f"Unequal lengths: {axon_length = }, {axon_diameter = }, {axon_nseg = }"

        user_defined_axon_hoc = \
            f"""
        proc replace_axon() {{
            
          // get rid of the old axon
          forsec axonal{{
            delete_section()
            
          create axon[{len(axon_length)}]
            
          }}
            """

        for idx, (al, ad, an) in enumerate(zip(axon_length, axon_diameter, axon_nseg)):

            assert 0 < al < 500e-6, "Please make sure you specify axon length in SI units (meters)"
            assert 0 < ad < 10e-6, "Please make sure you specify axon diameter in SI units (meters)"
            assert an % 2 == 1, f"{axon_nseg = } must contain odd integer"

            user_defined_axon_hoc += \
            f"""
          axon[{idx}] {{
            L = {al*1e6}
            diam = {ad*1e6}
            nseg = {an}
            all.append()
            axonal.append()
          }}
          """

        user_defined_axon_hoc += f"""
              nSecAxonal = {len(axon_length)}
        """

        for idx, al in enumerate(axon_length):

            if idx == 0:
                user_defined_axon_hoc += f"""
              soma[0] connect axon[0](0), 1
            """
            else:
                user_defined_axon_hoc += f"""
              axon[{idx-1}] connect axon[{idx}](0), 1
                """

        user_defined_axon_hoc += "}"

        return user_defined_axon_hoc

    def relocate_NEURON(self, sim):

        # TODO: Work in progress

        """ This function can only be called once per neuron, it updates the NEURON coordinates
            to match the Snudda coordinates. """

        if self.rotation is None:
            self.rotation = np.eye(3)

        if self.position is None:
            # Only relocate if a position is given
            return

        assert not self.is_relocated, f"You can only relocate a neuron once"
        self.is_relocated = True

        # Save the old NEURON soma center, subtract it from all coordinates
        # then apply the rotation matrix, and save the NEURON compartment coordinates

        num_points = int(self.icell.soma[0].n3d())
        x_center = sum(self.icell.soma[0].x3d(i) for i in range(num_points)) / num_points
        y_center = sum(self.icell.soma[0].y3d(i) for i in range(num_points)) / num_points
        z_center = sum(self.icell.soma[0].z3d(i) for i in range(num_points)) / num_points

        for sec in self.icell.soma[0].wholetree():
            n_points = int(sec.n3d())
            sim.neuron.h.pt3dstyle(0, sec=sec)  # We don't want the points to translate when parent point is moved

            for i in range(n_points):
                center_pos = np.array([sec.x3d(i) - x_center, sec.y3d(i) - y_center, sec.z3d(i) - z_center])
                rot_pos = np.matmul(self.rotation, center_pos)   # Ska det vara T på pos eller resultatet=
                new_pos = rot_pos + self.position*1e6  # Translated
                sim.neuron.h.pt3dchange(i, new_pos[0], new_pos[1], new_pos[2], sec.diam3d(i), sec=sec)

        #  print(f"Relocated neuron to {self.position}, with rotation {self.rotation}")
        