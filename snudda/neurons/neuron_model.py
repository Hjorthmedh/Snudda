""" Based on BluePyOpt exampel code by Werner van Geit, 
modified by Johannes Hjorth """

import os
import json
import numpy as np

import bluepyopt.ephys as ephys

script_dir = os.path.dirname(__file__)
config_dir = os.path.join(script_dir, 'config')


##############################################################################


def define_mechanisms(mechanism_config=None):
    """Define mechanisms"""

    assert (mechanism_config is not None)
    # print("Using mechanmism config: " + mechanism_config)

    with open(os.path.join(config_dir, mechanism_config), "r") as f:
        mech_definitions = json.load(f)

    if "modpath" in mech_definitions:
        mod_path = os.path.join(script_dir, mech_definitions["modpath"])
        print(f"mod_path set to {mod_path} (not yet implemented)")
    else:
        mod_path = None

    mechanisms = []
    for sectionlist in mech_definitions:

        channels = mech_definitions[sectionlist]

        # This allows us to specify a modpath in the file
        if sectionlist == "modpath":
            continue

        seclist_loc = ephys.locations.NrnSeclistLocation(
            sectionlist,
            seclist_name=sectionlist)
        for channel in channels:
            mechanisms.append(ephys.mechanisms.NrnMODMechanism(
                name='%s.%s' % (channel, sectionlist),
                mod_path=mod_path,
                suffix=channel,
                locations=[seclist_loc],
                preloaded=True))

    return mechanisms


##############################################################################


def define_parameters(parameter_config=None):
    """Define parameters"""

    assert (parameter_config is not None)

    # print("Using parameter config: " + parameter_config)

    with open(os.path.join(config_dir, parameter_config), "r") as f:
        param_configs = json.load(f)
    parameters = []

    for param_config in param_configs:
        if 'value' in param_config:
            frozen = True
            value = param_config['value']
            bounds = None
        elif 'bounds':
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
            elif param_config['dist_type'] == 'exp':
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

def define_morphology(replace_axon=True, morph_file=None):
    """Define morphology. Handles SWC and ASC."""

    assert (morph_file is not None)

    # print("Using morphology: " + morph_file)

    return ephys.morphologies.NrnFileMorphology(
        os.path.join(
            script_dir,
            morph_file,
        ),
        do_replace_axon=replace_axon)


##############################################################################

def create(param_config=None, morph_file=None, mechanisms_file=None, cell_name="Unknown"):
    """Create cell model"""

    cell = ephys.models.CellModel(
        cell_name,
        morph=define_morphology(replace_axon=True, morph_file=morph_file),
        mechs=define_mechanisms(mechanism_config=mechanisms_file),
        params=define_parameters(param_config))

    return cell


def find_dend_compartment(neuron, synapse_xyz, loc_type, sim):
    """Locate where on dend sections each synapse is"""

    dend_loc = []

    sec_lookup = {}

    n_points = 0
    # Find out how many points we need to allocate space for
    for sec in neuron.icell.dend:
        for seg in sec:
            # There must be a cleaner way to get the 3d points
            # when we already have the section
            n_points = n_points + int(sim.neuron.h.n3d(sec=sec))

    sec_points = np.zeros(shape=(n_points, 5))  # x,y,z,isec,arclen
    point_ctr = 0

    # Create secPoints with a list of all segment points
    for isec, sec in enumerate(neuron.icell.dend):
        sec_lookup[isec] = sec  # Lookup table
        print("Parsing ", sec)

        # TODO: Check if the seg loop is redundant!!
        for seg in sec:
            for i in range(int(sim.neuron.h.n3d(sec=sec))):
                sec_len = sim.neuron.h.arc3d(int(sim.neuron.h.n3d(sec=sec) - 1), sec=sec)
                # We work in SI units, so convert units from neuron
                sec_points[point_ctr, :] = [sim.neuron.h.x3d(i, sec=sec) * 1e-6,
                                            sim.neuron.h.y3d(i, sec=sec) * 1e-6,
                                            sim.neuron.h.z3d(i, sec=sec) * 1e-6,
                                            isec,
                                            (sim.neuron.h.arc3d(i, sec=sec) / sec_len)]
                point_ctr = point_ctr + 1

    # Loop through all (axon-dendritic) synapse locations and find matching compartment
    # type 1 = axon-dendritic
    for row, l_type in zip(synapse_xyz, loc_type):  # [locType == 1]:
        if l_type == 1:
            dist = np.sum((sec_points[:, 0:3] - row) ** 2, axis=-1)
            min_idx = np.argmin(dist)

            # Just to double check, the synapse coordinate and the compartment
            # we match it against should not be too far away
            assert (dist[min_idx] < 5e-6)

            min_info = sec_points[min_idx, :]
            # Save section and distance within section
            dend_loc.insert(len(dend_loc), [sec_lookup[min_info[3]], min_info[4]])

        # axo-somatic synapses (type = 2)
        if l_type == 2:
            dend_loc.insert(len(dend_loc), [neuron.icell.soma[0], 0.5])

        # For gap junctions (locType == 3) see how Network_simulate.py
        # creates a list of coordinates and calls the function
        if l_type == 3:
            dist = np.sum((sec_points[:, 0:3] - row) ** 2, axis=-1)
            min_idx = np.argmin(dist)

            # Just to double check, the synapse coordinate and the compartment
            # we match it against should not be too far away
            assert (dist[min_idx] < 5e-6)

            min_info = sec_points[min_idx, :]
            dend_loc.insert(len(dend_loc), [sec_lookup[min_info[3]], min_info[4]])

    # Currently only support axon-dend, and axon-soma synapses, not axon-axon
    # check that there are none in indata
    assert (all(x == 1 or x == 2 or x == 4 for x in loc_type))

    return dend_loc

# Add a generic synapse to location
# def addSynapse(neuron, location):
