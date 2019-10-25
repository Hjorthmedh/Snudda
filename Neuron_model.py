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

    assert(mechanism_config is not None)
    # print("Using mechanmism config: " + mechanism_config)
        
    mech_definitions = json.load(
        open(
            os.path.join(
                config_dir,
                mechanism_config)))

    if("modpath" in mech_definitions):
      mod_path = os.path.join(script_dir,mech_definitions["modpath"])
      print("mod_path set to " + mod_path + " (not yet implemented)")
    else:
      mod_path = None

      
    mechanisms = []
    for sectionlist in mech_definitions:

        channels = mech_definitions[sectionlist]
        
        # This allows us to specify a modpath in the file
        if(sectionlist == "modpath"):
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

    assert(parameter_config is not None)

    # print("Using parameter config: " + parameter_config)
        
    param_configs = json.load(open(os.path.join(config_dir, parameter_config)))
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

def define_morphology(replaceAxon=True,morph_file=None):
    """Define morphology. Handles SWC and ASC."""

    assert(morph_file is not None)

    # print("Using morphology: " + morph_file)
        
    return ephys.morphologies.NrnFileMorphology(
        os.path.join(
            script_dir,
            morph_file,
        ),
        do_replace_axon=replaceAxon) 

##############################################################################

def create(param_config=None,morph_file=None,mechanisms_file=None,cell_name="Unknown"):
    """Create cell model"""

    cell = ephys.models.CellModel(
        cell_name,
        morph=define_morphology(replaceAxon=True,morph_file=morph_file),
        mechs=define_mechanisms(mechanism_config=mechanisms_file),
        params=define_parameters(param_config))

    return cell


def findDendCompartment(neuron,synapse_xyz,locType,sim):
    """Locate where on dend sections each synapse is"""

    dendLoc = []
    
    secLookup = {}
    
    nPoints = 0
    # Find out how many points we need to allocate space for
    for sec in neuron.icell.dend:
      for seg in sec:
        # There must be a cleaner way to get the 3d points
        # when we already have the section
        nPoints = nPoints + int(sim.neuron.h.n3d(sec=sec))
            
    secPoints = np.zeros(shape=(nPoints,5)) # x,y,z,isec,arclen
    pointCtr = 0

    # Create secPoints with a list of all segment points
    for isec, sec in enumerate(neuron.icell.dend):
      secLookup[isec] = sec # Lookup table        
      print("Parsing ", sec)
      for seg in sec:
        for i in range(int(sim.neuron.h.n3d(sec=sec))):
          secLen = sim.neuron.h.arc3d(int(sim.neuron.h.n3d(sec=sec)-1), sec=sec)
          # We work in SI units, so convert units from neuron
          secPoints[pointCtr,:] = [sim.neuron.h.x3d(i,sec=sec)*1e-6,
                                   sim.neuron.h.y3d(i,sec=sec)*1e-6,
                                   sim.neuron.h.z3d(i,sec=sec)*1e-6,
                                   isec,
                                   (sim.neuron.h.arc3d(i,sec=sec)/secLen)]
          pointCtr = pointCtr + 1

    # Loop through all (axon-dendritic) synapse locations and find matching compartment
    # type 1 = axon-dendritic
    for row,lType in zip(synapse_xyz,locType): #[locType == 1]:
      if(lType == 1):
        dist = np.sum((secPoints[:,0:3] - row)**2,axis=-1)
        minIdx = np.argmin(dist)

        # Just to double check, the synapse coordinate and the compartment
        # we match it against should not be too far away
        assert(dist[minIdx] < 5e-6)

        minInfo = secPoints[minIdx,:]
        # Save section and distance within section
        dendLoc.insert(len(dendLoc), [secLookup[minInfo[3]], minInfo[4]])

        
      # axo-somatic synapses (type = 2)
      if(lType == 2):
        dendLoc.insert(len(dendLoc), [neuron.icell.soma[0], 0.5])

    # For gap junctions (locType == 3) see how Network_simulate.py
    # creates a list of coordinates and calls the function
      if(lType == 3):
        dist = np.sum((secPoints[:,0:3] - row)**2,axis=-1)
        minIdx = np.argmin(dist)

        # Just to double check, the synapse coordinate and the compartment
        # we match it against should not be too far away
        assert(dist[minIdx] < 5e-6)

        minInfo = secPoints[minIdx,:]       
        dendLoc.insert(len(dendLoc),[secLookup[minInfo[3]], minInfo[4]])
        
    # Currently only support axon-dend, and axon-soma synapses, not axon-axon
    # check that there are none in indata
    assert(all(x == 1 or x == 2 or x == 4 for x in locType))
      
    return dendLoc

# Add a generic synapse to location
#def addSynapse(neuron, location):
    
    
