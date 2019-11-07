""" Based on BluePyOpt exampel code by Werner van Geit, 
modified by Johannes Hjorth """ 



import os
import json
import numpy as np

import bluepyopt.ephys as ephys

class NeuronModel(ephys.models.CellModel):

  def __init__(self,
         cell_name="Unknown",
         morph_file=None,
         mech_file=None,
         param_file=None):

    self.script_dir = os.path.dirname(__file__)
    self.config_dir = os.path.join(self.script_dir, 'config')

    # morph=self.define_morphology(replaceAxon=False,morph_file=morph_file)    
    morph=self.define_morphology(replaceAxon=True,morph_file=morph_file)
    mechs=self.define_mechanisms(mechanism_config=mech_file)
    params=self.define_parameters(param_file)
    
    super(NeuronModel, self).__init__(name=cell_name,morph=morph,
                      mechs=mechs,params=params)      
    self.synlist = []
    self.sectionLookup = None

    #import pdb
    #pdb.set_trace()
    
##############################################################################

  # Helper function

  def define_mechanisms(self, mechanism_config=None):
    """Define mechanisms"""

    assert(mechanism_config is not None)
    # print("Using mechanmism config: " + mechanism_config)
  
    mech_definitions = json.load(open(mechanism_config))

    if("modpath" in mech_definitions):
      mod_path = os.path.join(self.script_dir,mech_definitions["modpath"])
      print("mod_path set to " + mod_path + " (not yet implemented)")
    else:
      mod_path = None


    #import pdb
    #pdb.set_trace()
      
    mechanisms = []
    for sectionlist in mech_definitions:
      channels = mech_definitions[sectionlist]
      
      # This allows us to specify a modpath in the file
      if(sectionlist == "modpath"):
        continue
    
      seclist_loc = \
          ephys.locations.NrnSeclistLocation(sectionlist,
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

  # Helper function

  def define_parameters(self, parameter_config=None):
    """Define parameters"""

    assert(parameter_config is not None)

    # print("Using parameter config: " + parameter_config)
    
    param_configs = json.load(open(parameter_config))
    parameters = []

    for param_config in param_configs:
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

  # Helper function

  def define_morphology(self, replaceAxon=True,morph_file=None):
    """Define morphology. Handles SWC and ASC."""

    assert(morph_file is not None)

    # print("Using morphology: " + morph_file)
    
    return ephys.morphologies.NrnFileMorphology(morph_file, do_replace_axon=replaceAxon)
  # OLD BUGFIX FOR segment pop
  #,replace_axon_hoc=self.getReplacementAxon()) 

##############################################################################

  # Neuron_morphology defines sectionID, these must match what this returns
  # so that they point to the same compartment.
  #
  # Soma is 0
  # axons are negative values (currently all set to -1) in Neuron_morphology
  # dendrites are 1,2,3,4,5... ie one higher than what Neuron internally
  # uses to index the dendrites (due to us wanting to include soma)
  

  def mapIDtoCompartment(self,sectionID):

    if(self.sectionLookup is None):
      
      self.sectionLookup = dict([])

      # Soma is zero
      self.sectionLookup[0] = self.icell.soma[0]

      # Dendrites are consequtive numbers starting from 1
      # Ie neurons dend(0) is in pos 1, dend(99) is in pos 100
      # This so we dont need to special treat soma (pos 0)
      
      for ic,c in enumerate(self.icell.dend):
        self.sectionLookup[ic+1] = c
      
      # Negative numbers for axon
      for ic,c in enumerate(self.icell.axon):
        self.sectionLookup[-ic-1] = c

    try:
      sec = [self.sectionLookup[x] for x in sectionID]
    except:
      print("Missing section ID?")
      print("sectionID = " + str(sectionID))
      print("sectionLookup = " + str(self.sectionLookup))
      import traceback
      tstr = traceback.format_exc()
      print(tstr)
      import pdb
      pdb.set_trace()
        
    return sec

    

############################################################################

  def findDendCompartment(self,synapse_xyz,locType,sim):
    """Locate where on dend sections each synapse is"""

    assert False, "Depricated use mapIDtoCompartment"
    
    dendLoc = []
  
    secLookup = {}
  
    nPoints = 0
    # Find out how many points we need to allocate space for
    for sec in self.icell.dend:
      for seg in sec:
        # There must be a cleaner way to get the 3d points
        # when we already have the section
        nPoints = nPoints + int(sim.neuron.h.n3d(sec=sec))
      
    secPoints = np.zeros(shape=(nPoints,5)) # x,y,z,isec,arclen
    pointCtr = 0

    # Create secPoints with a list of all segment points
    for isec, sec in enumerate(self.icell.dend):
      secLookup[isec] = sec # Lookup table    
      # print("Parsing ", sec)
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


    for row,lType in zip(synapse_xyz,locType): #[locType == 1]:
      # type 1 = axon-dendritic
      if(lType == 1):
        dist = np.sum((secPoints[:,0:3] - row)**2,axis=-1)
        minIdx = np.argmin(dist)

        # Just to double check, the synapse coordinate and the compartment
        # we match it against should not be too far away
        assert dist[minIdx] < 5e-6, \
          "Coords should be close to compartment"
          
        minInfo = secPoints[minIdx,:]
        # Save section and distance within section
        dendLoc.insert(len(dendLoc), [secLookup[minInfo[3]], minInfo[4]])

    
      # axo-somatic synapses (type = 2)
      if(lType == 2):
        dendLoc.insert(len(dendLoc), [self.icell.soma[0], 0.5])

      # For gap junctions (locType == 3) see how Network_simulate.py
      # creates a list of coordinates and calls the function
      if(lType == 4):
        dist = np.sum((secPoints[:,0:3] - row)**2,axis=-1)
        minIdx = np.argmin(dist)

        # Just to double check, the synapse coordinate and the compartment
        # we match it against should not be too far away
        assert(dist[minIdx] < 5e-6)

        minInfo = secPoints[minIdx,:]     
        dendLoc.insert(len(dendLoc),[secLookup[minInfo[3]], minInfo[4]])
    
    # Currently only support axon-dend, and axon-soma synapses, not axon-axon
    # check that there are none in indata
    try:
      assert(all(x == 1 or x == 2 or x == 4 for x in locType))
    except:
      print("Bad locTypes:")
      print("locType: " + str(locType))
      import pdb
      pdb.set_trace()
      
    return dendLoc

  ############################################################################

  # OVERRIDE the create_empty_template from CellModel
  
  @staticmethod
  def create_empty_template(
      template_name,
      seclist_names=None,
      secarray_names=None):
    '''create an hoc template named template_name for an empty cell'''

    objref_str = 'objref this, CellRef, synlist'
    newseclist_str = ''
    
    if seclist_names:
      for seclist_name in seclist_names:
        objref_str += ', %s' % seclist_name
        newseclist_str += \
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
              newseclist_str=newseclist_str,
              create_str=create_str)

    # print(">>>> begin template")
    # print(str(template))
    # print(">>>> end template")
    
    return template

  ############################################################################

  # OVERRIDE the create_empty_template from CellModel

  @staticmethod
  def create_empty_cell(
      name,
      sim,
      seclist_names=None,
      secarray_names=None):
    """Create an empty cell in Neuron"""

    # TODO minize hardcoded definition
    # E.g. sectionlist can be procedurally generated
    #hoc_template = ephys.models.CellModel.create_empty_template(name,
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
      
  def getReplacementAxon(self):

    assert False, "Old bugfix for segmetn stack pop, should not be needed anymore"
    
    newAxonHoc = \
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
    
    return newAxonHoc
