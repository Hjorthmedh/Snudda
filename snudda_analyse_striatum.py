# This script uses the functions defined in Network_analyse.py
#
# Performs analysis of the Striatum
#
#

import numpy as np
import sys
import os



from snudda_analyse import SnuddaAnalyse

class SnuddaAnalyseStriatum(SnuddaAnalyse):

  def __init__(self,simDir):

    if(os.path.isfile(simDir)):
      # We allow the user to also send in a hdf5 file as simDir...
      hdf5File = simDir
      self.simDir = os.path.basename(simDir)
    else:
      self.simDir = simDir
      hdf5File = simDir + "/network-pruned-synapses.hdf5"
    
      if(not os.path.exists(hdf5File)):
        althdf5File = simDir + "/network-connect-voxel-pruned-synapse-file.hdf5"

        if(os.path.exists(althdf5File)):
          hfd5File = althdf5File

    print("Loading " + str(hdf5File))
        
    super().__init__(hdf5File=hdf5File,loadCache=True)

if __name__ == "__main__":

  if(len(sys.argv) > 1):
    simDir = sys.argv[1]
    print("Reading network from " + str(simDir))
  else:
    print("Please specify which directory the striatum network files is in")
    exit(-1)

  nas = SnuddaAnalyseStriatum(simDir)


  #import pdb
  #pdb.set_trace()
  #
  #nas.plotNeurons(0,showSynapses=True)
  
  if(False):
    nas.plotNumSynapsesPerPair("FSN","dSPN")
    nas.plotNumSynapsesPerPair("FSN","iSPN")  
    nas.plotNumSynapsesPerPair("dSPN","dSPN")
    nas.plotNumSynapsesPerPair("dSPN","iSPN")    
    nas.plotNumSynapsesPerPair("iSPN","dSPN")
    nas.plotNumSynapsesPerPair("iSPN","iSPN")    
  
  plotHenrike = True
  plotChIN = True
  plotLTS = True

  dist3D = False
  #dist3D = True

   
  
  if(True):
    nas.plotConnectionProbability("LTS","ChIN", \
                                  dist3D=dist3D )

    # ALSO ADD GAP JUNCTIONS PLOT!!!
    # No exp data for this -- Gittis,...,Kreitzer 2010 (p2228) -- 7/12 (and 3/4 reciprocal) -- distance?
    # FS->FS synapses weaker, 1.1 +/- 1.5nS
    
    nas.plotConnectionProbability("FSN","FSN", \
                                  dist3D=dist3D, \
                                  expMaxDist=[250e-6],\
                                  expData=[7/12.0],
                                  expDataDetailed=[(7,12)] )


    
  # This plots figures for the article

  # 100e-6 from Planert 2010, and 250e-6 data from Gittis 2010
  # 150e-6 from Gittis 2011 (actually 100 +/- 50 micrometers)
  
  # MS <-> MS

  if(plotHenrike):

    yMaxH = 0.5

    nas.plotConnectionProbability("dSPN","iSPN", \
                                  dist3D=dist3D, \
                                  expMaxDist=[50e-6,100e-6],\
                                  expData=[3/47.0,3/66.0],
                                  expDataDetailed=[(3,47),(3,66)],
                                  yMax=yMaxH)
    nas.plotConnectionProbability("dSPN","dSPN", \
                                  dist3D=dist3D, \
                                  expMaxDist=[50e-6,100e-6],\
                                  expData=[5/19.0,3/43.0],
                                  expDataDetailed=[(5,19),(3,43)],
                                  yMax=yMaxH)    
    nas.plotConnectionProbability("iSPN","dSPN", \
                                  dist3D=dist3D, \
                                  expMaxDist=[50e-6,100e-6],\
                                  expData=[13/47.0,10/80.0],
                                  expDataDetailed=[(13,47),(10,80)],
                                  yMax=yMaxH)
    nas.plotConnectionProbability("iSPN","iSPN", \
                                  dist3D=dist3D, \
                                  expMaxDist=[50e-6,100e-6],\
                                  expData=[14/39.0,7/31.0],
                                  expDataDetailed=[(14,39),(7,31)],
                                  yMax=yMaxH)
      
    # FS -> MS
  
    nas.plotConnectionProbability("FSN","iSPN", \
                                  dist3D=dist3D, \
                                  expMaxDist=[100e-6, 150e-6, 250e-6],
                                  expData=[6/9.0, 21/54.0, 27/77.0],
                                  expDataDetailed=[(6,9),(21,54),(27,77)],
                                  yMax=1.0)

    nas.plotConnectionProbability("FSN","dSPN", \
                                  dist3D=dist3D, \
                                  expMaxDist=[100e-6, 150e-6, 250e-6],
                                  expData=[8/9.0, 29/48.0, 48/90.0],
                                  expDataDetailed=[(8,9),(29,48),(48,90)],
                                  yMax=1.0)


  if(plotChIN):
    # "I Janickova et al. 2017 så har de 2018 varicosities i en area på 655 um²,
    # deras slices är 70 um tjocka och om man antar att det inte är några
    # varicositites som täcker varandra så är volym-densiteten/mm³: 4.4*10⁷/mm3"
    # 1.7e6/24*0.01 = 708 ChIN per mm3
    # 4.4e7 / 708 = 62000 varicosities per ChIN
    #
    # 325 ChIN synapser per MS
    # 2-5 ChIN per MS
    # --> 65-160 synapser between a ChIN-MS pair
    # --> Each ChIN connect to 400 - 950 MS
    #
    # Number of MS within 350 micrometer radius 4*pi*(350e-6)^3/3*1.76e6/24e-9
    # --> 13100 MS reachable by ChIN at most (or rather number of MS somas
    # within radius of axonal arbour)
    # -->  3-7% connectivity probability??
     
    nas.plotConnectionProbability("ChIN","iSPN", \
                                  dist3D=dist3D,
                                  expMaxDist=[200e-6],
                                  expData=[62/89.0],
                                  expDataDetailed=[(62,89)],
                                  yMax=1.0)
    nas.plotConnectionProbability("ChIN","dSPN", \
                                  dist3D=dist3D,
                                  expMaxDist=[200e-6],
                                  expData=[62/89.0],
                                  expDataDetailed=[(62,89)],
                                  yMax=1.0)
    nas.plotConnectionProbability("ChIN","FSN", \
                                  dist3D=dist3D,
                                  yMax=1.0)

    # 2-5 ChIN should connect to each MS (approx) --- ref?
    nas.plotIncomingConnections(neuronType="dSPN",preType="ChIN")            
    nas.plotIncomingConnections(neuronType="iSPN",preType="ChIN")            

    # Om vi antar 2000 MS skulle kunna nå varje ChIN, när 10% aktiverade
    # så är 200 MS aktiva, om 75% av ChIN känner av MS input
    # (1-p)^200 = 0.25 --> 0.7 %
    nas.plotConnectionProbability("dSPN","ChIN", \
                                  dist3D=dist3D)
    nas.plotConnectionProbability("iSPN","ChIN", \
                                  dist3D=dist3D)

    # A MS neuron receives 1e4 assymetrical synapses (Kincaid et al 1998),
    # and 2500 symmetrical synapses (Ingham et al 1998). Symmetrical synapses
    # can be dopaminergic, cholinergic or GABAergic, with dopaminergic
    # being 13% (Roberts et al 2002). Assuming that cholinergic inputs are
    # a similar percentage, 650 symmetrical synapses per MS are not GABAergic.
    #
    # --> 0.13*2500 = 325 ChIN inputs to MS
    nas.plotIncomingConnections(neuronType="ChIN",preType="dSPN")
    nas.plotIncomingConnections(neuronType="ChIN",preType="iSPN")  

  if(True):
    nas.plotConnectionProbability("FSN","FSN", \
                                  dist3D=dist3D ,
                                  connectionType="gapjunctions",
                                  expMaxDist=[250e-6,250e-6],
                                  expData=[2/6.0,3/7.0],
                                  expDataDetailed=[(2,6),(3,7)],)

    nas.plotNumSynapsesPerPair("FSN","FSN",connectionType="gapjunctions")
   
    nas.plotIncomingConnections(neuronType="FSN",preType="FSN",
                                connectionType="gapjunctions")


    
      
  if(plotLTS):

    # 3/21 LTS->MS, Basal Ganglia book --- distance??
    # Ibanez-Sandoval, ..., Tepper  2011 3/21 -- if patching around visual axon
    # but 2/60 when patching blind
    nas.plotConnectionProbability("LTS","dSPN", \
                                  dist3D=dist3D,
                                  expMaxDist=[250e-6],
                                  expData=[2/60.0],
                                  expDataDetailed=[(2,60)])

    nas.plotConnectionProbability("LTS","iSPN", \
                                  dist3D=dist3D,
                                  expMaxDist=[250e-6],
                                  expData=[2/60.0],
                                  expDataDetailed=[(2,60)])

  
    # Silberberg et al 2013, 2/12 FS-> LTS connected --- distance??
    nas.plotConnectionProbability("FSN","LTS", \
                                  dist3D=dist3D,
                                  expMaxDist=[250e-6],
                                  expData=[2.0/12],
                                  expDataDetailed=[(2,12)])
    
    nas.plotConnectionProbability("ChIN","LTS", \
                                  dist3D=dist3D)
    nas.plotConnectionProbability("ChIN","iSPN", \
                                  dist3D=dist3D)
    nas.plotConnectionProbability("ChIN","dSPN", \
                                  dist3D=dist3D)


    
    nas.nearestPreNeighbourDistance("LTS","dSPN")
    nas.nearestPreNeighbourDistance("LTS","iSPN")

    
