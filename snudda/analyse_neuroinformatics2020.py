# This script uses the functions defined in Network_analyse.py
#
# Performs analysis of the Striatum
#
#

import numpy as np
import sys
import os

import matplotlib.pyplot as plt

from snudda.analyse import SnuddaAnalyse

class SnuddaAnalyseStriatum(SnuddaAnalyse):

  def __init__(self,simDir,volumeType="cube",sideLen=300e-6):

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
        
    super().__init__(hdf5File=hdf5File,loadCache=True,
                     volumeType=volumeType,
                     sideLen=sideLen)

  ############################################################################

  # Validation. How well do the synapse location for dSPN and iSPN match
  # the experimental data from Straub,..., Sabatini 2016
  
  def plotFSLTScumDist(self,plotFS=True,plotLTS=True):

    pairListList = [[("FSN","dSPN"),("LTS","dSPN")],
                    [("FSN","iSPN"),("LTS","iSPN")]]
    figureNameList = ["synapseCumulativeDistance-FSN-and-LTS-to-dSPN.png",
                      "synapseCumulativeDistance-FSN-and-LTS-to-iSPN.png"]
    figureColourList = [(6./255,31./255,85./255),
                         (150./255,63./255,212./255)]
    fillRange = [[0,100e-6],[50e-6,250e-6]]

    plotFlag = (plotFS,plotLTS)

    assert plotFS or plotLTS, "You must plot either FS or LTS, or both"
    
    if(not plotFS):
      figureNameList = [x.replace("FSN-and-","") for x in figureNameList]
    if(not plotLTS):
      figureNameList = [x.replace("and-LTS-","") for x in figureNameList]
    
    for pairList,figName \
        in zip(pairListList,figureNameList):

      
      plt.rcParams.update({'font.size': 22})      
      fig = plt.figure()
      ax = plt.subplot(111)
      #fig.tight_layout()
      fig.subplots_adjust(bottom=0.15,left=0.15)
      
      for pair,figCol,fillR,plotMeFlag \
          in zip(pairList,figureColourList,fillRange,plotFlag):

        if(not plotMeFlag):
          continue

        try:
          pairID = tuple([self.allTypes.index(x) for x in pair])
        except:
          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          continue
        
        cumDist = np.cumsum(self.dendPositionBin[pairID])  \
                    /np.sum(self.dendPositionBin[pairID])

        # Dont plot the full range
        endIdx = np.where(self.dendPositionEdges <= 400e-6)[0][-1]

        ax.plot(self.dendPositionEdges[:endIdx]*1e6, cumDist[:endIdx],
                color=figCol,label=pair[0],linewidth=3)

        fillIdx = np.where(np.logical_and(fillR[0]
                                          <= self.dendPositionEdges,
                                          self.dendPositionEdges
                                          <= fillR[1]))[0]
        fillStart = fillIdx[0]
        fillEnd = fillIdx[-1]
        
        # Add the area marking
        ax.fill_between(self.dendPositionEdges[fillIdx]*1e6,
                        np.zeros((len(fillIdx),)),
                        cumDist[fillIdx],alpha=0.95,color=figCol,
                        label=None)
        
        ax.set_xlabel('Distance from soma ($\mu$m)')
        ax.set_ylabel('Cumulative distrib.')

        if(plotFS and plotLTS):
          ax.set_title("Synapse locations onto " + pair[1])
        else:
          ax.set_title("Synapses " + pair[0] + " to " + pair[1])

      if(plotFS and plotLTS):
        #Only do legend if both are in figure
        ax.legend(loc="lower right")
        
      plt.ion()
      plt.show()
      plt.draw()
      plt.pause(0.0001)

      self.saveFigure(plt,figName)

      
  ############################################################################
    
if __name__ == "__main__":

  if(len(sys.argv) > 1):
    simDir = sys.argv[1]
    print("Reading network from " + str(simDir))
  else:
    print("Please specify which directory the striatum network files is in")
    exit(-1)

  nas = SnuddaAnalyseStriatum(simDir,volumeType="cube")

  dist3D = False

  yMaxH = None #0.5

  nas.plotIncomingConnections(neuronType="dSPN",preType="iSPN",nBins=20)

  
  nas.plotConnectionProbability("iSPN","dSPN", \
                                dist3D=dist3D, \
                                expMaxDist=[50e-6,100e-6],\
                                expData=[13/47.0,10/80.0],
                                expDataDetailed=[(13,47),(10,80)],
                                yMax=yMaxH)


  nas.plotNumSynapsesPerPair("iSPN","dSPN")
  
  # nas.plotSynapseCumDist()
        

  
