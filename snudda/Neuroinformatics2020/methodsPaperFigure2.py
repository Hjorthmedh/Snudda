import numpy as np
import sys
import os

import matplotlib
import matplotlib.pyplot as plt

from snudda.load import SnuddaLoad

from snudda.analyse import SnuddaAnalyse

class MethodsPaperFigure2(object):

  def __init__(self,networkFileList,legend):

    self.analysis = []

    self.networkFileList = networkFileList
    
    for nf in self.networkFileList:
      print(f"Loading {nf}")
      self.analysis.append(SnuddaAnalyse(hdf5File=nf,
                                         volumeType="cube"))
    
    self.sideLen = 250e-6

    self.legend = legend
    
    self.plotStyle = ['-', '-','--','--',':']
    self.plotColour = [(0,0,0),
                       (0.5,0.5,0.5),
                       (0,0,0),
                       (0.5,0.5,0.5),
                       (0,0,0)]
    
  ############################################################################
    
  def makeConnectionProbabilitySummary(self,
                                       preType,
                                       postType):

    fig,ax = plt.subplots(1)
    matplotlib.rcParams.update({'font.size': 24})

    xMax = 0.0
    yMax = 0.0
    nBins = 86
    dist3D = False

    connectionType = "synapses"
    
    for idx,a in enumerate(self.analysis):

      preID = a.populations[preType]
      postID = a.populations[postType]

      # We can in principle use all pairs, but here we restrict to just the
      # pairs who's post partner are in the centre
      postID = a.getSubPop(volumeType=a.volumeType,
                           volumePart="centre",
                           sideLen=self.sideLen,
                           neuronID=postID)      

      
      (dist,Pcon,countCon,countAll) = \
        a.connectionProbability(preID,postID,nBins,dist3D=dist3D,
                                connectionType=connectionType)

      dHalfStep = (dist[1]-dist[0])/2
      plt.plot((dist+dHalfStep)*1e6 ,Pcon,
               color=self.plotColour[idx],
               linestyle=self.plotStyle[idx],
               linewidth=2)

      try:
        xMax = np.maximum((dist[-1]+dHalfStep)*1e6,xMax)
        yMax = np.maximum(np.nanmax(Pcon),yMax)
      except:

        import traceback
        tstr = traceback.format_exc()
        print(tstr)
        
        import pdb
        pdb.set_trace()
        
    plt.xticks(fontsize=14, rotation=0)
    plt.yticks(fontsize=14, rotation=0)

    labelSize = 22
      
    plt.xlabel("Distance ($\mu$m)",fontsize=labelSize)
    plt.ylabel("Con Prob (%)",fontsize=labelSize)

    # Override max
    xMax = 400 
    
    plt.xlim([0, xMax])
    plt.ylim([0, yMax])

    locs,labels = plt.yticks()
    newLabels = ["{0:g}".format(yy*100) for yy in locs]

    if(locs[0] < 0):
      locs = locs[1:]
      newLabels = newLabels[1:]

    if(locs[-1] > plt.ylim()[1]):
      locs = locs[:-1]
      newLabels = newLabels[:-1]

    plt.yticks(locs,newLabels)

    fontP = matplotlib.font_manager.FontProperties()
    fontP.set_size('xx-small')
    plt.legend(self.legend,prop=fontP)
    
    plt.title(self.analysis[0].neuronName(preType) \
              + " to " + self.analysis[0].neuronName(postType))
    plt.tight_layout()
    plt.ion()
    plt.draw()
    plt.show()
    plt.pause(0.001)

    figName = 'Summary-pruning-dist-dep-connection-probability-' \
               + str(preType) + "-to-" + str(postType) \
                + "-" + str(connectionType)

    self.analysis[0].saveFigure(plt,figName)

  ############################################################################

  def makeNumSynapsesSummaryFigure(self,preType,postType):

      
    fig,ax = plt.subplots(1)
    matplotlib.rcParams.update({'font.size': 24})

    connectionType = "synapses"
    
    for idx,a in enumerate(self.analysis):      

      preID = a.populations[preType]
      postID = a.populations[postType]

      if(connectionType == "synapses"):
        conMat = a.connectionMatrix
      elif(connectionType == "gapjunctions"):
        conMat = a.connectionMatrixGJ
      else:
        print("Unknown connectionType: " + str(connectionType))
        print("Please use 'synapses' or 'gapjunctions'")
        import pdb
        pdb.set_trace()
      
      # We can in principle use all pairs, but here we restrict to just the
      # pairs who's post partner are in the centre
      postID = a.getSubPop(volumeType=a.volumeType,
                           volumePart="centre",
                           sideLen=self.sideLen,
                           neuronID=postID)      

      print("Calculating max synapses")
      maxSynapses = conMat[preID,:][:,postID].max()
      
      # The prune tuning func might set data to 0, we want to exclude those
      meanSynapses = float(np.sum(conMat[preID,:][:,postID].data))\
                     / np.sum(conMat[preID,:][:,postID].data!=0)

      con = conMat[preID,:][:,postID]

      #con = con.toarray()
      #existingCon = con[con != 0]

      existingCon = con[np.nonzero(con)].transpose()

      # Any connections? Otherwise skip plot
      if((type(existingCon) == np.matrixlib.defmatrix.matrix \
          and len(existingCon) == 0)
         or (type(existingCon) != np.matrixlib.defmatrix.matrix \
             and existingCon.getnnz() == 0)):
        continue

      print("Plotting " + str(existingCon.shape[0]) + " connections")


      matplotlib.rcParams.update({'font.size': 22})

      plt.hist(existingCon,
               #range(0,1+maxSynapses),
               range(0,21),
               density=False,
               align="left",
               color=self.plotColour[idx],
               linestyle=self.plotStyle[idx],
               linewidth=2,
               histtype=u"step")

      plt.xlabel("Number of " + connectionType)
      plt.ylabel('Count')

      plt.title(self.analysis[0].neuronName(preType) \
                + " to " + self.analysis[0].neuronName(postType))
  
    plt.tight_layout()

    plt.ion()
    plt.draw()
    plt.show()
    plt.pause(0.001)
    
    figName = "Summary-network-number-of-" + connectionType + "-from-" \
              + preType + "-to-" + postType + "-per-cell"

    self.analysis[0].saveFigure(plt,figName)

  ############################################################################

  def summaryPlotCumDist(self,preType,postType):

    fig,ax = plt.subplots(1)
    matplotlib.rcParams.update({'font.size': 24})

    pair = ("iSPN","dSPN")
    connectionType = "synapses"
    
    for idx,a in enumerate(self.analysis):

      pairID = tuple([a.allTypes.index(x) for x in pair])
      
      cumDist = np.cumsum(a.dendPositionBin[pairID])  \
        /np.sum(a.dendPositionBin[pairID])

      # Dont plot the full range
      endIdx = np.where(a.dendPositionEdges <= 400e-6)[0][-1]

      ax.plot(a.dendPositionEdges[:endIdx]*1e6, cumDist[:endIdx],
              color=self.plotColour[idx],label=self.legend[idx],
              linestyle=self.plotStyle[idx],linewidth=3)

      
    ax.set_xlabel('Distance from soma ($\mu$m)')
    ax.set_ylabel('Cumulative distrib.')

    plt.ion()
    plt.show()
    plt.draw()
    plt.pause(0.0001)

    figName = "Summary-cumDist-of-" + connectionType + "-from-" \
              + preType + "-to-" + postType + "-per-cell"

    self.analysis[0].saveFigure(plt,figName)

      
    
  ############################################################################

if __name__ == '__main__':

  files = ['Neuroinformatics2020/Net10062-var-1/network-pruned-synapses.hdf5',
           'Neuroinformatics2020/Net10062-var-2/network-pruned-synapses.hdf5',
           'Neuroinformatics2020/Net10062-var-3/network-pruned-synapses.hdf5',
           'Neuroinformatics2020/Net10062-var-4/network-pruned-synapses.hdf5',
           'Neuroinformatics2020/Net10062-var-5/network-pruned-synapses.hdf5']

  legends = ['No pruning',
             'DP',
             'DP, f1',
             'DP, f1, SM',
             'DP, f1, SM, mu2']
  
  mpf = MethodsPaperFigure2(files,
                            legends)


  mpf.makeConnectionProbabilitySummary('iSPN','dSPN')

  mpf.makeNumSynapsesSummaryFigure('iSPN','dSPN')

  mpf.summaryPlotCumDist('iSPN','dSPN')
  
  import pdb
  pdb.set_trace()

  
