import numpy as np
from .load import SnuddaLoad

class SnuddaExportConnectionMatrix(object):

  def __init__(self,inFile,outFile):
      
    self.sl = SnuddaLoad(inFile)

    self.outFile = outFile
    self.outFileMeta = outFile + "-meta"

    data = self.sl.data
    

    conMat = self.createConMat()
    neuronType = [x["type"] for x in data["neurons"]]
    pos = data["neuronPositions"]

    print("Writing " + self.outFile + " (row = src, column=dest)")
    np.savetxt(self.outFile, conMat, delimiter=",",fmt="%d")

    print("Writing " + self.outFileMeta)
    with open(self.outFileMeta,"w") as fOutMeta:
      for i,(nt,p) in enumerate(zip(neuronType,pos)):
        s = "%d,%s,%f,%f,%f\n" % (i,nt,p[0],p[1],p[2])
        fOutMeta.write(s)
      fOutMeta.close()
        
    #import pdb
    #pdb.set_trace()

  ############################################################################

  def createConMat(self):

    nNeurons = self.sl.data["nNeurons"]

    conMat = np.zeros((nNeurons,nNeurons),dtype=int)
    cnt = 0
    pre,post = 0,0
    
    for synChunk in self.sl.synapse_iterator(data_type="synapses"):
      for syn in synChunk:
        p1 = syn[0]
        p2 = syn[1]

        if(p1 == pre and p2 == post):
          cnt += 1
        else:
          conMat[pre,post] += cnt
          pre = p1
          post = p2
          cnt = 1
          
    conMat[pre,post] += cnt
    cnt = 0

    assert np.sum(np.sum(conMat)) == self.sl.data["nSynapses"], \
      "Synapse numbers in connection matrix does not match"
    
    return conMat
    
    
if __name__ == "__main__":

  from argparse import ArgumentParser

  parser = ArgumentParser(description="Export connection matrix to CSV file")
  parser.add_argument("inFile",help="Snudda HDF5 file with network")
  parser.add_argument("outFile",help="CSV output file")
  args = parser.parse_args()

  secm = SnuddaExportConnectionMatrix(args.inFile,args.outFile)

  
  
