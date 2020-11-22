import os
import h5py
from snudda.load import SnuddaLoad

class InspectInput(object):

  def __init__(self,networkFile,inputFile):

    self.networkFile = networkFile
    self.inputFile = inputFile

    # We just need this to identify which neuron is which
    self.network = SnuddaLoad(self.networkFile, load_synapses=False)

    self.inputData = h5py.File(inputFile,'r')


  def getMorphologies(self):

    return [os.path.basename(c["morphology"]) \
            for c in self.network.data["neurons"]]


  def getInputTypes(self,cellID):

    sCellID = [str(c) for c in cellID]
    inputTypes = set()
    
    for scID in sCellID:
      inpT = set([inp for inp in self.inputData["input"][scID]])
      inputTypes = inputTypes.union(inpT)
      
    return list(inputTypes)

  
  def checkInputRatio(self, neuronType, verbose=True):

    print(f"Counting inputs for {neuronType}")
    
    cellID = self.network.get_cell_id_of_type(neuronType)
    cellIDstr = [str(c) for c in cellID]
    
    inputCount = dict()

    cellMorph = self.getMorphologies()

    uniqueMorph = set([cellMorph[c] for c in cellID])

    inputTypeList = self.getInputTypes(cellID)
    
    for inp in inputTypeList:
      inputCount[inp] = dict()
      for um in uniqueMorph:
        inputCount[inp][um] = 0

    morphCounter = dict()

    for um in uniqueMorph:
      morphCounter[um] = 0
        
      
    # !!! TODO: We should split this by morphology also...
      
    for cID in self.inputData['input']:
      if cID in cellIDstr:
        morph = cellMorph[int(cID)]
        morphCounter[morph] += 1
                          
        for inp in inputTypeList:
          if inp in self.inputData['input'][cID]:
            nInput = len(self.inputData['input'][cID][inp]['nSpikes'])
            inputCount[inp][morph] += nInput

    for inp in inputTypeList:
      for um in uniqueMorph:
        avgInp = round(inputCount[inp][um]/morphCounter[um],1)
        print(f"{inp} morphology {um}: {avgInp} inputs")
        
    return inputCount


if __name__ == '__main__':

  from argparse import ArgumentParser

  parser = ArgumentParser(description="Inspect input")
  parser.add_argument("networkFile", help="Network file (hdf5)",type=str)
  parser.add_argument("inputFile", help="Input file (hdf5)",type=str)
  parser.add_argument("neuronType", help="Neuron type",type=str)

  args = parser.parse_args()

  inspector = InspectInput(networkFile=args.networkFile,
                           inputFile=args.inputFile)

  inspector.checkInputRatio(neuronType=args.neuronType)
