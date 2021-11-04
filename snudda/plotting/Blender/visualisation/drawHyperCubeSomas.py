# import bpy

import h5py
import numpy as np

# Start by loading the positions from file

#fi = h5py.File("../FullStriatumPos/network-neuron-positions.hdf5","r")
#fi = h5py.File("/home/hjorth/HBP/StriatumNetwork/model/Article-cube-2160/network-neuron-positions.hdf5","r")

# 2174 neurons, density 80500... 0.3x0.3x0.3
# hFile = "/home/hjorth/HBP/StriatumNetwork/model/Article2019HyperVoxel/network-neuron-positions.hdf5"
hFile = "/home/hjorth/HBP/Snudda/snudda/examples/networks/Net10062/network-neuron-positions.hdf5"
fi = h5py.File(hFile,"r")

keepAll = True
centreNeurons = True

xm = 0.005
ym = 0.0035
zm = 0.0027
wm = 0.0003

keepIdx = []

print("Loaded")

pos = fi["network/neurons/position"]
nColour = []

centre = np.mean(pos,axis=0)

for ir,row in enumerate(pos):
  if(ir % 50000 == 0):
    print(str(ir) + "/" + str(pos.shape[0]))
    
  if(keepAll or (xm <= row[0] and row[0] <= xm + wm \
                 and ym <= row[1] and row[1] <= ym + wm \
                 and zm <= row[2] and row[2] <= zm + wm)):
    keepIdx.append(ir)

    nType = fi["network/neurons/name"][ir].decode().split("_")[0]
    if(nType == "dSPN" or nType == "MSD1"):
      # col = [1.0, 0.2, 0.2]
      col = (77./255,151./255,1.0)
    elif(nType == "iSPN" or nType == "MSD2"):
      #col = [0.2, 0.2, 1.0]
      col = (67./255,55./255,181./255)
    elif(nType == "FS"):
      col = (6./255,31./255,85./255)
    elif(nType == "ChIN"):
      col = (252./266,102./255,0.0)
    elif(nType == "LTS"):
      col = (150./255,63./255,212./255)
    else:
      col = [0.4, 0.4, 0.4]
      
    nColour.append(col)

if(centreNeurons):
  insideCoords = pos[keepIdx,:] - centre
else:
  insideCoords = pos[keepIdx,:]
neuronColours = np.array(nColour)
posCol = np.concatenate((insideCoords,neuronColours),axis=1)
np.savetxt("hypervoxelcoords.txt",posCol,delimiter=",")

print("Wrote " + str(len(keepIdx)) + " neurons")

print("Next render... drawHyperCubeSomasRender.py")

#import pdb
#pdb.set_trace()
      
