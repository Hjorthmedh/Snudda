# import bpy

import h5py
import numpy as np

# Start by loading the positions from file

fi = h5py.File("../FullStriatumPos/network-neuron-positions.hdf5","r")

xm = 0.005
ym = 0.0035
zm = 0.0027
wm = 0.0003

keepIdx = []

pos = fi["network/neurons/position"]
nColour = []

for ir,row in enumerate(pos):
  if(ir % 50000 == 0):
    print(str(ir) + "/" + str(pos.shape[0]))
  if(xm <= row[0] and row[0] <= xm + wm \
     and ym <= row[1] and row[1] <= ym + wm \
     and zm <= row[2] and row[2] <= zm + wm):
    keepIdx.append(ir)

    nType = fi["network/neurons/name"][ir].decode().split("_")[0]
    if(nType == "MSD1"):
      col = [1.0, 0.2, 0.2]
    elif(nType == "MSD2"):
      col = [0.2, 0.2, 1.0]
    elif(nType == "LTS"):
      col = [0.7, 0.0, 0.7]
    elif(nType == "FS"):
      col = [0.1, 0.1, 0.1]
    elif(nType == "ChIN"):
      col = [0.9, 0.7, 0.1]      
    else:
      col = [0.4, 0.4, 0.4]
      
    nColour.append(col)
    
insideCoords = pos[keepIdx,:]
neuronColours = np.array(nColour)
posCol = np.concatenate((insideCoords,neuronColours),axis=1)
np.savetxt("hypervoxelcoords2.txt",posCol,delimiter=",")

print("Next render... drawHyperCubeSomasRender.py")

import pdb
pdb.set_trace()
      
