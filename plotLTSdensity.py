import numpy as np



# fileName = "morphology/network/v7/LTS-20190212/LTS_papercorrected_savedcorrected_axon.swc"
fileName = "cellspecs/lts/LTS_Experiment-9862_20181211/Experiment-9862corrected-cor-rep.swc-withAxon"

# Load SWC file
with open(fileName,"r") as f:
  swcLines = f.readlines()    

binWidth = 10e-6
dendArborHist = np.zeros((100,))
radius = np.arange(0,dendArborHist.shape[0])*binWidth + binWidth/2
  

swcVals = np.zeros(shape=(len(swcLines),7))

comp_type = { 1: "soma", 2: "axon", 3: "dend", 4: "apic" }

# Parse morphology
nComps = 0
for line in swcLines:
  if(line[0] != "#"):
    swcVals[nComps,:] = [float(s) for s in line.split()]
    nComps += 1

swcVals = swcVals[:nComps,:]
    
# Subtract 1 from ID and parentID, so we get easier indexing
swcVals[:,0] -= 1
swcVals[:,6] -= 1

swcVals[:,2:6] *= 1e-6 # Convert to meter x,y,z, radie

# Plot neuron
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#plt.plot(swcVals[:,2],swcVals[:,3],swcVals[:,4],'*')

assert swcVals[0,1] == 1, "First row should be soma"

somaCoord = swcVals[0,2:5]

for row in swcVals:
  if row[1] == 2:
    col = "r"
  elif(row[1] in [3,4]):
    col = "b"
  else:
    continue
    
  parentCoords = swcVals[int(row[6]),2:5] - somaCoord
  coords = row[2:5] - somaCoord
  centre = (coords + parentCoords) / 2.0
    
  plt.plot([coords[0],parentCoords[0]],
           [coords[1],parentCoords[1]],
           [coords[2],parentCoords[2]], col + '-')





########################################################################

nPoints = 100000
rng = [-100e-6,1500e-6,-200e-6,200e-6,-200e-6,200e-6]

x = rng[0] + (rng[1] - rng[0]) * np.random.rand(nPoints,1)
y = rng[2] + (rng[3] - rng[2]) * np.random.rand(nPoints,1)
z = rng[4] + (rng[5] - rng[4]) * np.random.rand(nPoints,1)

x1 = 200e-6 # Inner bump
xw1 = 100e-6 #100e-6
x2 = 300e-6 # Tube
xw2 = 300e-6
x3 = 700e-6 # Outer bump
xw3 = 100e-6

Pkeep = 0.25*np.exp(-(((x-x1)/xw1)**2 + ((y-0)/50e-6)**2 + ((z-0)/30e-6)**2)) \
    + 1*np.exp(-(((x-x2)/xw2)**2 + ((y-0)/15e-6)**2 + ((z-0)/10e-6)**2)) \
+ 1*np.exp(-(((x-x3)/xw3)**2 + ((y-0)/15e-6)**2 + ((z-0)/15e-6)**2))

Pkeep = Pkeep / np.max(Pkeep)

idx = np.where(Pkeep > np.random.rand(nPoints,1))[0]
xkeep = x[idx]
ykeep = y[idx]
zkeep = z[idx]

dall = np.sqrt(xkeep**2 + ykeep**2 + zkeep**2)

#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

#fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

ax.scatter(0,0,0,'*',c="g")
ax.scatter(xkeep,ykeep,zkeep,'*',c='k')
ax.set_aspect("equal")


#ax.set_xlim([rng[0],rng[1]])

plt.ion()
plt.show()

plt.savefig("figures/LTS-axon-density-tuning-3D.png")

# Read the data from Tepper
tepperData = "DATA/Tepper2011extracted.csv"
teppersholldata = np.genfromtxt(tepperData,delimiter=",")

xsholl = teppersholldata[1:,0]
axonsholl = teppersholldata[1:,1]
dendsholl = teppersholldata[1:,2]



fig,ax = plt.subplots()
ax.hist(dall)
#plt.scatter(xsholl*1e-6, axonsholl,c="r")
ax2 = ax.twinx()

ax2.plot(xsholl*1e-6, axonsholl,"*",c="r")
ax.set_xlabel("Distance from soma")
ax.set_ylabel("Axon count")
ax2.set_ylabel("Axon sholl count")

plt.ion()
plt.show()

plt.savefig("figures/LTS-axon-density-tuning-histogram.png")

import pdb
pdb.set_trace()
