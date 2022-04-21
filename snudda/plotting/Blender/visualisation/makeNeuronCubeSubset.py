import bpy
import mathutils
import numpy as np
import h5py

from snudda.utils.snudda_path import snudda_parse_path

# The script loads this position file and renders it
# Note that you might need to reposition the camera manually
# if your neurons are centred around a different position

# posFile = "../SmallTest2/network-neuron-positions.hdf5"
#posFile = "/home/hjorth/HBP/StriatumNetwork/model/Article-cube-2160/network-neuron-positions.hdf5"
posFile = "/home/hjorth/HBP/StriatumNetwork/model/Article2019HyperVoxel/network-neuron-positions.hdf5"

fi = h5py.File(posFile,"r")

neuronID = fi["network/neurons/neuronID"].value
pos = fi["network/neurons/position"].value
rot = fi["network/neurons/rotation"].value
morph = fi["network/neurons/morphology"]
name = fi["network/neurons/name"]

# Remove the start cube
bpy.ops.object.delete() 

bpy.data.scenes['Scene'].render.engine = 'CYCLES'
world = bpy.data.worlds['World']
world.use_nodes = True

# changing these values does affect the render.
bg = world.node_tree.nodes['Background']
bg.inputs[0].default_value[:3] = (1.0, 1.0, 1.0)
bg.inputs[1].default_value = 1.0

matMSD1 = bpy.data.materials.new("PKHG")
# matMSD1.diffuse_color = (1.0,0.2,0.2)
matMSD1.diffuse_color = (77./255,151./255,1.0)
matMSD2 = bpy.data.materials.new("PKHG")
# matMSD2.diffuse_color = (0.2,0.2,1.0)
matMSD2.diffuse_color = (67./255,55./255,181./255)

matFS = bpy.data.materials.new("PKHG")
matFS.diffuse_color = (6./255,31./255,85./255)
matChIN = bpy.data.materials.new("PKHG")
matChIN.diffuse_color = (252./255,102./255,0.0)
matLTS = bpy.data.materials.new("PKHG")
matLTS.diffuse_color = (150./255,63./255,212./255)

matOther = bpy.data.materials.new("PKHG")
matOther.diffuse_color = (0.4,0.4,0.4)

objLookup = dict([])

# Possible solution for extra speed:
# https://blender.stackexchange.com/questions/61094/copying-an-object-with-python-without-actually-creating-a-new-object
#

if(True):
  subset = np.random.permutation(2174)[:200]
  nTotal = len(subset)
else:
  subset = None
  nTotal = -1

ctr = 0


for ps,rt,mo,nm,nID in zip(pos,rot,morph,name,neuronID):

  if subset is not None and nID not in subset:
    continue

  ctr += 1
  
  if(mo in objLookup):
    print("Reusing obj for " + str(mo))
    obj = objLookup[mo].copy()
    #obj.data = objLookup[mo].data
    
    for o in objLookup[mo].children:
      oc = o.copy()
      # oc.data = o.data.copy()
      oc.parent = obj

      bpy.context.scene.objects.link(oc)
      
    bpy.context.scene.objects.link(obj)
    
  else:
    print("Loading morphology " + str(mo))
    bpy.ops.import_mesh.swc(filepath=snudda_parse_path(mo))
    obj = bpy.context.selected_objects[0]
    # obj = bpy.data.objects[-1]
    objLookup[mo] = obj

  eRot = mathutils.Matrix(rt.reshape(3,3)).to_euler()    
  obj.rotation_euler = eRot

  # Draw this neuron (the SWC import scales from micrometers to mm), the
  # positions in the simulation are in meters, need to scale it to mm for
  # blender to have same units.
  
  print("Setting position: " + str(ps*1e3))
  obj.location = ps*1e3

  nType = nm.decode().split("_")[0]
  if(nType == "dSPN"):
    print(str(nID) + " dSPN" + " - " + str(ctr) + "/" + str(nTotal))
    mat = matMSD1
  elif(nType == "iSPN"):
    print(str(nID) + " iSPN" + " - " + str(ctr) + "/" + str(nTotal))
    mat = matMSD2
  elif(nType == "FS"):
    print(str(nID) + " FS" + " - " + str(ctr) + "/" + str(nTotal))
    mat = matFS
  elif(nType == "ChIN"):
    print(str(nID) + " ChIN" + " - " + str(ctr) + "/" + str(nTotal))
    mat = matChIN
  elif(nType == "LTS"):
    print(str(nID) + " LTS" + " - " + str(ctr) + "/" + str(nTotal))
    mat = matLTS
  
  else:
    print(str(nID) + " other")
    mat = matOther


  for ch in obj.children:
    ch.active_material = mat

  obj.select = False
  
#import pdb
#pdb.set_trace()

# Add a light source

lampData = bpy.data.lamps.new(name="Sun", type='SUN')
lampObject = bpy.data.objects.new(name="Sun", object_data=lampData)
bpy.context.scene.objects.link(lampObject)

# Place lamp to a specified location
lampObject.location = (1000.0, 1000.0, 1000.0)

# Reposition camera
cam = bpy.data.objects["Camera"]
cam.location = (2.98,2.68,4.96)
cam.rotation_euler = (1.59,0,-0.26)

bpy.context.scene.update()


if(False):
  bpy.ops.render.render()
  bpy.data.images['Render Result'].save_render(filepath="cube-of-neurons-with-neurites.png")

