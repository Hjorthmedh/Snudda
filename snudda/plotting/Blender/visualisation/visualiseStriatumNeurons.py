# blender -b -P visualiseStriatumNeurons.py

import bpy
import mathutils
import numpy as np
import h5py
from snudda.utils.snudda_path import snudda_parse_path

# The script loads this position file and renders it
# Note that you might need to reposition the camera manually
# if your neurons are centred around a different position

# posFile = "../SmallTest2/network-neuron-positions.hdf5"
if(False):
   posFile = "/home/hjorth/HBP/StriatumNetwork/model/Visualise17300/network-neuron-positions.hdf5"
   outFile = "Visualise17300.blend"
elif(False):
   posFile = "/home/hjorth/HBP/StriatumNetwork/model/Visualise1000/network-neuron-positions.hdf5"
   outFile = "Visualise1000.blend"
else:
   posFile = "/home/hjorth/HBP/StriatumNetwork/model/Vis100/network-neuron-positions.hdf5"
   outFile = "Visualise100.blend"

   #posFile = "/home/hjorth/HBP/StriatumNetwork/model/Vis30/network-neuron-positions.hdf5"
   #outFile = "Visualise30.blend"

pngFile = outFile.replace(".blend",".png")
assert pngFile != outFile, "Something went wrong with renaming"

   
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

addSphere = []

for ps,rt,mo,nm,nID in zip(pos,rot,morph,name,neuronID):

  if(mo in objLookup):
    print("Reusing obj for " + str(mo))
    templateObj = objLookup[mo]
    
    obj = templateObj.copy()

    for ch in templateObj.children:
      chCopy = ch.copy()
      chCopy.parent = obj
      bpy.context.scene.objects.link(chCopy)       

    if(templateObj.data is not None):
       obj.data = templateObj.data.copy()
       
    obj.animation_data_clear()
    bpy.context.scene.objects.link(obj)
      
    
  else:
    print("Loading morphology " + str(mo))
    bpy.ops.import_mesh.swc(filepath=snudda_parse_path(mo))
    # obj = bpy.data.objects[-1]
    obj = bpy.context.selected_objects[0] 

    # print("Morphology object: " + str(obj))
    objLookup[mo] = obj
      
  eRot = mathutils.Matrix(rt.reshape(3,3)).to_euler()    
  obj.rotation_euler = eRot
     
  # Draw this neuron (the SWC import scales from micrometers to mm), the
  # positions in the simulation are in meters, need to scale it to mm for
  # blender to have same units.
  
  print("Setting position: " + str(ps*1e3))
  #print("Offset: " + str(-templateOffset*1e3))
  #print("Final location: " + str((ps-templateOffset)*1e3))
  
  obj.location = ps*1e3
  # somaObj.location = ps*1e3
  
  # Loop through all children also
  #for ch in obj.children:
  #   ch.rotation_euler = eRot
  #   ch.location = (ps-templateOffset)*1e3

  
  nType = nm.decode().split("_")[0]
  if(nType == "dSPN"):
    print(str(nID) + " dSPN")
    mat = matMSD1
  elif(nType == "iSPN"):
    print(str(nID) + " iSPN")
    mat = matMSD2
  elif(nType == "FS"):
    print(str(nID) + " FS")
    mat = matFS
  elif(nType == "ChIN"):
    print(str(nID) + " ChIN")
    mat = matChIN
  elif(nType == "LTS"):
    print(str(nID) + " LTS")
    mat = matLTS
  
  else:
    print(str(nID) + " other")
    mat = matOther


  for ch in obj.children:
    ch.active_material = mat

  # somaObj.active_material = mat

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
# cam.location = (2.98,2.68,4.96) -- zoomed out
cam.location = (3.23,3.62,4.98)
# cam.location = (-0.0631,-8.76,4.733)
cam.rotation_euler = (1.59,0,-0.26)

try:
  bpy.data.scenes['Scene'].render.resolution_x = 3000
  bpy.data.scenes['Scene'].render.resolution_y = 3000
except:
   print("Something wrong with setting render resolution")
   import pdb
   pdb.set_trace()


if(True):
  bpy.ops.render.render()
  bpy.data.images['Render Result'].save_render(filepath=pngFile)

  
bpy.ops.wm.save_as_mainfile(filepath=outFile)
