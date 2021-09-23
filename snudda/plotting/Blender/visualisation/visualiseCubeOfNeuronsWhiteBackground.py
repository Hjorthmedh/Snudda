# This script visualises two neurons and the location of the synapses
# connecting them.

print("!!! THIS CODE IS REALLY SLOW !!! -- use makeNeuronCube.py if you only want neurons without synapses")

import bpy
import mathutils
import numpy as np
import h5py
from snudda.utils.snudda_path import snudda_parse_path

# networkFile = "/home/hjorth/HBP/StriatumNetwork/model/SmallTest2/network-connect-voxel-pruned-synapse-file.hdf5"
networkFile = "/home/hjorth/HBP/StriatumNetwork/model/SmallTest2/network-connect-voxel-pruned-synapse-file.hdf5"

outFile = "cubeOfNeurons-synapses-white.blend"

visualiseID = range(0,1800) #[20,1217]

# Load the neuron positions
fi = h5py.File(networkFile,"r")

neuronID = fi["network/neurons/neuronID"].value
keepIdx = np.array([(x in visualiseID) for x in neuronID])
neuronID = neuronID[keepIdx]

pos = fi["network/neurons/position"][keepIdx,:]
rot = fi["network/neurons/rotation"][keepIdx,:]
morph = fi["network/neurons/morphology"][keepIdx]
name = fi["network/neurons/name"][keepIdx]

origo = fi["meta/simulationOrigo"].value
voxelSize = fi["meta/voxelSize"].value


synapses = fi["network/synapses"]


# Remove the start cube
bpy.ops.object.delete() 

bpy.data.scenes['Scene'].render.engine = 'CYCLES'
world = bpy.data.worlds['World']
world.use_nodes = True

whiteBackground = True

# changing these values does affect the render.
bg = world.node_tree.nodes['Background']
if(whiteBackground): # Set to True for white background
  bg.inputs[0].default_value[:3] = (1.0, 1.0, 1.0)
  bg.inputs[1].default_value = 1.0
else:
  bg.inputs[0].default_value[:3] = (0.0, 0.0, 0.0)
  bg.inputs[1].default_value = 0.0

  
matMSD1 = bpy.data.materials.new("PKHG")
#matMSD1.diffuse_color = (1.0,0.2,0.2)
matMSD1.diffuse_color = (77./255,151./255,1.0)

matMSD2 = bpy.data.materials.new("PKHG")
#matMSD2.diffuse_color = (0.2,0.2,1.0)
matMSD2.diffuse_color = (67./255,55./255,181./255)

matFS = bpy.data.materials.new("PKHG")
matFS.diffuse_color = (6./255,31./255,85./255)
matChIN = bpy.data.materials.new("PKHG")
matChIN.diffuse_color = (252./255,102./255,0.0)
matLTS = bpy.data.materials.new("PKHG")
matLTS.diffuse_color = (150./255,63./255,212./255)

matOther = bpy.data.materials.new("PKHG")
matOther.diffuse_color = (0.4,0.4,0.4)

matSynapse = bpy.data.materials.new("PKHG")
if(whiteBackground):
  matSynapse.diffuse_color = (0.8,0.0,0.0)
else:
  matSynapse.diffuse_color = (1.0,1.0,0.9)
  
#matSynapse.use_transparency = True
matSynapse.use_nodes = True
#import pdb
#pdb.set_trace()


if(not whiteBackground):
  emissionStrength = 5.0
  
  # Make synapses glow
  emission = matSynapse.node_tree.nodes.new('ShaderNodeEmission')
  emission.inputs['Strength'].default_value = emissionStrength

  material_output = matSynapse.node_tree.nodes.get('Material Output')
  matSynapse.node_tree.links.new(material_output.inputs[0], emission.outputs[0])


for ps,rt,mo,nm,nID in zip(pos,rot,morph,name,neuronID):

  eRot = mathutils.Matrix(rt.reshape(3,3)).to_euler()
  bpy.ops.import_mesh.swc(filepath=snudda_parse_path(mo))
  obj = bpy.context.selected_objects[0]
  
  obj.rotation_euler = eRot
  print("Setting position: " + str(ps*1e3))
  obj.location = ps*1e3

  nType = nm.decode().split("_")[0]
  if(nType == "dSPN" or nType == "MSD1"):
    print(str(nID) + " dSPN")
    mat = matMSD1
  elif(nType == "iSPN" or nType == "MSD2"):
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

  obj.select = False


# Draw the synapses

drawSynapses = False

if(drawSynapses):

  nSynapses = 0

  for ob in bpy.context.selected_objects:
    ob.select = False

  for syn in synapses:
    preID = syn[0]
    postID = syn[1]

    if(preID in visualiseID and postID in visualiseID):
      # Draw this neuron (the SWC import scales from micrometers to mm), the
      # positions in the simulation are in meters, need to scale it to mm for
      # blender to have same units.
      x = (origo[0] + voxelSize*syn[2]) * 1e3
      y = (origo[1] + voxelSize*syn[3]) * 1e3
      z = (origo[2] + voxelSize*syn[4]) * 1e3

      #bpy.ops.mesh.primitive_cube_add(location=(x,y,z))
      #bpy.context.scene.objects.active.dimensions = \
        #  (voxelSize, voxelSize, voxelSize)

      bpy.ops.mesh.primitive_uv_sphere_add(location=(x,y,z),size=0.001*4)
      ob = bpy.context.selected_objects[0] 
      ob.active_material = matSynapse
      ob.select = False
      
      nSynapses += 1
    
      print("Added synapse at " + str([x,y,z]))

  print("nSynapses = " + str(nSynapses))
    
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
#cam.location = (2.98,2.68,4.96)
#cam.rotation_euler = (1.59,0,-0.26)

cam.location = (3.19,3.46,4.96)
cam.rotation_euler = (96.7*np.pi/180,0,-14.3*np.pi/180)


if(True):
  bpy.ops.render.render()
  bpy.data.images['Render Result'].save_render(filepath="cube-of-neurons-with-neurites.png")

bpy.ops.wm.save_as_mainfile(filepath=outFile)
