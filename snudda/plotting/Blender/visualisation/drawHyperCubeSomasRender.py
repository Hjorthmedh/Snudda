import bpy
import numpy as np

poscol = np.genfromtxt('/home/hjorth/HBP/StriatumNetwork/model/article/hypervoxelcoords.txt', delimiter=',')
# poscol = np.genfromtxt('hypervoxelcoords2.txt', delimiter=',')


# Remove the start cube
bpy.ops.object.delete()


bpy.data.scenes['Scene'].render.engine = 'CYCLES'
world = bpy.data.worlds['World']
world.use_nodes = True

# changing these values does affect the render.
bg = world.node_tree.nodes['Background']
bg.inputs[0].default_value[:3] = (1.0, 1.0, 1.0)
bg.inputs[1].default_value = 1.0

mats = dict([])

for row in poscol:
  x = row[0]*1e4
  y = row[1]*1e4
  z = row[2]*1e4
  col = row[3:6]
  tcol = tuple(col)

  if(tcol in mats):
    m = mats[tcol]
  else:
    m = bpy.data.materials.new("PKHG")  
    m.diffuse_color = col
    mats[tcol] = m

  bpy.ops.mesh.primitive_uv_sphere_add(location=(x,y,z),size=0.01*10)
  o = bpy.context.selected_objects[0] 
  o.active_material = m

# Add a light source

lampData = bpy.data.lamps.new(name="Sun", type='SUN')
lampObject = bpy.data.objects.new(name="Sun", object_data=lampData)
bpy.context.scene.objects.link(lampObject)


# Place lamp to a specified location
lampObject.location = (5.0, 5.0, 15.0)

bpy.context.scene.camera.data.clip_end = 10000


# Reposition camera
cam = bpy.data.objects["Camera"]
#cam.location = (55.965,24.934,36.768)
cam.location = (8.89, 5.23, 3.51 )
cam.rotation_euler = (71.2*np.pi/180,0,-239*np.pi/180)
  
