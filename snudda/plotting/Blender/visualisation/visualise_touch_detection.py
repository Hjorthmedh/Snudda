# This script visualises two neurons and the location of the synapses
# connecting them.

import bpy
import mathutils
import numpy as np
import h5py
from snudda.utils.snudda_path import snudda_parse_path, get_snudda_data
import os

#network_file = "/home/hjorth/HBP/example/network/delme/network-synapses.hdf5"
network_file = "/home/hjorth/HBP/example/network/delme2/network-synapses.hdf5"
network_path = os.path.basename(network_file)

snudda_data = get_snudda_data(network_path=network_path)

out_file = "two-neurons-synapses.blend"

visualise_id = [0, 3]

# Load the neuron positions
fi = h5py.File(network_file, "r")

neuron_id = fi["network/neurons/neuronID"].value
keep_idx = np.array([(x in visualise_id) for x in neuron_id])
neuron_id = neuron_id[keep_idx]

pos = fi["network/neurons/position"][keep_idx, :]
rot = fi["network/neurons/rotation"][keep_idx, :]
morph = fi["network/neurons/morphology"][keep_idx]
name = fi["network/neurons/name"][keep_idx]

origo = fi["meta/simulationOrigo"].value
voxelSize = fi["meta/voxelSize"].value

synapses = fi["network/synapses"]

# Remove the start cube
bpy.ops.object.delete()

bpy.data.scenes['Scene'].render.engine = 'CYCLES'
world = bpy.data.worlds['World']
world.use_nodes = True

# changing these values does affect the render.
bg = world.node_tree.nodes['Background']
if False:  # Set to True for white background
    bg.inputs[0].default_value[:3] = (1.0, 1.0, 1.0)
    bg.inputs[1].default_value = 1.0
else:
    bg.inputs[0].default_value[:3] = (0.0, 0.0, 0.0)
    bg.inputs[1].default_value = 0.0

mat_dspn = bpy.data.materials.new("PKHG")
# matMSD1.diffuse_color = (1.0,0.2,0.2)
mat_dspn.diffuse_color = (77. / 255, 151. / 255, 1.0)

mat_ispn = bpy.data.materials.new("PKHG")
# matMSD2.diffuse_color = (0.2,0.2,1.0)
mat_ispn.diffuse_color = (67. / 255, 55. / 255, 181. / 255)

mat_fs = bpy.data.materials.new("PKHG")
mat_fs.diffuse_color = (6. / 255, 31. / 255, 85. / 255)
mat_chin = bpy.data.materials.new("PKHG")
mat_chin.diffuse_color = (252. / 255, 102. / 255, 0.0)
mat_lts = bpy.data.materials.new("PKHG")
mat_lts.diffuse_color = (150. / 255, 63. / 255, 212. / 255)

mat_other = bpy.data.materials.new("PKHG")
mat_other.diffuse_color = (0.4, 0.4, 0.4)

mat_synapse = bpy.data.materials.new("PKHG")
# matSynapse.diffuse_color = (0.3,1.0,0.3)
mat_synapse.diffuse_color = (1.0, 1.0, 0.9)
# matSynapse.use_transparency = True
mat_synapse.use_nodes = True
# import pdb
# pdb.set_trace()

# Make synapses glow
emission = mat_synapse.node_tree.nodes.new('ShaderNodeEmission')
emission.inputs['Strength'].default_value = 5.0

material_output = mat_synapse.node_tree.nodes.get('Material Output')
mat_synapse.node_tree.links.new(material_output.inputs[0], emission.outputs[0])

for ps, rt, mo, nm, n_id in zip(pos, rot, morph, name, neuron_id):

    e_rot = mathutils.Matrix(rt.reshape(3, 3)).to_euler()
    bpy.ops.import_mesh.swc(filepath=snudda_parse_path(mo, snudda_data))
    obj = bpy.context.selected_objects[0]

    obj.rotation_euler = e_rot
    print(f"Setting position: {ps * 1e3}")
    obj.location = ps * 1e3

    n_type = nm.decode().split("_")[0]
    if n_type == "dSPN" or n_type == "MSD1":
        print(f"{n_id} dSPN")
        mat = mat_dspn
    elif n_type == "iSPN" or n_type == "MSD2":
        print(f"{n_id} iSPN")
        mat = mat_ispn
    elif n_type == "FS":
        print(f"{n_id} FS")
        mat = mat_fs
    elif n_type == "ChIN":
        print(f"{n_id} ChIN")
        mat = mat_chin
    elif n_type == "LTS":
        print(f"{n_id} LTS")
        mat = mat_lts
    else:
        print(f"{n_id} other")
        mat = mat_other

    for ch in obj.children:
        ch.active_material = mat

    obj.select = False

# Draw the synapses

nSynapses = 0

for ob in bpy.context.selected_objects:
    ob.select = False

for syn in synapses:
    pre_id = syn[0]
    post_id = syn[1]

    if pre_id in visualise_id and post_id in visualise_id:
        # Draw this neuron (the SWC import scales from micrometers to mm), the
        # positions in the simulation are in meters, need to scale it to mm for
        # blender to have same units.
        x = (origo[0] + voxelSize * syn[2]) * 1e3
        y = (origo[1] + voxelSize * syn[3]) * 1e3
        z = (origo[2] + voxelSize * syn[4]) * 1e3

        # bpy.ops.mesh.primitive_cube_add(location=(x,y,z))
        # bpy.context.scene.objects.active.dimensions = \
        #  (voxelSize, voxelSize, voxelSize)

        bpy.ops.mesh.primitive_uv_sphere_add(location=(x, y, z), size=0.001)
        ob = bpy.context.selected_objects[0]
        ob.active_material = mat_synapse
        ob.select = False

        nSynapses += 1

        print(f"Added synapse at {[x, y, z]}")

print(f"nSynapses = {nSynapses}")

# import pdb
# pdb.set_trace()

# Add a light source

lamp_data = bpy.data.lamps.new(name="Sun", type='SUN')
lamp_object = bpy.data.objects.new(name="Sun", object_data=lamp_data)
bpy.context.scene.objects.link(lamp_object)

# Place lamp to a specified location
lamp_object.location = (1000.0, 1000.0, 1000.0)

# Reposition camera
cam = bpy.data.objects["Camera"]
cam.location = (2.98, 2.68, 4.96)
cam.rotation_euler = (1.59, 0, -0.26)

if False:
    bpy.ops.render.render()
    bpy.data.images['Render Result'].save_render(filepath="cube-of-neurons-with-neurites.png")

bpy.ops.wm.save_as_mainfile(filepath=out_file)
