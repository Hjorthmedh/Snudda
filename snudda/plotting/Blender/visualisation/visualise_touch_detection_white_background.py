# This script visualises two neurons and the location of the synapses
# connecting them.

import bpy
import mathutils
import numpy as np
import h5py
from snudda.load import SnuddaLoad
from snudda.utils.snudda_path import snudda_parse_path

network_file = "/home/hjorth/HBP/Snudda/snudda/examples/networks/Net10062/networ-synapses.hdf5"
output_file = "twoNeurons-synapses-white-2020.blend"

visualise_id = [443, 7]

# Load the neuron positions
sl = SnuddaLoad(network_file)

neuron_id = sl.data["neuronID"]
keep_idx = np.array([(x in visualise_id) for x in neuron_id])
neuron_id = neuron_id[keep_idx]

neurons = [sl.data["neurons"][x] for x in keep_idx]
origo = sl.data["simulationOrigo"]
voxel_size = sl.data["voxelSize"]


# Remove the start cube
bpy.ops.object.delete()

bpy.data.scenes['Scene'].render.engine = 'CYCLES'
world = bpy.data.worlds['World']
world.use_nodes = True

white_background = True

# changing these values does affect the render.
bg = world.node_tree.nodes['Background']
if white_background:  # Set to True for white background
    bg.inputs[0].default_value[:3] = (1.0, 1.0, 1.0)
    bg.inputs[1].default_value = 1.0
else:
    bg.inputs[0].default_value[:3] = (0.0, 0.0, 0.0)
    bg.inputs[1].default_value = 0.0

mat_msd1 = bpy.data.materials.new("PKHG")
mat_msd1.diffuse_color = (77. / 255, 151. / 255, 1.0)

mat_msd2 = bpy.data.materials.new("PKHG")
mat_msd2.diffuse_color = (67. / 255, 55. / 255, 181. / 255)

mat_fs = bpy.data.materials.new("PKHG")
mat_fs.diffuse_color = (6. / 255, 31. / 255, 85. / 255)
mat_chin = bpy.data.materials.new("PKHG")
mat_chin.diffuse_color = (252. / 255, 102. / 255, 0.0)
mat_lts = bpy.data.materials.new("PKHG")
mat_lts.diffuse_color = (150. / 255, 63. / 255, 212. / 255)

mat_other = bpy.data.materials.new("PKHG")
mat_other.diffuse_color = (0.4, 0.4, 0.4)

mat_synapse = bpy.data.materials.new("PKHG")
if white_background:
    mat_synapse.diffuse_color = (0.8, 0.0, 0.0)
else:
    mat_synapse.diffuse_color = (1.0, 1.0, 0.9)

# matSynapse.use_transparency = True
mat_synapse.use_nodes = True


if not white_background:
    emissionStrength = 5.0

    # Make synapses glow
    emission = mat_synapse.node_tree.nodes.new('ShaderNodeEmission')
    emission.inputs['Strength'].default_value = emissionStrength

    material_output = mat_synapse.node_tree.nodes.get('Material Output')
    mat_synapse.node_tree.links.new(material_output.inputs[0], emission.outputs[0])

material_lookup = { "dspn": mat_msd1,
                    "ispn": mat_msd2,
                    "fsn:" mat_fs,
                    "fs": mat_fs,
                    "chin": mat_chin,
                    "lts": mat_lts,
                    "other": mat_other}

for neuron in neurons:

    e_rot = mathutils.Matrix(neuron["rotation"].reshape(3, 3)).to_euler()
    bpy.ops.import_mesh.swc(filepath=snudda_parse_path(mo))
    obj = bpy.context.selected_objects[0]

    obj.rotation_euler = e_rot
    print(f"Setting position: {neuron['position'] * 1e3}")
    obj.location = neuron["position"] * 1e3

    n_type = neuron["type"].lower()
    if n_type in material_lookup:
      mat = material_lookup[n_type]
    else:
      mat = material_lookup["other"]

    for ch in obj.children:
        ch.active_material = mat

    obj.select = False

# Draw the synapses

n_synapses = 0

for ob in bpy.context.selected_objects:
    ob.select = False

for vis_pre_id in visualise_id:
    for vis_post_id in visualise_id:
        synapses = sl.find_synapses(pre_id=vis_pre_id, post_id=vis_post_id)

        for syn in synapses:
            pre_id = syn[0]
            post_id = syn[1]

            assert pre_id == vis_pre_id and post_id == vis_post_id  # Just sanity check, should be true

            # Draw this neuron (the SWC import scales from micrometers to mm), the
            # positions in the simulation are in meters, need to scale it to mm for
            # blender to have same units.
            x = (origo[0] + voxel_size * syn[2]) * 1e3
            y = (origo[1] + voxel_size * syn[3]) * 1e3
            z = (origo[2] + voxel_size * syn[4]) * 1e3

            bpy.ops.mesh.primitive_uv_sphere_add(location=(x, y, z), size=0.001 * 4)
            ob = bpy.context.selected_objects[0]
            ob.active_material = mat_synapse
            ob.select = False

            n_synapses += 1

            print(f"Added synapse at {[x, y, z]}")

print(f"nSynapses = {n_synapses}")

# Add a light source

lamp_data = bpy.data.lamps.new(name="Sun", type='SUN')
lamp_object = bpy.data.objects.new(name="Sun", object_data=lamp_data)
bpy.context.scene.objects.link(lamp_object)

# Place lamp to a specified location
lamp_object.location = (1000.0, 1000.0, 1000.0)

# Reposition camera
cam = bpy.data.objects["Camera"]
# cam.location = (2.98,2.68,4.96)
# cam.rotation_euler = (1.59,0,-0.26)

cam.location = (3.19, 3.46, 4.96)
cam.rotation_euler = (96.7 * np.pi / 180, 0, -14.3 * np.pi / 180)

if False:
    bpy.ops.render.render()
    bpy.data.images['Render Result'].save_render(filepath="cube-of-neurons-with-neurites.png")

bpy.ops.wm.save_as_mainfile(filepath=output_file)
