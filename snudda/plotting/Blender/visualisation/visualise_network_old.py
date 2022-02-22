# This script visualises two neurons and the location of the synapses
# connecting them.

import bpy
import os
import mathutils
import numpy as np
from snudda.utils.load import SnuddaLoad
from snudda.utils.snudda_path import snudda_parse_path


class VisualiseNetworkOld(object):

    # You need to provide neuron
    def __init__(self, network_path, blender_save_file=None, blender_output_image=None,
                 network_json=None):

        self.network_path = network_path

        if network_json:
            self.network_json = network_json
            self.network_file = None
        else:
            self.network_json = None
            self.network_file = os.path.join(network_path, "network-synapses.hdf5")

        if blender_save_file:
            self.blender_save_file = blender_save_file
        else:
            self.blender_save_file = os.path.join(network_path, "visualise-network.blend")

        self.blender_output_image = blender_output_image

        self.neuron_cache = dict([])

        # Load the neuron positions
        if self.network_file:
            self.sl = SnuddaLoad(self.network_file)
            self.data = self.sl.data
        elif self.network_json:
            from snudda.utils.fake_load import FakeLoad
            self.sl = FakeLoad()
            self.sl.import_json(self.network_json)
            self.data = self.sl.data

    def visualise(self, neuron_id=None, blender_output_image=None, white_background=True,
                  show_synapses=True):

        if neuron_id:
            neurons = [self.data["neurons"][x] for x in neuron_id]
        else:
            neurons = self.data["neurons"]
            neuron_id = self.data["neuronID"]

        if blender_output_image:
            self.blender_output_image = blender_output_image

        origo = self.data["simulationOrigo"]
        voxel_size = self.data["voxelSize"]

        # Remove the start cube
        # bpy.ops.object.delete()
        VisualiseNetworkOld.clean_scene()

        bpy.data.scenes['Scene'].render.engine = 'CYCLES'
        world = bpy.data.worlds['World']
        world.use_nodes = True

        # changing these values does affect the render.
        bg = world.node_tree.nodes['Background']
        if white_background:  # Set to True for white background
            bg.inputs[0].default_value[:3] = (1.0, 1.0, 1.0)
            bg.inputs[1].default_value = 1.0
        else:
            bg.inputs[0].default_value[:3] = (0.0, 0.0, 0.0)
            bg.inputs[1].default_value = 0.0

        # Define materials
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

        material_lookup = { "dspn": mat_msd1,
                            "ispn": mat_msd2,
                            "fsn": mat_fs,
                            "fs": mat_fs,
                            "chin": mat_chin,
                            "lts": mat_lts,
                            "synapse": mat_synapse,
                            "other": mat_other}

        if white_background:
            mat_synapse.diffuse_color = (0.8, 0.0, 0.0)
        else:
            mat_synapse.diffuse_color = (1.0, 1.0, 0.9)

        # matSynapse.use_transparency = True
        mat_synapse.use_nodes = True

        if not white_background:
            emission_strength = 5.0

            # Make synapses glow
            emission = mat_synapse.node_tree.nodes.new('ShaderNodeEmission')
            emission.inputs['Strength'].default_value = emission_strength

            material_output = mat_synapse.node_tree.nodes.get('Material Output')
            mat_synapse.node_tree.links.new(material_output.inputs[0], emission.outputs[0])

        for neuron in neurons:

            e_rot = mathutils.Matrix(neuron["rotation"].reshape(3, 3)).to_euler()

            if neuron["name"] in self.neuron_cache:
                # If we already have the object in memory, copy it.
                obj = self.neuron_cache[neuron["name"]].copy()

                if self.neuron_cache[neuron["name"]].data:
                    obj.data = self.neuron_cache[neuron["name"]].data.copy()

                VisualiseNetworkOld.copy_children(self.neuron_cache[neuron["name"]], obj)
                obj.animation_data_clear()

                # Will return None if there is no obj named CUBe
                obj.name = f"{neuron['name']}-{neuron['neuronID']}"
                VisualiseNetworkOld.link_object(obj)
            else:
                bpy.ops.import_mesh.swc(filepath=snudda_parse_path(neuron["morphology"]))
                obj = bpy.context.selected_objects[0]
                obj.name = f"{neuron['name']}-{neuron['neuronID']}"

                self.neuron_cache[neuron["name"]] = obj

            obj.rotation_euler = e_rot
            print(f"Setting neuron {neuron['neuronID']} ({neuron['name']}) position: {neuron['position'] * 1e3}")
            obj.location = neuron["position"] * 1e3

            n_type = neuron["type"].lower()
            if n_type in material_lookup:
                mat = material_lookup[n_type]
            else:
                mat = material_lookup["other"]

            for ch in obj.children:
                ch.active_material = mat

            obj.select = False

        if show_synapses:
            print("Adding synapses...")

            # Draw the synapses
            n_synapses = 0

            for ob in bpy.context.selected_objects:
                ob.select = False

            synapse_obj = None

            for vis_pre_id in neuron_id:
                for vis_post_id in neuron_id:
                    synapses, synapse_coords = self.sl.find_synapses(pre_id=vis_pre_id, post_id=vis_post_id)

                    if synapses is None:
                        # No synapses between pair
                        continue

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

                        if synapse_obj:
                            obj = synapse_obj.copy()
                            if synapse_obj.data:
                                obj.data = synapse_obj.data.copy()
                            obj.animation_data_clear()
                            obj.location = (x, y, z)
                            obj.name = f"synapse-{n_synapses}"
                            VisualiseNetworkOld.link_object(obj)

                        else:
                            bpy.ops.mesh.primitive_uv_sphere_add(location=(x, y, z), size=0.001 * 4)
                            obj = bpy.context.selected_objects[0]
                            obj.active_material = mat_synapse
                            obj.select = False
                            synapse_obj = obj

                        n_synapses += 1

                        # print(f"Added synapse #{n_synapses} at {[x, y, z]}")
                        if n_synapses % 5000 == 0:
                            print(f"Synapses added so far: {n_synapses}")

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

        # Is this needed?
        bpy.context.scene.update()

        bpy.ops.wm.save_as_mainfile(filepath=self.blender_save_file)

        if self.blender_output_image:
            print("Rendering image.")
            bpy.ops.render.render()
            bpy.data.images['Render Result'].save_render(filepath=self.blender_output_image)

    @staticmethod
    def copy_children(parent, parent_copy):
        for child in parent.children:
            child_copy = child.copy()
            child_copy.parent = parent_copy
            VisualiseNetworkOld.link_object(child_copy)
            VisualiseNetworkOld.copy_children(child, child_copy)

    @staticmethod
    def link_object(obj):
        try:
            bpy.context.scene.objects.link(obj)  # Blender 2.7
        except:
            print("Blender 2.7 failed, switch over to only use Blender 2.8 syntax")
            bpy.context.collection.objects.link(obj)  # Blender 2.8

    @staticmethod
    def clean_scene():
        # TODO: This does not seem to remove everything. Had some leftover synapses present still in notebook.
        print("Cleaning the scene.")
        del_list = bpy.context.copy()
        del_list['selected_objects'] = list(bpy.context.scene.objects)
        bpy.ops.object.delete(del_list)