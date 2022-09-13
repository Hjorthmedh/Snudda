# Visualises neurons and the location of the synapses connecting them.
# Tested using Blender 2.93. Should not be used with older version of Blender.

import bpy
from bpy_extras.io_utils import unpack_list

import os
import mathutils
import numpy as np
from snudda.utils.load import SnuddaLoad
from snudda.utils.snudda_path import snudda_parse_path, get_snudda_data
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation


class VisualiseNetwork(object):

    # You need to provide neuron
    def __init__(self, network_path, blender_save_file=None, blender_output_image=None,
                 network_json=None, simulation_output_file_name=None):

        self.network_path = network_path
        self.snudda_data = get_snudda_data(network_path=network_path)
        self.scale_f = 1000  # factor to downscale the data

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

        if simulation_output_file_name:
            self.slns = SnuddaLoadNetworkSimulation(simulation_output_file_name)
            self.spike_times = self.slns.get_spikes()
        else:
            self.spike_times = None

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

    def visualise(self,
                  neuron_id=None,
                  blender_output_image=None,
                  white_background=True,
                  show_synapses=True,
                  draw_meshes=True,
                  full_meshes=None,
                  camera_location=None,
                  camera_rotation=None,
                  camera_scale=None,
                  detail_level=1):

        """
            Visualise network in blender.

            Args:
                neuron_id (list) : Neuron ID to visualise, default None means all
                blender_output_image
                white_background
                show_synapses
                camera_location
                camera_rotation
                camera_scale
                detail_level (int) : 1 = full morphologies, 2 = reduced morphologies, 3 = soma only

        """

        if neuron_id:
            neurons = [self.data["neurons"][x] for x in neuron_id]
        else:
            neurons = self.data["neurons"]
            neuron_id = self.data["neuronID"]

        if blender_output_image:
            self.blender_output_image = blender_output_image

        if camera_location is None:
            camera_location = (2.98, 2.68, 4.96)

        if camera_rotation is None:
            camera_rotation = (1.59, 0, -0.26)

        if camera_scale is None:
            camera_scale = (1, 1, 1)

        origo = self.data["simulationOrigo"]
        voxel_size = self.data["voxelSize"]

        # Remove the start cube
        VisualiseNetwork.clean_scene()

        # Add a light source
        # TODO: Add choice of adding a light source, instead of True here.
        if True:
            sun_location = (10, 10, 10)
            lamp_data = bpy.data.lights.new(name="Sun", type='SUN')

        # Create new object, pass the light data
            sun_object = bpy.data.objects.new(name="sun_object", object_data=lamp_data)

        # Link object to collection in context
            bpy.context.collection.objects.link(sun_object)

        # Change light position
            sun_object.location = sun_location

        bpy.data.scenes['Scene'].render.engine = 'CYCLES'
        world = bpy.data.worlds['World']
        world.use_nodes = True

        # changing these values does affect the render.
        bg = world.node_tree.nodes['Background']
        if white_background:  # Set to True for white background
            bg.inputs[0].default_value[:3] = (1.0, 1.0, 1.0)
            bg.inputs[1].default_value = 1.0
            bpy.context.scene.view_settings.view_transform = 'Standard'

        else:
            bg.inputs[0].default_value[:3] = (0.0, 0.0, 0.0)
            bg.inputs[1].default_value = 0.0

        # Define materials
        mat_dspn = bpy.data.materials.new("PKHG")
        mat_dspn.diffuse_color = (77. / 255, 151. / 255, 1.0, 0.5)
        mat_ispn = bpy.data.materials.new("PKHG")
        mat_ispn.diffuse_color = (67. / 255, 55. / 255, 181. / 255, 0.5)
        mat_fs = bpy.data.materials.new("PKHG")
        mat_fs.diffuse_color = (6. / 255, 31. / 255, 85. / 255, 1.0)
        mat_chin = bpy.data.materials.new("PKHG")
        mat_chin.diffuse_color = (252. / 255, 102. / 255, 0.0, 1.0)
        mat_lts = bpy.data.materials.new("PKHG")
        mat_lts.diffuse_color = (150. / 255, 63. / 255, 212. / 255, 1.0)

        mat_snr = bpy.data.materials.new("PKHG")

        # PURPLE

        mat_snr.diffuse_color = (102 / 255, 0 / 255, 102 / 255, 1.0)
        mat_snr.diffuse_color = (102 / 255, 0 / 255, 102 / 255, 1.0)

        # BLUE-PURPLE

        mat_arky = bpy.data.materials.new("PKHG")
        mat_arky.diffuse_color = (43 / 255, 0 / 255, 255 / 255, 1.0)

        # GREEN
        mat_proto = bpy.data.materials.new("PKHG")
        mat_proto.diffuse_color = (0 / 255, 130 / 255, 0 / 255, 1.0)

        mat_other = bpy.data.materials.new("PKHG")
        mat_other.diffuse_color = (0.4, 0.4, 0.4, 1.0)
        mat_synapse = bpy.data.materials.new("PKHG")

        material_lookup = {"dspn": mat_dspn,
                           "ispn": mat_ispn,
                           "fsn": mat_fs,
                           "fs": mat_fs,
                           "chin": mat_chin,
                           "lts": mat_lts,
                           "snrneurons": mat_snr,
                           "proto": mat_proto,
                           "arky": mat_arky,
                           "synapse": mat_synapse,
                           "other": mat_other}

        if white_background:
            mat_synapse.diffuse_color = (0.8, 0.0, 0.0, 1.0)
        else:
            mat_synapse.diffuse_color = (1.0, 1.0, 0.9, 1.0)

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

                VisualiseNetwork.copy_children(self.neuron_cache[neuron["name"]], obj)
                obj.animation_data_clear()

                # Will return None if there is no obj named CUBe
                obj.name = f"{neuron['name']}-{neuron['neuronID']}"
                VisualiseNetwork.link_object(obj)
            else:
                self.read_swc_data(filepath=snudda_parse_path(neuron["morphology"], self.snudda_data), detail_level=detail_level)
                obj = bpy.context.selected_objects[0]
                obj.name = f"{neuron['name']}-{neuron['neuronID']}"

                self.neuron_cache[neuron["name"]] = obj

            obj.rotation_euler = e_rot

            print(f"Setting neuron {neuron['neuronID']} ({neuron['name']}) position: {neuron['position']}")
            obj.location = neuron["position"] * self.scale_f

            n_type = neuron["type"].lower()

            if n_type in material_lookup:
                mat = material_lookup[n_type]
            else:
                mat = material_lookup["other"]

            if self.spike_times:
                rest_color = mat.diffuse_color[:]
                # if animating spike times we need to make a fresh material per neuron
                mat_spikes = bpy.data.materials.new("PKHG")
                mat_spikes.diffuse_color = rest_color
                mat_spikes.keyframe_insert(data_path="diffuse_color", frame=1.0, index=-1)
                if str(neuron['neuronID']) in self.spike_times.keys():
                    spikes = self.spike_times[str(neuron['neuronID'])]

                    # convert 'time' to Blender frames; factor of ~100 works nicely
                    spike_frames = np.round(100 * np.array(spikes))
                    for t in spike_frames:
                        mat_spikes.diffuse_color = rest_color
                        # need to add an instruction to remain at rest colour at a pre-spike time
                        # so that the colour change is not gradual but quasi-instantaneous
                        mat_spikes.keyframe_insert(data_path="diffuse_color", frame=t - 1, index=-1)
                        mat_spikes.diffuse_color = (1, 1, 1, 1)
                        mat_spikes.keyframe_insert(data_path="diffuse_color", frame=t, index=-1)
                        mat_spikes.diffuse_color = rest_color
                        # change back to rest colour
                        mat_spikes.keyframe_insert(data_path="diffuse_color", frame=t + 5, index=-1)
                print("Color......")
                for ch in obj.children:
                    ch.active_material = mat_spikes
            else:
                print("Color......")
                for ch in obj.children:
                    ch.active_material = mat

            obj.select_set(False)

        if show_synapses:
            print("Adding synapses...")
            self.syn_coll = bpy.data.collections.new(name="Synapses") 
            bpy.context.scene.collection.children.link(self.syn_coll)
            # Draw the synapses
            n_synapses = 0

            for ob in bpy.context.selected_objects:
                ob.select_set(False)
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
                        x = (origo[0] + voxel_size * syn[2]) * self.scale_f
                        y = (origo[1] + voxel_size * syn[3]) * self.scale_f
                        z = (origo[2] + voxel_size * syn[4]) * self.scale_f

                        if synapse_obj:
                            obj = synapse_obj.copy()
                            if synapse_obj.data:
                                obj.data = synapse_obj.data.copy()
                            obj.animation_data_clear()
                            obj.location = (x, y, z)
                            obj.name = f"synapse-{n_synapses}"
                            VisualiseNetwork.link_object(obj)

                        else:
                            bpy.ops.mesh.primitive_uv_sphere_add(radius=0.001 * 4, location=(x, y, z), scale=(1, 1, 1),
                                                                 segments=16, ring_count=8)
                            obj = bpy.context.selected_objects[0]
                            obj.active_material = mat_synapse
                            obj.select_set(False)
                            synapse_obj = obj

                        n_synapses += 1

                        # print(f"Added synapse #{n_synapses} at {[x, y, z]}")
                        if n_synapses % 5000 == 0:
                            print(f"Synapses added so far: {n_synapses}")
            bpy.ops.object.select_all( action='DESELECT' )
            bpy.ops.object.select_pattern(pattern="synapse*")
            s_objs = bpy.context.selected_objects
            for s in s_objs:
                self.syn_coll.objects.link(s)
            print(f"nSynapses = {n_synapses}")

        if draw_meshes:
            self.struct_coll = bpy.data.collections.new(name="Structures") #create new coll in data
            bpy.context.scene.collection.children.link(self.struct_coll) #add new coll to the scene
            self.add_all_meshes()

        if full_meshes:
            for struct, mesh_file in full_meshes.items():
                self.add_mesh_structure(mesh_file=snudda_parse_path(mesh_file, self.snudda_data), colour=(0.1, 0.1, 0.1),
                                    alpha=0.1)

        bpy.ops.object.camera_add(enter_editmode=False, align='VIEW',
                                  location=camera_location,
                                  rotation=camera_rotation,
                                  scale=camera_scale)
        cam = bpy.data.objects["Camera"]
        print(bpy.context.selected_objects)
        bpy.context.scene.camera = cam

        bpy.ops.wm.save_as_mainfile(filepath=self.blender_save_file)

        if self.blender_output_image:
            # When testing (e.g. adjusting camera position) you can skip this and just look at the blender file directly
            print("Rendering image.")
            bpy.ops.render.render()
            bpy.data.images['Render Result'].save_render(filepath=self.blender_output_image)

    def add_mesh_structure(self, mesh_file, colour, alpha):

        mat = bpy.data.materials.new("PKHG")
        mat.diffuse_color = (colour[0], colour[1], colour[2], alpha)
        mat.use_nodes = True
        mat.node_tree.nodes["Principled BSDF"].inputs['Alpha'].default_value = alpha
        mat.node_tree.nodes["Principled BSDF"].inputs['Base Color'].default_value = (colour[0], colour[1], colour[2], alpha)

        structure_object = bpy.ops.import_scene.obj(filepath=mesh_file, axis_up="Z", axis_forward="Y")
        o = bpy.context.selected_objects[0]
        # scale_f = 1000
        o.scale[0] = 1 / self.scale_f
        o.scale[1] = 1 / self.scale_f
        o.scale[2] = 1 / self.scale_f
        o.active_material = mat
        self.struct_coll.objects.link(o)
    def add_all_meshes(self):

        for name, structure in self.sl.config["Volume"].items():
            self.add_mesh_structure(mesh_file=snudda_parse_path(structure["meshFile"], self.snudda_data), colour=(0.1, 0.1, 0.1),
                                    alpha=0.1)

    @staticmethod
    def copy_children(parent, parent_copy):
        for child in parent.children:
            child_copy = child.copy()
            child_copy.parent = parent_copy
            VisualiseNetwork.link_object(child_copy)
            VisualiseNetwork.copy_children(child, child_copy)

    @staticmethod
    def link_object(obj):
        try:
            bpy.context.collection.objects.link(obj)
        except:
            print("Blender 2.8 failed. Likely due to 2.7 syntax.")
            # bpy.context.scene.objects.link(obj)  # Blender 2.7

    @staticmethod
    def clean_scene():
        # TODO: This does not seem to remove everything. Had some leftover synapses present still in notebook.
        print("Cleaning the scene.")
        del_list = bpy.context.copy()
        del_list['selected_objects'] = list(bpy.context.scene.objects)
        bpy.ops.object.delete(del_list)

    def read_swc_data(self, filepath, detail_level=1):

        """
            Read SWC file

            Args:
                filepath (str) : Path to SWC file
                detail_level (int) : Detail level 1 = full detail, 2 = reduced quality, 3 = soma only

        """

        ''' read swc file '''
        print(filepath)
        f = open(filepath)
        lines = f.readlines()
        f.close()

        ''' find starting point '''
        x = 0
        while lines[x][0] == '#':
            x += 1

        ''' Create a dictionary with the first item '''
        data = lines[x].strip().split(' ')
        soma_id = int(data[0])
        soma_type = float(data[1])
        soma_x = float(data[2])
        soma_y = float(data[3])
        soma_z = float(data[4])
        soma_r = float(data[5])
        soma_parent = int(data[6])

        # We centre the neuron
        neuron = {soma_id: [soma_type, 0.0, 0.0, 0.0, soma_r, soma_parent]}
        x += 1

        ''' Read the rest of the lines to the dictionary '''
        for ls in lines[x:]:
            data = ls.strip().split(' ')
            comp_id = int(data[0])
            comp_type = float(data[1])
            comp_x = float(data[2])
            comp_y = float(data[3])
            comp_z = float(data[4])
            comp_r = float(data[5])
            comp_parent = int(data[6])

            # Centre neuron, so soma is at 0,0,0
            neuron[comp_id] = [comp_type, comp_x - soma_x, comp_y - soma_y, comp_z - soma_z, comp_r, comp_parent]

        bpy.ops.object.empty_add(type='ARROWS',
                                 location=(
                                     neuron[1][1] / self.scale_f, neuron[1][2] / self.scale_f,
                                     neuron[1][3] / self.scale_f),
                                 rotation=(0, 0, 0))
        a = bpy.context.selected_objects[0]
        a.name = 'neuron_swc'

        last = -10.0

        line_points = []
        line_radius = []

        ''' Create object '''
        for key, value in neuron.items():
            # value contains: 0: type, 1: x, 2: y, 3: z, 4: r, 5: parent

            if value[0] == 1:

                if detail_level > 1:
                    segments = 10
                    ring_count = 5
                else:
                    segments = 32
                    ring_count = 16

                # This is the soma, add it
                soma_radie = value[-2]
                bpy.ops.mesh.primitive_uv_sphere_add(segments=segments,
                                                     ring_count=ring_count,
                                                     location=(
                                                         value[1] / self.scale_f, value[2] / self.scale_f,
                                                         value[3] / self.scale_f),
                                                     radius=soma_radie / self.scale_f)
                soma_obj = bpy.context.selected_objects[0]
                soma_obj.parent = a

                print(f"Adding soma {value}")

            if value[-1] == -1:
                continue

            if value[0] == 10:
                continue

            if value[-1] == last:
                line_points.append(
                    [neuron[value[-1]][1] / self.scale_f, neuron[value[-1]][2] / self.scale_f,
                     neuron[value[-1]][3] / self.scale_f])
                line_radius.append(neuron[value[-1]][4] / self.scale_f)
            else:
                # Add Bezier curve for previous data
                VisualiseNetwork.add_bezier(curve_parent=a,
                                            line_points=line_points,
                                            line_radius=line_radius,
                                            detail_level=detail_level)
                line_points = []
                line_radius = []

            last = key

        # Add the last line
        VisualiseNetwork.add_bezier(curve_parent=a,
                                    line_points=line_points,
                                    line_radius=line_radius,
                                    detail_level=detail_level)
        line_points = []
        line_radius = []

        a.select_set(True)

        return {'FINISHED'}

    @staticmethod
    def add_bezier(curve_parent, line_points, line_radius, detail_level=1):

        if len(line_points) == 0:
            return

        if detail_level > 2:
            # Only soma
            return

        if detail_level == 2:
            # Keep only end points and middle point
            if len(line_radius) > 3:
                line_points = [line_points[0],
                               line_points[int(len(line_radius) / 2)],
                               line_points[-1]]
                line_radius = [line_radius[0],
                               line_radius[int(len(line_radius) / 2)],
                               line_radius[-1]]

        tracer = bpy.data.curves.new('tracer', 'CURVE')
        tracer.dimensions = '3D'
        spline = tracer.splines.new('BEZIER')
        spline.bezier_points.add(len(line_points) - 1)

        curve = bpy.data.objects.new('curve', tracer)
        curve.data.use_fill_caps = True  # Added 2019-06-17
        curve.parent = curve_parent

        bpy.context.scene.collection.objects.link(curve)

        # render ready curve
        tracer.resolution_u = 8
        tracer.bevel_resolution = 8
        tracer.fill_mode = 'FULL'
        tracer.bevel_depth = 1.0

        # move nodes to objects
        # spline.bezier_points.foreach_set("co", unpack_list(line_points))
        # spline.bezier_points.foreach_set("radius", unpack_list(line_radius))
        # spline.bezier_points.foreach_set("handle_right_type", unpack_list(["VECTOR"] * len(line_radius)))
        # spline.bezier_points.foreach_set("handle_left_type", unpack_list(["VECTOR"] * len(line_radius)))

        for p, r, co in zip(spline.bezier_points, line_radius, line_points):
            p.radius = r
            p.co = co
            p.handle_right_type = "VECTOR"
            p.handle_left_type = "VECTOR"

# TODO: Look for speedup -- https://blender.stackexchange.com/questions/7358/python-performance-with-blender-operators
# TODO: https://blenderartists.org/t/python-slowing-down-over-time/569534/8
