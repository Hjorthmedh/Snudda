# Visualises neurons and the location of the synapses connecting them.
# Tested using Blender 2.93. Should not be used with older version of Blender.
#from snudda.plotting.Blender.visualisation.visualise_network import VisualiseNetwork2
from snudda.plotting.Blender.visualisation.visualise_network2 import VisualiseNetwork2
import os
import bpy

####################Adjust to your own network location###########
network_path = os.path.join("networks","example100")
##################################################################

blender_output_image = os.path.join(network_path, "neuron-rendering.png")
vn = VisualiseNetwork2(network_path=network_path, blender_output_image=blender_output_image)
vn.visualise(neuron_id=[0,1,2,7])#Visualize selected neurons
#vn.visualise() #Visualize all neurons in your network. Might take a while (Around 0.75h for 100neurons on laptop)

