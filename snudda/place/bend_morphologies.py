import numpy as np

from snudda.neurons import NeuronMorphologyExtended
from snudda.place.region_mesh_redux import RegionMeshRedux


class BendMorphologies:

    def __init__(self, region_mesh: RegionMeshRedux):

        self.region_mesh = region_mesh

        pass

    def check_if_inside(self, morphology: NeuronMorphologyExtended):

        coords = morphology.morphology_data["neuron"].geometry
        inside_flag = self.region_mesh.check_inside(points=coords[:, :3])

        return inside_flag

    def bend_morphology(self, morphology: NeuronMorphologyExtended):

        # Iterate over all parts of neuron
        for section in morphology.section_iterator():
            # We need to track the rotation of each point, in particular save rotations
            # for each branch point, so that all children can start with that rotation

            # Loop over all points in section

            # Gradient is calculated based on distance to mesh
            # Calculate delta_gradient = gradient_self - gradient_parent
            # Amount of angle to bend f(delta_gradient)
            # Randomize the rotation matrix

            # After section is done, store the last rotation matrix, so its children can get their parent rotation
            # store in dictionary?




            pass