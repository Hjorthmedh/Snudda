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

    def bend_morphology(self, morphology: NeuronMorphologyExtended, inside_flag=None):

        if inside_flag is None:
            inside_flag = self.check_if_inside(morphology=morphology)

        # Iterate over all parts of neuron
        for section in morphology.section_iterator():
            # We can use the inside_flag to determine if we need to rotate the branch anything
            pass