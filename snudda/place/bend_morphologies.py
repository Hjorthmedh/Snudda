import numpy as np
import warnings
from scipy.spatial.transform import Rotation

from snudda.neurons import NeuronMorphologyExtended
from snudda.neurons.morphology_data import MorphologyData, SectionMetaData
from snudda.place.region_mesh_redux import RegionMeshRedux


class BendMorphologies:

    def __init__(self, region_mesh: RegionMeshRedux, rng):

        if type(region_mesh) == str:
            region_mesh = RegionMeshRedux(mesh_path=region_mesh)

        self.region_mesh = region_mesh
        self.rng = rng

    def check_if_inside(self, morphology: NeuronMorphologyExtended):

        coords = morphology.morphology_data["neuron"].geometry
        inside_flag = self.region_mesh.check_inside(points=coords[:, :3])

        return inside_flag

    def bend_morphology(self, morphology: NeuronMorphologyExtended,
                        k_dist=30e-6, max_angle=0.1,  # angle in radians
                        n_random=5, random_seed=None):

        # k -- how early will the neuron start bending when it approaches the border

        # print(f"random_seed = {random_seed}")

        if random_seed is not None:
            rng = np.random.default_rng(random_seed)
        else:
            rng = self.rng

        # Check distance to border for all points, negative distance means inside
        all_original_dist = self.region_mesh.distance_to_border(morphology.geometry[:, :3])
        if (all_original_dist < 0).all():
            # Morphology entirely inside mesh, nothing to do
            return None, False

        # We randomize n_random points, and store candidates here, and pick "best" one
        candidate_pos = np.zeros((n_random, 3))

        parent_direction = dict()

        old_rotation_representation = self.get_full_rotation_representation(morphology=morphology)
        new_rotation_representation = dict()
        morphology_changed = False

        for section in morphology.section_iterator():
            if (section.section_id, section.section_type) in parent_direction:
                parent_dir, parent_point, parent_dist, parent_moved = parent_direction[section.section_id, section.section_type]
            else:
                if morphology.rotation is not None:
                    parent_dir = np.matmul(morphology.rotation, np.array([[0], [0], [1]])).T
                else:
                    parent_dir = np.array([[0, 0, 1]])

                if morphology.position is not None:
                    parent_point = morphology.position
                else:
                    parent_point = np.zeros((3, ))

                parent_dist = self.region_mesh.distance_to_border(points=parent_point.reshape((1, 3)))[0]
                parent_moved = False

            rot_rep = old_rotation_representation[section.section_id, section.section_type]
            new_rot_rep = []

            section_dist = all_original_dist[section.point_idx]

            # Loop over all points in section
            for idx, (rotation, length) in enumerate(rot_rep):

                try:
                    segment_direction = rotation.apply(parent_dir)
                    putative_point = segment_direction * length + parent_point
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

                # TODO: This check can be cached (assuming no parent segment have been rotated, track that)

                # Check if point is too close to edge
                if parent_moved:
                    dist = self.region_mesh.distance_to_border(points=putative_point)[0]
                else:
                    # Parent has not moved, so use stored original distance
                    dist = section_dist[idx]

                P_move = 1 / (1 + np.exp(-dist/k_dist))

                # Cache the random numbers for segments in the section...
                if dist > parent_dist and rng.uniform() < P_move:

                    morphology_changed = True
                    parent_moved = True

                    # We need to randomize new rotation matrix
                    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html
                    angles = rng.uniform(size=(n_random, 3), low=-max_angle, high=max_angle)  # Angles in radians
                    avoidance_rotations = Rotation.from_euler(seq="XYZ", angles=angles)

                    for idx2, av_rot in enumerate(avoidance_rotations):
                        candidate_pos[idx2, :] = parent_point + length * (av_rot*rotation).apply(vectors=parent_dir)

                    candidate_dist = self.region_mesh.distance_to_border(points=candidate_pos)

                    # We want the smallest (or most negative) distance
                    picked_idx = np.argsort(candidate_dist)[0]
                    dist = candidate_dist[picked_idx]

                    new_rot = avoidance_rotations[picked_idx]*rotation
                    new_rot_rep.append((new_rot, length))
                    segment_direction = new_rot.apply(parent_dir)

                else:
                    new_rot_rep.append((rotation, length))

                parent_point = segment_direction * length + parent_point
                parent_dir = segment_direction
                parent_dist = dist

            for child_id, child_type in section.child_section_id.T:
                parent_direction[child_id, child_type] = (parent_dir, parent_point, parent_dist, parent_moved)

            new_rotation_representation[section.section_id, section.section_type] = new_rot_rep

        return new_rotation_representation, morphology_changed

    def get_full_rotation_representation(self, morphology: MorphologyData):

        rotation_representation = dict()
        parent_direction = dict()

        for section in morphology.section_iterator():

            if (section.section_id, section.section_type) in parent_direction.keys():
                parent_dir = parent_direction[section.section_id, section.section_type]
            else:
                parent_dir = np.array([[0, 0, 1]])

            try:
                rot_and_len, last_direction = self.rotation_representation(section=section, parent_direction=parent_dir)
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

            rotation_representation[section.section_id, section.section_type] = rot_and_len

            for child_id, child_type in section.child_section_id.T:
                parent_direction[child_id, child_type] = last_direction

        return rotation_representation

    def apply_rotation(self, morphology: MorphologyData, rotation_representation):

        parent_direction = dict()

        # new_coords = np.zeros((morphology.geometry.shape[0], 3))
        new_coords = np.full((morphology.geometry.shape[0], 3), np.nan)  # Use nan as fill value, to see problems

        for section in morphology.section_iterator():
            if (section.section_id, section.section_type) in parent_direction:
                parent_dir, parent_pos = parent_direction[section.section_id, section.section_type]
            else:
                if morphology.rotation is not None:
                    parent_dir = np.matmul(morphology.rotation, np.array([[0, 0, 1]]).T).T
                else:
                    parent_dir = np.array([[0, 0, 1]])

                if morphology.position is not None:
                    parent_pos = morphology.position
                else:
                    parent_pos = np.zeros((3, ))

            if section.section_type == section.parent_section_type or section.section_type == 1:
                # The section includes the parent point, if both sections are of the same type
                # or if the section is the soma (since then it has no rotations, so position gets passed as parent point)
                include_parent_point = True
            else:
                include_parent_point = False

            rot_rep = rotation_representation[section.section_id, section.section_type]
            coords, last_dir = self.coordinate_representation(rotation_representation=rot_rep,
                                                              parent_direction=parent_dir,
                                                              parent_point=parent_pos,
                                                              return_last_direction=True,
                                                              include_parent_point=include_parent_point)

            new_coords[section.point_idx, :3] = coords

            if np.isnan(coords).any():
                raise ValueError(f"NaN coordinates calculated.")

            for child_id, child_type in section.child_section_id.T:
                parent_direction[child_id, child_type] = (last_dir, coords[-1, :3])

        # TODO: This just returns the coords for now, add option to update coords in morphology?
        #       OBS! Then rotation should also be reset, since it is now included in the coordinates

        if np.isnan(new_coords).any():
            print(f"Some coordinates were not assigned!")
            bad_idx = np.where(np.sum(np.isnan(new_coords), axis=1))[0]
            print(f"missing_idx = {bad_idx}")
            import pdb
            pdb.set_trace()

        return new_coords

    def rotation_representation(self, section: SectionMetaData, parent_direction=None):

        """ Represent each section as a series of length of segment, and rotations relative the parent segment."""

        if parent_direction is None:
            # parent_direction = np.array([[1, 0, 0]])
            parent_direction = np.array([0, 0, 1])

        rotations_and_length = []
        parent_direction = parent_direction / np.linalg.norm(parent_direction)
        segment_direction = parent_direction  # To handle soma...

        # If parent compartment is of a different type (e.g. soma parent, for the dendrite) then we need
        # to make sure that the first point is also included. So we need to artificially add the soma.
        # However, the soma itself has no parent, so in that case we should not do it.

        if section.section_type != section.parent_section_type and section.parent_section_type != -1:
            coords = np.zeros((len(section.point_idx) + 1, 3))
            coords[0, :3] = section.morphology_data.geometry[section.parent_point_idx, :3]
            coords[1:, :3] = section.morphology_data.geometry[section.point_idx, :3]
        else:
            coords = section.morphology_data.geometry[section.point_idx, :3]

        delta = np.diff(coords, axis=0)
        delta_length = np.linalg.norm(delta, axis=1)
        delta_direction = delta / delta_length.reshape((delta.shape[0], 1))

        for segment_direction, segment_length in zip(delta_direction, delta_length):
            # segment_direction = segment_direction.reshape((1, 3))
            try:
                if np.allclose(parent_direction, segment_direction) or 1.00001 > np.dot(parent_direction, segment_direction) > 1:
                    # Sometimes numerical precision gives us values larger than 1
                    rotation = Rotation.from_euler("xyz", [0, 0, 0])
                else:

                    with warnings.catch_warnings(record=True) as w:

                        # Calculate axis and angle of rotation
                        axis = np.cross(parent_direction, segment_direction)
                        axis /= np.linalg.norm(axis)
                        angle = np.arccos(np.dot(parent_direction.flatten(), segment_direction))
                        rotation = Rotation.from_rotvec(angle * axis)

                        # Verify that we got same point back if we apply rotation
                        # if not np.allclose(segment_direction, rotation.apply(parent_direction), atol=1e-3):
                        #     print(f"sd: {segment_direction},  rotated parent: {rotation.apply(parent_direction)}")
                        #     import pdb
                        #     pdb.set_trace()

                    if w:
                        for wm in w:
                            print(wm.message)
                        import pdb
                        pdb.set_trace()

            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

            rotations_and_length.append((rotation, segment_length))
            parent_direction = segment_direction

        return rotations_and_length, segment_direction

    def coordinate_representation(self, rotation_representation,
                                  parent_direction=None,
                                  parent_point=None,
                                  return_last_direction=False,
                                  include_parent_point=True):

        if parent_direction is None:
            parent_direction = np.array([1, 0, 0])

        if parent_point is None:
            parent_point = np.zeros(3)

        parent_direction = parent_direction / np.linalg.norm(parent_direction)

        coords = np.zeros((len(rotation_representation)+include_parent_point, 3))

        if include_parent_point:
            coords[0, :] = parent_point

        for idx, (rotation, length) in enumerate(rotation_representation):
            segment_direction = rotation.apply(parent_direction)
            parent_point = coords[idx+include_parent_point, :] = segment_direction * length + parent_point
            parent_direction = segment_direction

        if return_last_direction:
            return coords, parent_direction

        return coords

    def write_neuron(self, neuron: NeuronMorphologyExtended, output_file):

        morphology = neuron.morphology_data["neuron"]
        comment = f"{neuron.name} located at {neuron.position}"

        self.write_swc(morphology=morphology, output_file=output_file, comment=comment)

    def write_swc(self, morphology: MorphologyData, output_file, comment=None):
        # We need to write file with micrometer units
        # 0: compartment number (start from 1)
        # 1: compartment type (1-soma, 2-axon, 3-dendrite)
        # 2,3,4: x,y,z
        # 5: r
        # 6: parent compartment

        swc_data = np.zeros((morphology.section_data.shape[0], 7))
        swc_data[:, 0] = np.arange(1, swc_data.shape[0]+1)   # id, start from 1
        swc_data[:, 1] = morphology.section_data[:, 2]       # type
        swc_data[:, 2:5] = (morphology.geometry[:, :3] - morphology.geometry[0, :3]) * 1e6  # x,y,z,
        swc_data[:, 5] = morphology.geometry[:, 3] * 1e6     # radie in micrometer
        swc_data[:, 6] = morphology.section_data[:, 3] + 1   # parent compartment
        swc_data[0, 6] = -1

        # There is a special case, when the first point after the soma is a branch point
        # which could lead to a 1 point section, to handle those Snudda set section_type to 0 for that point
        # https://github.com/neuronsimulator/nrn/blob/5038de0b79ddf7da9b536639989da4c10dbae7f7/share/lib/hoc/import3d/read_swc.hoc#L304
        # We need to find those points and set the section_type to that of the child

        bad_idx_list = np.where(morphology.section_data[:, 2] == 0)[0]
        for bad_idx in bad_idx_list:
            child_idx = np.where(morphology.section_data[:, 3] == bad_idx)[0]
            s_type = morphology.section_data[child_idx, 2]
            assert (s_type == s_type[0]).all(), f"Children of different type in {morphology.swc_file} {s_type}"
            swc_data[bad_idx, 1] = s_type[0]
            print(f"Setting {output_file} row {bad_idx} type to {s_type[0]}.")

        with open(output_file, "wt") as f:
            if comment:
                f.write(f"#{comment}\n")

            for row in swc_data:
                f.write(f"{row[0]:.0f} {row[1]:.0f} {row[2]:.5f} {row[3]:.5f} {row[4]:.5f} {row[5]:.5f} {row[6]:.0f}\n")

        print(f"Wrote {output_file}")

    def edge_avoiding_morphology(self, swc_file, new_file, original_position, original_rotation,
                                 k_dist=30e-6, max_angle=0.1, n_random=5,
                                 random_seed=None):

        md = MorphologyData(swc_file=swc_file)
        md.place(rotation=original_rotation, position=original_position)
        rot_rep, morphology_changed = self.bend_morphology(md,
                                                           k_dist=k_dist, max_angle=max_angle,
                                                           n_random=n_random,
                                                           random_seed=random_seed)

        if morphology_changed:
            new_coord = self.apply_rotation(md, rot_rep)
            md.geometry[:, :3] = new_coord
            self.write_swc(morphology=md, output_file=new_file)
            return new_file

        # Returns None if morphology was not changed
        return None


def test_rotation_representation():

    file_path = "../data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508/WT-0728MSN01-cor-rep-ax.swc"

    md = MorphologyData(swc_file=file_path)

    bm = BendMorphologies(None, rng=np.random.default_rng())
    # sec = md.sections[3][0]
    # rot, _ = bm.rotation_representation(sec)
    # coords = bm.coordinate_representation(rotation_representation=rot, parent_point=sec.position[0, :])

    rot_rep = bm.get_full_rotation_representation(morphology=md)
    coords = bm.apply_rotation(morphology=md, rotation_representation=rot_rep)

    if not (np.abs(md.geometry[:, :3] - coords) < 1e-6).all():
        print(f"Geometry mismatch, problem with representation.")
        import pdb
        pdb.set_trace()
    else:
        print(f"Geometry matches")


def test_write():

    file_path = "../data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508/WT-0728MSN01-cor-rep-ax.swc"
    md = MorphologyData(swc_file=file_path)
    bm = BendMorphologies(None, rng=np.random.default_rng())
    bm.write_swc(morphology=md, output_file="test-write.swc")

    import pdb
    pdb.set_trace()


def test_bending():

    file_path = "../data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508/WT-0728MSN01-cor-rep-ax.swc"
    # file_path = "delme3.swc"
    mesh_path = "../data/mesh/Striatum-d-right.obj"

    nm = NeuronMorphologyExtended(swc_filename=file_path)

    # md = MorphologyData(swc_file=file_path)
    bm = BendMorphologies(mesh_path, rng=np.random.default_rng(1))

    pos = np.array([0.006, 0.004, 0.00205])

    before = nm.clone(position=pos, rotation=np.eye(3))
    after = nm.clone(position=pos, rotation=np.eye(3))

    new_rot_rep, _ = bm.bend_morphology(after.get_morphology())
    new_coord = bm.apply_rotation(after.get_morphology(), new_rot_rep)
    after.get_morphology().geometry[:, :3] = new_coord

    change = np.sum(np.abs(before.get_morphology().geometry[:, :3] - after.get_morphology().geometry[:, :3]))
    print(f"Change = {change}")

    bm.write_neuron(neuron=before, output_file="before-neuron.swc")
    bm.write_neuron(neuron=after, output_file="after-neuron.swc")

    # before.plot_neuron()
    # after.plot_neuron()

    # bm.region_mesh.plot(neurons=[before], show_axis=True, show_faces=False)
    # bm.region_mesh.plot(neurons=[after], show_axis=True, show_faces=False)

    # import pdb
    # pdb.set_trace()


if __name__ == "__main__":

    import cProfile

    # test_write()
    test_rotation_representation()

    profiler = cProfile.Profile()
    profiler.enable()

    test_bending()

    profiler.disable()
    profiler.print_stats(sort='cumulative')

    # import pdb
    # pdb.set_trace()
