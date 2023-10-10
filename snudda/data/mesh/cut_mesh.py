import open3d as o3d
import numpy as np


def cut_mesh(path, save=False):
    mesh = o3d.io.read_triangle_mesh(path)
    x, y, z = np.array(mesh.vertices).T
    left = z > np.mean(z)  # apparent mediolateral axis
    left_idxs = np.where(left)[0]
    left_mesh = mesh.select_by_index(left_idxs)
    right_idxs = np.where(~left)[0]
    right_mesh = mesh.select_by_index(right_idxs)

    if save:
        base = path.strip('.obj')
        o3d.io.write_triangle_mesh(filename=base+'-left.obj', mesh=left_mesh)
        o3d.io.write_triangle_mesh(filename=base+'-right.obj', mesh=right_mesh)

    return mesh, right_mesh, left_mesh

if __name__ == '__main__':
    meshes = ['Striatum-v.obj', 'Striatum-d.obj', 'GPe.obj', 'GPi.obj', 'STN.obj', 'SNr.obj', 'SNc.obj']
    for mesh_name in meshes:
        cut_mesh(mesh_name, save=True)
