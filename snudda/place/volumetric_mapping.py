#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:56:04 2022

@author: wst
"""

import json
import os

import numpy as np
from scipy.spatial.distance import pdist, squareform

from snudda.place.region_mesh import RegionMesh


class VolumetricMap:

    def __init__(self, meshfile_a, meshfile_b, n=100, inversion=False, initial=True, predefined_points=None):
        self.meshfile_a = meshfile_a
        self.meshfile_b = meshfile_b
        self.n = n
        self.inversion = inversion
        self.initial = initial

        if predefined_points:  # user defined nested list of points to map between
            self.a_points = np.array(predefined_points[0])
            self.b_points = np.array(predefined_points[1])
        else:  # highjack RegionMesh to sample the volume
            self.a = RegionMesh(meshfile_a)
            self.b = RegionMesh(meshfile_b)
            self.a_points = self.a.place_neurons(self.n)
            self.b_points = self.b.place_neurons(self.n)

    def distance(self, coords):

        return squareform(pdist(coords))

    def initial_map(self):

        bcopy = self.b_points.copy()
        if self.inversion:  # invert ventrodorsal
            bcopy[:, 2] *= -1

        mapping = []
        for i in range(len(self.a_points)):
            j = np.linalg.norm(self.a_points[i] - bcopy, axis=1).argmin()
            k = np.argwhere(self.b_points == bcopy[j])[0, 0]
            mapping.append(k)
            bcopy = np.delete(bcopy, j, axis=0)

        return mapping

    def energy(self, mapping):

        a_dist = self.distance(self.a_points)
        if self.inversion:  # invert ventrodorsal topography
            inverted_points = self.b_points.copy()
            inverted_points[:, 2] *= -1
            # inverted_points[:,1]*= 1
            b_dist = self.distance(inverted_points)
        else:
            b_dist = self.distance(self.b_points)

        return np.sum(np.power(a_dist - b_dist[:, mapping][mapping, :], 2))

    def optimize_map(self, n_iter):

        if self.initial:
            self.mapping = self.initial_map()
        else:
            self.mapping = np.arange(0, len(self.a_points), 1)

        history = []  # track optimization
        e1 = self.energy(self.mapping)
        e = e1
        history.append(e)

        m = np.random.randint(low=0, high=self.n, size=(n_iter, 2))
        p_val = np.random.uniform(size=n_iter)

        for i, ((m1, m2), px) in enumerate(zip(m, p_val)):

            new = self.mapping.copy()
            new[m1] = self.mapping[m2]
            new[m2] = self.mapping[m1]
            e_new = self.energy(new)
            e_change = e_new - e

            p_swap = 1 / (1 + np.exp(50 * e_change / e1 * 100))

            if px < p_swap:
                self.mapping = new
                e = e_new
            history.append(e)

        self.opt_history = np.array(history)
        self.a_points = self.a_points
        self.b_points = self.b_points[self.mapping]

        return

    def plot_energy(self):

        import matplotlib.pyplot as plt
        if self.opt_history is None:
            self.optimize_map(1000)

        fig, ax = plt.subplots()
        ax.plot(self.opt_history)
        ax.set_xlabel('Iterations')
        ax.set_ylabel('Energy')
        title = '{} points, {} iterations'.format(self.n, len(self.opt_history) - 1)
        ax.set_title(title)

        return

    def get_colour(self, coords, idx):
        min_coord = np.min(coords)
        max_coord = np.max(coords)

        return np.divide(coords[idx, :] - min_coord, max_coord - min_coord)

    def visualize_map(self, lines=True, points=True):

        import trimesh
        from vedo import embedWindow
        embedWindow(backend=False)
        from vedo import Mesh
        from vedo import Line
        import brainrender
        brainrender.settings.SHADER_STYLE = 'plastic'
        brainrender.settings.SHOW_AXES = False
        from brainrender.scene import Scene

        a_tri = trimesh.load_mesh(self.meshfile_a)
        b_tri = trimesh.load_mesh(self.meshfile_b)

        scene = Scene(root=True)
        scene.jupyter = True

        scene.add(Mesh(a_tri, c='red', alpha=0.1))
        scene.add(Mesh(b_tri, c='blue', alpha=0.1))

        if lines:
            for i in range(self.n):
                col = self.get_colour(self.a_points[:], i)
                scene.add(Line(self.a_points[i, :] * 1e6, self.b_points[i, :] * 1e6, lw=2.5, c=col))
        if points:
            for i in range(self.n):
                col = self.get_colour(self.a_points[:], i)
                scene.add(Point(self.a_points[i, :] * 1e6, radius=25, color=col, alpha=0.4))
                scene.add(Point(self.b_points[i, :] * 1e6, radius=25, color=col, alpha=0.4))

        scene.slice(scene.atlas.get_plane(norm=[0, 0, 1], sx=15000, sy=15000), close_actors=False)
        scene.render(camera='top')

        return

    def write_json(self, projection_name, outfile):

        mapping = {projection_name: {'source': self.a_points.tolist(), 'destination': self.b_points.tolist()}}

        with open(outfile, 'w+') as f:
            json.dump(mapping, f, indent=4)

        return


if __name__ == "__main__":
    source_file = os.path.join("$SNUDDA_DATA", "mesh", "striatum-d.obj")
    target_file = os.path.join("$SNUDDA_DATA", "mesh", "SNr.obj")
    outfile = os.path.join("$SNUDDA_DATA", "projections", "direct.json")
    projection_name = 'DirectPathway'

    vm = VolumetricMap(source_file, target_file, n=50, inversion=True)
    vm.optimize_map(10000)
    vm.plot_energy()
    vm.visualize_map()
    vm.write_json(projection_name=projection_name, outfile=outfile)
