#!/usr/bin/env python3

import numpy as np
import copy
from scipy.spatial.distance import pdist, squareform


class PointCluster:

    def __init__(self):
        self.coords = None
        self.dist = None
        self.num = 0

    def generate_cube(self, centre, side, num):
        self.coords = side * np.random.uniform(size=(num, 3))
        self.coords[:, 0] += centre[0]
        self.coords[:, 1] += centre[1]
        self.coords[:, 2] += centre[2]

        self.num = num
        self.calc_distance()

    def generate_disk(self, centre, side, num):
        self.coords = side * np.random.uniform(size=(num, 3))
        self.coords[:, 0] += centre[0]
        self.coords[:, 1] += centre[1]
        self.coords[:, 2] = centre[2]

        self.num = num
        self.calc_distance()

    def calc_distance(self):
        self.dist = squareform(pdist(self.coords))

    def get_colour(self, idx):
        min_coord = np.min(self.coords)
        max_coord = np.max(self.coords)

        return np.divide(self.coords[idx, :] - min_coord, max_coord - min_coord)


class PointMap:

    def __init__(self, source, dest):

        self.source = source
        self.dest = dest

        self.map = np.arange(source.num)

        assert source.num == dest.num, f"source and dest should have equal number of points"

        self.start_energy = None

    def delta_energy_WRONG(self, a, b):

        # For each pair of neurons in the source structure, find their distance, and the distance of their projections.
        # The sum of the product of these distances is the total energy.
        # So the question is how does this sum change if we would swap the neurons.

        a_dest = self.map[a]
        b_dest = self.map[b]

        test_map = copy.deepcopy(self.map)
        test_map[a] = self.map[b]
        test_map[b] = self.map[a]

        energy_before = np.sum(np.multiply(self.source.dist[a, :], self.dest.dist[a_dest, :][test_map])) \
                        + np.sum(np.multiply(self.source.dist[b, :], self.dest.dist[b_dest, :][test_map]))

        energy_after = np.sum(np.multiply(self.source.dist[a, :], self.dest.dist[b_dest, :][test_map])) \
                       + np.sum(np.multiply(self.source.dist[b, :], self.dest.dist[a_dest, :][test_map]))

        return energy_after - energy_before

    def delta_energy(self, a, b):

        energy_before = self.energy()
        alt_map = copy.deepcopy(self.map)
        alt_map[a] = self.map[b]
        alt_map[b] = self.map[a]
        energy_after = self.energy(alt_map)

        return energy_after - energy_before

    def energy(self, alt_map=None):

        if alt_map is not None:
            return np.sum(np.power(self.source.dist - self.dest.dist[:, alt_map][alt_map, :], 2))

        return np.sum(np.power(self.source.dist - self.dest.dist[:, self.map][self.map, :], 2))

    def iterate(self, n_iter):

        if self.start_energy is None:
            self.start_energy = self.energy()

        # Pick a pair of points in source, check if their dest points should be swapped
        ab = np.random.randint(low=0, high=self.source.coords.shape[0], size=(n_iter, 2))
        p_val = np.random.uniform(size=n_iter)

        for idx, ((a, b), px) in enumerate(zip(ab, p_val)):
            energy_change = self.delta_energy(a, b)
            p_swap = 1 / (1 + np.exp(4 * energy_change / self.start_energy * 100))
            print(f"{idx}/{n_iter}  energy = {self.energy()}, p_swap = {p_swap}")
            if px < p_swap:
                # Swap the two if energy
                old_a = self.map[a]
                self.map[a] = self.map[b]
                self.map[b] = old_a

            if idx % 100000 == 0:
                self.plot_map()

        print(f"Energy = {self.energy()}")
        self.plot_map()

    def plot_map(self):

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for idx in range(self.source.num):
            col = self.source.get_colour(idx)
            sp = self.source.coords[idx, :]
            dp = self.dest.coords[self.map[idx], :]
            ax.plot(xs=sp[0], ys=sp[1], zs=sp[2], color=col, marker='o')
            ax.plot(xs=dp[0], ys=dp[1], zs=dp[2], color=col, marker='o')

        fig.show()
        plt.pause(0.1)


if __name__ == "__main__":
    pc_a = PointCluster()
    pc_b = PointCluster()

    pc_a.generate_disk(centre=np.array([0, 0, 0]), side=1000e-6, num=100)
    pc_b.generate_disk(centre=np.array([0, 0, 5e-3]), side=1000e-6, num=100)

    pm = PointMap(source=pc_a, dest=pc_b)

    pm.iterate(10000)

    import pdb

    pdb.set_trace()
