"""Exact finite-cutoff neighbors without dense all-pairs arrays.

KD-tree results are sorted so candidate event clocks are consumed in particle
order. Cost depends on the number of actual local contacts; a deliberately
overcrowded configuration can still contain quadratically many contacts.
"""

import itertools
import math

import numpy as np
from scipy.spatial import cKDTree


def neighbor_pairs(positions, cutoff, tree=None):
    """Return lexicographically sorted i<j pairs strictly inside ``cutoff``."""
    if len(positions) < 2 or cutoff <= 0:
        return np.empty((0, 2), dtype=np.intp)
    tree = cKDTree(positions) if tree is None else tree
    pairs = tree.query_pairs(float(cutoff), output_type="ndarray")
    if not pairs.size:
        return np.empty((0, 2), dtype=np.intp)
    delta = positions[pairs[:, 0]] - positions[pairs[:, 1]]
    pairs = pairs[np.einsum("ij,ij->i", delta, delta) < cutoff * cutoff]
    if not pairs.size:
        return np.empty((0, 2), dtype=np.intp)
    return pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]


class GrowingNeighbors:
    """Incremental spatial hash for overlap-free initial chain growth."""

    def __init__(self, cell_size, max_radius):
        self.cell_size = float(cell_size)
        self.max_radius = float(max_radius)
        self.cells = {}
        self.entries = {}

    def _key(self, position):
        return tuple(np.floor(np.asarray(position) / self.cell_size).astype(int))

    def add(self, index, position, radius):
        point = np.asarray(position, dtype=float).copy()
        key = self._key(point)
        self.entries[index] = (point, float(radius), key)
        self.cells.setdefault(key, []).append(index)

    def remove(self, index):
        _, _, key = self.entries.pop(index)
        self.cells[key].remove(index)
        if not self.cells[key]:
            del self.cells[key]

    def is_clear(self, position, radius, fraction=0.65):
        point = np.asarray(position)
        key = self._key(point)
        reach = int(math.ceil(fraction * (float(radius) + self.max_radius) / self.cell_size))
        for offset in itertools.product(range(-reach, reach + 1), repeat=3):
            nearby_key = tuple(key[d] + offset[d] for d in range(3))
            for index in self.cells.get(nearby_key, ()):
                other, other_radius, _ = self.entries[index]
                delta = point - other
                if np.dot(delta, delta) <= (fraction * (radius + other_radius)) ** 2:
                    return False
        return True
