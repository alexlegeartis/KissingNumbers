"""Bound the 27D norm over polar candidates with bodies on any Leech root ray."""

import argparse
from fractions import Fraction
from itertools import combinations, permutations
import json
from pathlib import Path

import numpy as np

from bw_refill import TAILS
from midpoint_switches import source_layers
from source import leech_shell


def determinant(matrices):
    """Batched exact 4-by-4 determinants. Entries are bounded by 128 here."""
    det = np.zeros(len(matrices), dtype=np.int64)
    for perm in permutations(range(4)):
        inversions = sum(perm[i] > perm[j] for i, j in combinations(range(4), 2))
        product = np.prod(matrices[:, np.arange(4), perm], axis=1)
        det += (-1) ** inversions * product
    return det


def largest_norm(profile, normals, active):
    # Variables are the root multiplier a and three Cartesian tail coordinates.
    matrix = np.vstack([np.c_[profile, normals], np.c_[profile, -normals],
                        [-3, 0, 0, 0], [4, 0, 0, 0]])
    rhs = np.r_[np.full(18, 64), -1, 3]
    constraints, target = matrix[active], rhs[active]
    det = determinant(constraints)
    keep = det != 0
    det, constraints, target = det[keep], constraints[keep], target[keep]
    coordinates = []
    for j in range(4):
        m = constraints.copy()
        m[:, :, j] = target
        coordinates.append(determinant(m))
    numerators = np.stack(coordinates, axis=1) * np.sign(det)[:, None]
    denominators = abs(det)
    assert np.all((constraints @ numerators[:, :, None])[:, :, 0] == denominators[:, None] * target)
    keep = (numerators @ matrix.T <= denominators[:, None] * rhs).all(1)
    values = (numerators[keep] ** 2) @ [128, 1, 1, 1]
    return max(Fraction(int(n), int(d) ** 2) for n, d in zip(values, denominators[keep]))


def run(directory):
    mid, labels, heads = source_layers()
    shell = 2 * leech_shell().astype(np.int16)
    products = np.memmap(directory / "cap.bin", mode="r", dtype="i1", shape=(len(mid), len(shell)))
    scores = []
    for i in range(6):
        rows = np.flatnonzero(labels == i)
        assert np.all(np.diff(rows) == 1)
        scores.append(16 * np.max(products[rows[0]:rows[-1] + 1], axis=0))
    for i in range(3):
        block = heads[heads[:, 24 + i] == 8, :24]
        assert len(block) == 48
        scores.append(np.max(shell.astype(np.float32) @ block.astype(np.float32).T, axis=1))
    profiles = np.unique(np.stack(scores, axis=1).astype(np.int64), axis=0)
    normals = np.r_[TAILS, 8 * np.eye(3, dtype=np.int64)]
    active = np.array(list(combinations(range(20), 4)))
    maximum = max(largest_norm(p, normals, active) for p in profiles)
    return {"root_rays": len(shell), "support_profiles": len(profiles),
            "multiplier_interval": ["1/3", "3/4"], "largest_squared_norm": str(maximum)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="output of contact_tail_input.py")
    print(json.dumps(run(parser.parse_args().directory), indent=2))
