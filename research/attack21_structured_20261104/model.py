"""Exact graph model for the line-weighted 21D construction.

Geometry: Kevin Li, https://hexagonmath.org/2610.00248v1, Section 3,
CC BY 4.0. Inputs were independently rebuilt by the preceding campaign.
"""
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'
OLD = HERE.parent / 'compute19_21_20261104' / 'data'
Z = np.load(OLD / 'component21.npz')
C = Z['component'].astype(int)
N = Z['neighbors'].astype(int)
P = Z['line_ids'].astype(int)
B = np.load(OLD / 'baseline21.npz')
PARITY = np.array([int(x).bit_count() % 2 for x in range(4096)], dtype=np.int8)
SIDE = np.array([int(B['words'][x]).bit_count() % 2 for x in C])
assert N.shape == (1024, 21)
assert np.array_equal(C[N], C[:, None] ^ P)
assert np.all(SIDE[N] != SIDE[:, None])


def labels_from_holes(holes, rng=None):
    holes = np.asarray(holes, dtype=bool)
    assert not holes[N[holes]].any(), 'holes must be independent'
    choices = holes[N]
    assert np.all(holes | choices.any(axis=1)), 'holes must dominate'
    if rng is None:
        labels = choices.argmax(axis=1)
    else:
        labels = (rng.random(N.shape) * choices).argmax(axis=1)
    labels[holes] = -1
    verify_labels(labels)
    return labels


def verify_labels(labels):
    labels = np.asarray(labels)
    assert labels.shape == (1024,)
    assert np.all((labels >= -1) & (labels < 21))
    for u in np.flatnonzero(labels >= 0):
        p = int(labels[u])
        assert labels[N[u, p]] == -1
        assert np.all(labels[N[u]] != p)
    return int(np.sum(labels >= 0))


def verify_component_geometry(labels):
    selected = np.flatnonzero(labels >= 0)
    rays = B['profiles'][labels[selected]] * (1 - 2 * B['code'][C[selected], :21])
    dots = rays.astype(np.int64) @ rays.astype(np.int64).T
    np.fill_diagonal(dots, -164)
    assert dots.max() <= 82
    return dict(selected=len(selected), unordered_pairs=len(selected)*(len(selected)-1)//2,
                maximum_dot=int(dots.max()), squared_norm=164)
