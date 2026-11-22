"""Prepare integer dot products, all Leech contacts and exact feasible tail profiles."""

import argparse
from itertools import combinations, product
from pathlib import Path

import numpy as np

from bw_refill import TAILS
from source import checked_archive, leech_shell, load_points, verify_recipe


def tail_profiles():
    normals = TAILS.astype(np.int64) // 4
    triples, transforms, determinants = [], [], []
    for indices in combinations(range(6), 3):
        a, b, c = normals[list(indices)]
        adj = np.stack([np.cross(b, c), np.cross(c, a), np.cross(a, b)], axis=1)
        det = int(a @ adj[:, 0])
        if det == 0:
            continue
        assert np.array_equal(normals[list(indices)] @ adj, det * np.eye(3, dtype=int))
        for signs in product((-1, 1), repeat=3):
            triples.append(indices)
            transforms.append(adj * np.sign(det) * signs)
            determinants.append(abs(det))
    triples = np.array(triples)
    transforms = np.array(transforms)
    determinants = np.array(determinants)
    minimal = []
    for profile in sorted(product(range(5), repeat=6), key=sum):
        bounds = np.array(profile)
        if minimal and np.any((np.array(minimal) <= bounds).all(1)):
            continue
        numerators = (transforms @ bounds[triples][:, :, None])[:, :, 0]
        in_polytope = (np.abs(numerators @ normals.T) <= bounds * determinants[:, None]).all(1)
        outside_sphere = (numerators * numerators).sum(1) >= 8 * determinants**2
        if np.any(in_polytope & outside_sphere):
            minimal.append(bounds)
    return np.array(minimal, dtype="<i4")


def prepare(directory):
    directory.mkdir(parents=True, exist_ok=True)
    shell = 2 * leech_shell().astype(np.int16)
    verify_recipe(load_points(), shell)
    coords = shell.astype(np.float32)
    with checked_archive("d27_certificate.npz") as archive:
        mids = archive["added_integer"][archive["midpoint_rows"]]
    first = np.argmax(mids[:, 24:] != 0, axis=1)
    mids = mids[mids[np.arange(len(mids)), 24 + first] > 0]
    line = {tuple(t): i for i, t in enumerate(TAILS)}
    labels = np.array([line[tuple(t)] for t in mids[:, 24:]])
    mids = mids[np.argsort(labels)]
    counts = np.r_[0, np.cumsum(np.bincount(labels))]
    assert len(mids) == 13248 and len(counts) == 7
    with (directory / "capmeta").open("wb") as f:
        np.array([len(shell), len(mids), *counts], dtype="<i4").tofile(f)
    # Integer products and absolute partial sums fit in binary32.
    with (directory / "cap.bin").open("wb") as f:
        for start in range(0, len(mids), 128):
            dots = (mids[start:start + 128, :24].astype(np.float32) @ coords.T).astype(np.int16)
            assert not np.any(dots % 16) and dots.max() <= 96 and dots.min() >= -96
            (dots // 16).astype(np.int8).tofile(f)
    degree = []
    with (directory / "live_neighbors").open("wb") as f:
        for start in range(0, len(coords), 256):
            rows = coords[start:start + 256]
            gram = rows @ coords[start:].T
            a, b = np.nonzero(np.triu(gram == 64, 1))
            (b + start).astype("<i4").tofile(f)
            degree.extend(np.bincount(a, minlength=len(rows)))
    offsets = np.r_[0, np.cumsum(degree)].astype("<i8")
    assert len(offsets) == len(shell) + 1 and offsets[-1] == 452088000
    offsets.tofile(directory / "live_offsets")
    profiles = tail_profiles()
    with (directory / "profiles").open("wb") as f:
        np.array([len(profiles)], dtype="<i4").tofile(f)
        profiles.tofile(f)
    print(f"{len(shell)} shell points, {offsets[-1]} contacts, {len(profiles)} minimal profiles")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    prepare(parser.parse_args().directory)
