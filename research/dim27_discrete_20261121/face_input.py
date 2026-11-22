"""Prepare the contact graph on the 13,248 deleted roots and exact tail-radius tables."""
import argparse
from itertools import product
from pathlib import Path

import numpy as np

from bw_refill import TAILS
from source import checked_archive, leech_shell, load_points, verify_recipe


def radius_table(k):
    required = 4 * k * (k - 1)
    limit = {3: 7, 4: 10}[k]
    digits = {3: 8, 4: 16}[k]
    # Each vertex solves three signed root-plane equalities, with determinant 1 or 2.
    twice = np.array(list(product(range(-2 * limit, 2 * limit + 1), repeat=3)), dtype=np.int32)
    bounds = (np.abs(twice @ (TAILS // 4).T) + 1) // 2
    keep = (bounds <= limit).all(1) & ((twice * twice).sum(1) >= 4 * required)
    table = np.zeros((limit + 1,) * 6, dtype=np.uint8)
    table[tuple(bounds[keep].T)] = 1
    for axis in range(6):
        table = np.maximum.accumulate(table, axis=axis)
    padded = np.zeros((digits,) * 6, dtype=np.uint8)
    padded[(slice(limit + 1),) * 6] = table
    return padded.T.copy()


def prepare(directory):
    directory.mkdir(parents=True, exist_ok=True)
    shell = 2 * leech_shell().astype(np.int16)
    verify_recipe(load_points(), shell)
    with checked_archive("d27_certificate.npz") as archive:
        mids = archive["added_integer"][archive["midpoint_rows"]]
        parents = np.lexsort(shell.T[::-1])[archive["midpoint_endpoints"]]
    keep = mids[np.arange(len(mids)), 24 + np.argmax(mids[:, 24:] != 0, axis=1)] > 0
    mids, parents = mids[keep], parents[keep]
    lookup = {tuple(t): i for i, t in enumerate(TAILS)}
    labels = np.array([lookup[tuple(t)] for t in mids[:, 24:]], dtype="<i4")
    order = np.argsort(labels)
    mids, parents, labels = mids[order], parents[order], labels[order]
    roots, inverse = np.unique(parents, return_inverse=True)
    assert len(roots) == len(mids) == 13248
    incident = np.argsort(inverse.ravel()).reshape(-1, 2) // 2
    assert np.all(np.bincount(inverse.ravel()) == 2)
    incident.astype("<i4").tofile(directory / "incident")
    labels.tofile(directory / "colors")
    coordinates = shell[roots].astype(np.float32)
    scores = np.memmap(directory / "products", mode="w+", dtype=np.int8, shape=(len(roots), len(mids)))
    transposed = mids[:, :24].astype(np.float32).T
    for i in range(0, len(roots), 128):
        gram = (coordinates[i:i + 128] @ transposed).astype(np.int16)
        assert not np.any(gram % 16) and np.max(np.abs(gram)) <= 96
        scores[i:i + 128] = gram // 16
    scores.flush()
    groups = [np.flatnonzero(labels == c) for c in range(6)]
    with (directory / "strongest").open("wb") as output:
        for row in scores:
            strongest = [g[np.argsort(-row[g], kind="stable")[:48]] for g in groups]
            np.array(strongest, dtype="<i4").T.copy().tofile(output)
    degrees = []
    with (directory / "edges").open("wb") as output, (directory / "root_products").open("wb") as products:
        for i in range(0, len(roots), 192):
            gram = (coordinates[i:i + 192] @ coordinates.T).astype(np.int16)
            assert not np.any(gram % 32)
            (gram // 32).astype(np.int8).tofile(products)
            a, b = np.nonzero(np.triu(gram[:, i:] == 64, 1))
            (b + i).astype("<i4").tofile(output)
            degrees.extend(np.bincount(a, minlength=len(gram)))
    offsets = np.r_[0, np.cumsum(degrees)].astype("<i4")
    assert offsets[-1] == 2966566
    offsets.tofile(directory / "offsets")
    for k in (3, 4):
        radius_table(k).tofile(directory / f"radius{k}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    prepare(parser.parse_args().directory)
