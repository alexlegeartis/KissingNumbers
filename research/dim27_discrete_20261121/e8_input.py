"""Prepare finite E8 head-frame and midpoint data for the exact degree bound."""
import argparse
from itertools import combinations, product
from pathlib import Path

import numpy as np

from fibres import native_fibres
from source import COMPLEMENT, OCTAD, RESIDUAL, leech_shell, load_points


def positive(rows):
    return rows[np.arange(len(rows)), np.argmax(rows != 0, axis=1)] > 0


def canonical(rows):
    return np.unique(rows * np.where(positive(rows), 1, -1)[:, None], axis=0).astype(np.int16)


def clear_exterior(rows, exterior):
    matrix = exterior.astype(np.float32)
    permitted = []
    for k in range(0, len(rows), 256):
        v = rows[k:k + 256].astype(np.float32)
        total = np.abs(v[:, :8] @ matrix[:, :8].T) + np.abs(v[:, 8:] @ matrix[:, 8:].T)
        permitted.extend((total <= 64).all(axis=1))
    return np.array(permitted)


def midpoint_options(roots, exterior):
    parent_lines = {tuple(row): i for i, row in enumerate(roots[positive(roots)])}
    parent_lines.update({tuple(-np.array(row)): i for row, i in list(parent_lines.items())})
    i, j = np.where(np.triu(roots @ roots.T == 64, 1))
    bodies = (roots[i] + roots[j]) // 2
    pos = positive(bodies)
    bodies = bodies[pos]
    parents = np.array([(parent_lines[tuple(roots[a])], parent_lines[tuple(roots[b])])
                        for a, b in zip(i[pos], j[pos])])
    tails = np.array([row for row in product((-4, 0, 4), repeat=3)
                      if sum(x != 0 for x in row) == 2 and next(x for x in row if x) > 0])
    options = np.c_[np.repeat(bodies, len(tails), axis=0), np.tile(tails, (len(bodies), 1))]
    parents = np.repeat(parents, len(tails), axis=0)
    keep = clear_exterior(options, exterior)
    return options[keep], parents[keep]


def head_frame_catalogue(shell, roots, kernel, exterior):
    _, native = native_fibres(shell)
    i, j = np.where(np.triu(roots @ roots.T == 0, 1))
    head_lines = canonical((roots[i] + roots[j]) // 2)
    axes = 8 * np.eye(3, dtype=np.int16)
    hadamard = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])
    catalogues = []
    for axis in range(3):
        tail = axes[axis]
        incumbent = canonical(kernel[(kernel[:, 8:] == tail).all(axis=1), :8])
        assert len(incumbent) == 8
        eligible = head_lines[clear_exterior(np.c_[head_lines, np.tile(tail, (len(head_lines), 1))], exterior)]
        allowed = {tuple(row) for row in eligible}
        found = {incumbent.tobytes(): incumbent}

        def add(frame):
            frame = canonical(frame)
            if all(tuple(row) in allowed for row in frame):
                assert np.array_equal(frame @ frame.T, 64 * np.eye(8, dtype=np.int16))
                found.setdefault(frame.tobytes(), frame)

        for frame in native:
            add(frame)
        for block in combinations(range(8), 4):
            for signs in ([1, 1, 1, 1], [1, 1, 1, -1]):
                frame = incumbent.copy()
                new = hadamard @ (frame[list(block)] * np.array(signs)[:, None])
                assert not np.any(new % 2)
                frame[list(block)] = new // 2
                add(frame)
        catalogues.append(np.array(list(found.values())))
    return catalogues


def build_input(output):
    shell = 2 * leech_shell().astype(np.int16)
    source = load_points()
    elsewhere = np.any(source[:, COMPLEMENT] != 0, axis=1)
    exterior = np.unique(source[elsewhere][:, RESIDUAL], axis=0)
    kernel = source[~elsewhere][:, RESIDUAL]
    roots = shell[~np.any(shell[:, COMPLEMENT], axis=1)][:, OCTAD]
    points, parents = midpoint_options(roots, exterior)
    frames = head_frame_catalogue(shell, roots, kernel, exterior)
    masks = []
    for axis, catalogue in enumerate(frames):
        dot = np.abs(points[:, :8].astype(np.float32) @ catalogue.reshape(-1, 8).astype(np.float32).T)
        dot = dot.reshape(len(points), len(catalogue), 8).max(axis=2)
        masks.append(((points[:, 8 + axis, None] == 0) | (dot <= 32)).T)
    assert [len(a) for a in masks] == [144, 275, 144]
    assert len(points) == 5696
    with output.open("wb") as file:
        np.array([len(points), *(len(a) for a in masks)], dtype="<i4").tofile(file)
        parents.astype("<i4").tofile(file)
        points.astype("<i4").tofile(file)
        for mask in masks:
            mask.astype("u1").tofile(file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    build_input(parser.parse_args().output)
