"""Rebuild equal-cardinality 27D witnesses and distinguish them by contact count."""

import argparse
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from fibres import native_fibres
from source import COMPLEMENT, keys, leech_shell, load_points


def make_heads(frame, axis):
    heads = np.zeros((64, 27), dtype=np.int16)
    heads[:, COMPLEMENT] = np.tile(frame, (2, 1))
    heads[:, 24 + axis] = np.repeat([8, -8], 32)
    return heads


def enumerate_frames(output=None):
    source = load_points()
    frames, _ = native_fibres(2 * leech_shell().astype(np.int16))
    assert len(np.unique(keys(frames.reshape(-1, 16)))) == 4320
    is_head = ((source[:, COMPLEMENT] ** 2).sum(1) == 64) & (
        (source[:, 24:] ** 2).sum(1) == 64
    )
    fixed = source[~is_head]
    removed = source[is_head]
    assert len(removed) == 192
    outside = fixed[:, COMPLEMENT].astype(np.float32).T
    costs = [{} for _ in range(3)]
    original = [None] * 3
    for index, frame in enumerate(frames):
        gram = frame @ frame.T
        assert np.all(np.diag(gram) == 64)
        assert np.all(np.sort(gram, axis=1) == [-64] + [0] * 30 + [64])
        # Only integer products of small coordinates occur, exactly representable in binary32.
        body = frame.astype(np.float32) @ outside
        for axis in range(3):
            offset = 8 * fixed[:, 24 + axis]
            if np.any(body + np.abs(offset) > 64):
                continue
            costs[axis][index] = int(np.sum(body + offset == 64) + np.sum(body - offset == 64))
            old = removed[removed[:, 24 + axis] != 0]
            if set(keys(old).tolist()) == set(keys(make_heads(frame, axis)).tolist()):
                original[axis] = index
    assert all(index is not None for index in original)

    def contacts(choices):
        # Within each axis, opposite tails never touch; each signed tail contributes 480 pairs.
        cross = sum(128 for a, b in combinations(choices, 2) if a == b)
        return 3 * 960 + cross + sum(costs[i][c] for i, c in enumerate(choices))

    choices = [(contacts(s), s) for s in product(*costs)]
    lowest, example = min(choices)
    if output:
        points = np.concatenate([fixed] + [make_heads(frames[j], i) for i, j in enumerate(example)])
        np.savez_compressed(output, points=points, frames=example)
    histogram = Counter(count for count, _ in choices)
    return {
        "points": len(source),
        "compatible_frames_per_axis": list(map(len, costs)),
        "frame_triples": len(choices),
        "distinct_contact_totals": len(histogram),
        "source_frames": original,
        "source_contacts_incident_to_heads": contacts(original),
        "minimum_contact_frames": example,
        "minimum_contacts_incident_to_heads": lowest,
        "maximum_contacts_incident_to_heads": max(histogram),
        "contact_histogram": dict(sorted(histogram.items())),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    print(json.dumps(enumerate_frames(parser.parse_args().out), indent=2))
