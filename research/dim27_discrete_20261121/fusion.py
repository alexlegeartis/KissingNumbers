"""Exact two-code recombination by bipartite matching."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import maximum_bipartite_matching

from source import COMPLEMENT, OCTAD, keys, load_points
from golay import golay24


def complement_group():
    code = golay24().astype(np.int16)
    shortened = code[~np.any(code[:, OCTAD], axis=1)]
    assert len(shortened) == 32
    labels = np.zeros(len(COMPLEMENT), dtype=np.int16)
    for word in shortened:
        new = 2 * labels + word[COMPLEMENT]
        if len(np.unique(new)) > len(np.unique(labels)):
            labels = new
        if len(np.unique(labels)) == len(COMPLEMENT):
            break
    assert set(labels.tolist()) == set(range(16))
    coordinate = dict(zip(labels.tolist(), COMPLEMENT.tolist()))
    codewords = set(keys(code).tolist())
    permutations = []
    for shift in range(16):
        order = np.arange(27)
        order[COMPLEMENT] = [coordinate[label ^ shift] for label in labels.tolist()]
        assert set(keys(code[:, order[:24]]).tolist()) == codewords
        permutations.append(order)
    signs = np.pad(1 - 2 * shortened, ((0, 0), (0, 3)), constant_values=1)
    return permutations, signs


def cross_conflicts(left, right):
    """Both inputs are valid and their zero-tail points belong to the same shell."""
    R = right.astype(np.float32)
    rows, cols = [], []
    tail_left = np.any(left[:, 24:] != 0, axis=1)
    all_right = np.arange(len(right))
    tail_right = all_right[np.any(right[:, 24:] != 0, axis=1)]
    for chosen_left, chosen_right in [(np.flatnonzero(tail_left), all_right),
                                     (np.flatnonzero(~tail_left), tail_right)]:
        part = R[chosen_right]
        for start in range(0, len(chosen_left), 256):
            ix = chosen_left[start:start + 256]
            dot = left[ix].astype(np.float32) @ part.T
            i, j = np.nonzero(dot > 64)
            rows.append(ix[i])
            cols.append(chosen_right[j])
    i, j = np.concatenate(rows), np.concatenate(cols)
    return csr_array((np.ones(len(i), dtype=np.int8), (i, j)), shape=(len(left), len(right)))


def fuse(left, right, left_keys=None):
    left_keys = keys(left) if left_keys is None else left_keys
    right_keys = keys(right)
    same = np.isin(left_keys, right_keys, assume_unique=True)
    A = left[~same]
    B = right[~np.isin(right_keys, left_keys, assume_unique=True)]
    if not len(B):
        return 0, len(A), 0, 0, None
    graph = cross_conflicts(A, B)
    match = maximum_bipartite_matching(graph, perm_type="row")
    matched_right = np.flatnonzero(match >= 0)
    paired_left = match[matched_right]
    assert len(np.unique(paired_left)) == len(matched_right)
    # Independently check the concrete matching using integer arithmetic.
    assert np.all(np.sum(A[paired_left] * B[matched_right], axis=1) > 64)
    gain = len(B) - len(matched_right)
    if not gain:
        return gain, len(A), len(B), graph.nnz, None

    # A matching and a minimum vertex cover furnish the optimum of the union.
    reach_left = np.ones(len(A), dtype=bool)
    reach_left[paired_left] = False
    reach_right = np.zeros(len(B), dtype=bool)
    work = np.flatnonzero(reach_left).tolist()
    while work:
        i = work.pop()
        for j in graph.indices[graph.indptr[i]:graph.indptr[i + 1]]:
            if not reach_right[j] and match[j] != i:
                reach_right[j] = True
                if match[j] >= 0 and not reach_left[match[j]]:
                    reach_left[match[j]] = True
                    work.append(match[j])
    i, j = graph.nonzero()
    assert not (reach_left[i] & ~reach_right[j]).any()
    points = np.r_[left[same], A[reach_left], B[~reach_right]]
    assert len(points) == len(left) + gain
    return gain, len(A), len(B), graph.nnz, points


def initialize():
    global SOURCE, SOURCE_KEYS, PERMUTATIONS, SIGNS
    SOURCE = load_points()
    SOURCE_KEYS = keys(SOURCE)
    PERMUTATIONS, SIGNS = complement_group()


def run(transform):
    shift, sign = transform
    other = SOURCE[:, PERMUTATIONS[shift]] * SIGNS[sign]
    gain, left_size, right_size, edges, points = fuse(SOURCE, other, SOURCE_KEYS)
    return shift, sign, gain, left_size, right_size, edges, points


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    transforms = list(product(range(16), range(32)))
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize) as pool:
        for shift, sign, gain, left, right, edges, points in pool.map(run, transforms):
            rows.append([shift, sign, gain, left, right, edges])
            if gain and args.out:
                np.savez_compressed(args.out.with_name(f"fusion_{shift}_{sign}.npz"), points=points)
    result = {
        "transforms": len(rows),
        "maximum_gain": max(row[2] for row in rows),
        "matched_right_side_in_every_union": all(row[2] == 0 for row in rows),
        "transforms_with_no_coordinate_change": sum(row[4] == 0 for row in rows),
        "largest_difference_on_one_side": max(row[4] for row in rows),
    }
    print(json.dumps(result, indent=2))
    if args.out:
        args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
