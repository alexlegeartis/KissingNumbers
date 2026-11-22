"""Exhaust all last two normal vectors for the certified ordered first pair."""
from collections import Counter
from pathlib import Path

import numpy as np


def main():
    with np.load(Path(__file__).with_name("ambient.npz"), allow_pickle=False) as data:
        x = data["section_points"].astype(np.int64)
        gram = data["gram"].astype(np.int64)
        first = data["normals"][:2].astype(np.int64)
    dots = (x @ gram) @ x.T
    histogram = Counter()
    best = None
    for i, row in enumerate(dots):
        zero = np.flatnonzero(row == 0)
        three = np.flatnonzero(row == 3)
        candidates = zero[zero > i]
        if not len(candidates):
            continue
        pairs0 = dots[np.ix_(zero, candidates)]
        pairs3 = dots[np.ix_(three, candidates)]
        counts = np.array([(pairs3 == 3).sum(0), (pairs3 == 0).sum(0), (pairs0 == 3).sum(0)])
        sums = counts.sum(0)
        histogram.update(map(int,sums))
        j = int(np.argmin(sums))
        result = (int(sums[j]), i, int(candidates[j]))
        if best is None or result < best:
            best = result
    defect, i, j = best
    normals = np.vstack([first, -x[i], -x[j]])
    assert np.array_equal(normals @ gram @ normals.T,
                          [[6,-3,-3,-3],[-3,6,0,0],[-3,0,6,0],[-3,0,0,6]])
    print("Orthogonal completions:", sum(histogram.values()))
    print("Completions attaining the best count:", histogram[defect])
    print("Minimum sum of the three shell counts:", defect)
    print("Configuration size:", 4_940_961-70*defect)
    print("Rows of the first maximizing pair in section_points:", i, j)


if __name__ == "__main__":
    main()
