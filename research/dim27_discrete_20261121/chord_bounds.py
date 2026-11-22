"""Check the whole chord between each pair of deleted contacting Leech roots."""

import argparse
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from bw_refill import TAILS


def possible(bounds, slopes, interval):
    """Test all tail-polytope vertices as affine functions on a rational interval."""
    normals = TAILS.astype(np.int64) // 4
    count = len(bounds)
    feasible = np.zeros(count, dtype=bool)
    for indices in combinations(range(6), 3):
        a, b, c = normals[list(indices)]
        adj = np.stack([np.cross(b, c), np.cross(c, a), np.cross(a, b)], axis=1)
        den = int(a @ adj[:, 0])
        if den == 0:
            continue
        adj *= np.sign(den)
        den = abs(den)
        for signs in product((-1, 1), repeat=3):
            transform = adj * signs
            x = bounds[:, indices] @ transform.T
            y = slopes[:, indices] @ transform.T
            # The vertex equals (x + lambda*y) / den. Intersect all twelve inequalities.
            lo, ld, hi, hd = [np.full(count, v, dtype=np.int64) for v in interval]
            valid = np.ones(count, dtype=bool)
            for sign in (-1, 1):
                intercept = den * bounds - sign * (x @ normals.T)
                slope = den * slopes - sign * (y @ normals.T)
                for g, s in zip(intercept.T, slope.T):
                    valid &= (s != 0) | (g >= 0)
                    lower = (s > 0) & (-g * ld > lo * s)
                    upper = (s < 0) & (g * hd < -hi * s)
                    lo[lower], ld[lower] = -g[lower], s[lower]
                    hi[upper], hd[upper] = g[upper], -s[upper]
            valid &= lo * hd <= hi * ld
            # Squared vertex norm minus 8*lambda*(1-lambda) is convex.
            for num, divisor in [(lo, ld), (hi, hd)]:
                vertex = x * divisor[:, None] + y * num[:, None]
                feasible |= valid & ((vertex**2).sum(1) >= 8 * num * (divisor - num) * den**2)
    return feasible


def check(directory):
    packed = np.fromfile(directory / "chords", dtype="V18")
    assert len(packed) == 2966566
    patterns, counts = np.unique(packed, return_counts=True)
    rows = patterns.view(np.uint8).reshape(-1, 6, 3).astype(np.int64)
    total = rows[:, :, 0] - 16
    feasible = np.zeros(len(rows), dtype=bool)
    for field, interval in [(1, (1, 3, 1, 2)), (2, (1, 2, 2, 3))]:
        first = rows[:, :, field] - 8
        second = total - first
        feasible |= possible(4 - second, second - first, interval)
    return {"chords": len(packed), "constraint_patterns": len(rows),
            "feasible_chords": int(counts[feasible].sum())}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    result = check(parser.parse_args().directory)
    print(json.dumps(result, indent=2))
    if result["feasible_chords"]:
        raise SystemExit(2)
