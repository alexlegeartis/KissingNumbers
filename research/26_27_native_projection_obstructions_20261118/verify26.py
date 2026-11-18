"""Verify a 26-dimensional certificate with integral and rational coordinates."""

import argparse
import json
import math
from pathlib import Path
import numpy as np
from verify_exact import verify


def verify_mixed(path):
    with np.load(path, allow_pickle=False) as z:
        points = z["integer_points"]
        numerators = z["fraction_numerators"]
        denominators = z["fraction_denominators"]
    geometry = verify(points)
    dimension = geometry["dimension"]
    if dimension != 26 or numerators.shape != (len(denominators), dimension):
        raise ValueError("Incorrect shape")
    if numerators.dtype.kind != "i" or denominators.dtype.kind != "i":
        raise ValueError("Rational coordinates must be supplied as integers")
    core = points.astype(np.int64)
    magnitude = max(abs(int(core.min())), abs(int(core.max())))
    for i, (row, denominator) in enumerate(zip(numerators, denominators)):
        n = list(map(int, row))
        d = int(denominator)
        if d <= 0 or sum(v * v for v in n) != 128 * d * d:
            raise ValueError("Rational norm or denominator failure")
        if dimension * magnitude * max(map(abs, n)) >= 2**63:
            raise ValueError("Rational/integer product exceeds this kernel's range")
        if int((core @ row.astype(np.int64)).max()) > 64 * d:
            raise ValueError("Rational/integer conflict")
        for j in range(i):
            dot = sum(a * int(b) for a, b in zip(n, numerators[j]))
            if dot > 64 * d * int(denominators[j]):
                raise ValueError("Rational/rational conflict")
    k = len(denominators)
    n = len(points) + k
    pairs = geometry["unordered_pairs"] + k * len(points) + math.comb(k, 2)
    if pairs != math.comb(n, 2):
        raise ValueError("Incomplete coverage")
    return {
        "dimension": dimension,
        "points": n,
        "rational_points": k,
        "integer_geometry": geometry,
        "unordered_pairs": pairs,
        "status": "PASS",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coordinates", type=Path)
    print(json.dumps(verify_mixed(parser.parse_args().coordinates), indent=2))
