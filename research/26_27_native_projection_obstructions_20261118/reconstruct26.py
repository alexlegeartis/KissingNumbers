"""Rebuild the alternative 26-dimensional configuration from the supplied exchange."""

import argparse
import hashlib
from pathlib import Path
import numpy as np

SOURCE = "0f0d9447269f6581057a4944308594a53d96b1f1310c80aedda20d7ecfc230f5"
TARGET = "56e4e1c18d328bfb7265e60fcb1bfd9c9b454e997090aa4ad36be2fb15c1b5b8"


def keys(a):
    a = np.ascontiguousarray(a, dtype="<i2")
    return a.view(np.dtype((np.void, 2 * a.shape[1]))).ravel()


def reconstruct(source, exchange, target):
    if hashlib.sha256(source.read_bytes()).hexdigest() != SOURCE:
        raise ValueError("Incorrect input archive")
    with np.load(source, allow_pickle=False) as z:
        old = z["integer_points"]
        n = z["fraction_numerators"]
        d = z["fraction_denominators"]
    with np.load(exchange, allow_pickle=False) as z:
        added, removed = z["inserted"], z["removed"]
    ko, kr, ka = map(keys, (old, removed, added))
    if len(np.unique(kr)) != 768 or len(np.unique(ka)) != 768:
        raise ValueError("Incorrect exchange cardinality")
    if np.any(np.isin(ka, ko)) or not np.all(np.isin(kr, ko)):
        raise ValueError("Exchange disagrees with source membership")
    out = np.vstack([old[~np.isin(ko, kr)], added])
    out = out[np.lexsort(out.T[::-1])].astype("<i2")
    if hashlib.sha256(out.tobytes()).hexdigest() != TARGET:
        raise ValueError("Incorrect integer coordinate digest")
    np.savez_compressed(
        target, integer_points=out, fraction_numerators=n, fraction_denominators=d
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["source", "exchange", "target"]:
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    reconstruct(args.source, args.exchange, args.target)
