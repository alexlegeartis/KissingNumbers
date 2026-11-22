"""Prepare the exact 11D kernel inequalities and incumbent for kernel_grid.cpp."""
import argparse
from pathlib import Path

import numpy as np

from source import COMPLEMENT, RESIDUAL, keys, load_points


def export(path):
    points = load_points()
    outside = np.any(points[:, COMPLEMENT], axis=1)
    normals = np.unique(points[outside][:, RESIDUAL], axis=0)
    incumbent = points[~outside][:, RESIDUAL]
    for rows in (normals, incumbent):
        if not np.array_equal(np.sort(keys(rows)), np.sort(keys(-rows))):
            raise ValueError("This search requires a centrally symmetric projection set")
    assert (len(normals), len(incumbent)) == (2895, 476)
    pos = incumbent[np.arange(len(incumbent)), np.argmax(incumbent != 0, axis=1)] > 0
    with path.open("w") as output:
        output.write(f"{len(normals)} {pos.sum()}\n")
        np.savetxt(output, np.r_[normals, incumbent[pos]], fmt="%d")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    export(parser.parse_args().output)
