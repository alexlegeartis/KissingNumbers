"""Rebuild the finite sparse conflict graph in METIS format."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
from tempfile import TemporaryFile

import numpy as np

HERE = Path(__file__).resolve().parent


def write_graph(points, destination):
    if points.dtype.kind not in "iu" or np.any(np.abs(points.astype(np.int64)) > 5):
        raise ValueError("The graph expects small integer coordinates")
    if points.ndim != 2 or points.shape[1] != 33:
        raise ValueError("The graph expects 33 coordinates")
    rows = points.astype(np.float32)
    # Each matrix entry is an integer and each partial sum has magnitude at most 825.
    edges = 0
    with TemporaryFile(mode="w+t") as body:
        for start in range(0, len(rows), 180):
            conflict = rows[start:start + 180] @ rows.T > 16
            local = np.arange(len(conflict))
            conflict[local, start + local] = False
            for row in conflict:
                neighbors = np.flatnonzero(row) + 1
                body.write(" ".join(map(str, neighbors)) + "\n")
                edges += len(neighbors)
        if edges % 2:
            raise ValueError("Expected a symmetric graph")
        body.seek(0)
        with destination.open("w") as output:
            output.write(f"{len(rows)} {edges // 2}\n")
            shutil.copyfileobj(body, output)
    return edges // 2


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("output", type=Path)
    cli.add_argument("--stage", choices=("initial", "repeat"), default="initial")
    args = cli.parse_args()
    record = json.loads((HERE / "manifest.json").read_text())["search_instances"][args.stage]
    path = HERE / record["path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
        raise ValueError("Search input digest mismatch")
    with np.load(path, allow_pickle=False) as archive:
        points = archive["v"]
    print(f"{len(points):,} vertices, {write_graph(points, args.output):,} edges")
