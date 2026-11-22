"""Read the registered 33D witness count."""

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def results():
    info = json.loads((HERE / "manifest.json").read_text())
    with np.load(HERE / "witness.npz", allow_pickle=False) as witness:
        sparse, dimension = witness["sparse"].shape
        count = (1 << len(witness["dense_generator"])) + 128 * len(witness["supports"]) + sparse
    if (dimension, count) != (info["dimension"], info["count"]):
        raise ValueError("Manifest does not match the witness")
    return [{"dimension": dimension, "total": count}]


if __name__ == "__main__":
    print(json.dumps(results()[0]))
