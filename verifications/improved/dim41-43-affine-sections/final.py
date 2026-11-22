"""Count the three affine configurations from their family and point arrays."""
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def results():
    records = []
    for dimension in (41, 42, 43):
        with np.load(HERE / f"configuration_{dimension}.npz", allow_pickle=False) as data:
            sizes = data["populations"]
            assert sizes.shape == data["radii"].shape and np.all(sizes > 0)
            assert len(data["tiny_bodies"]) == len(data["tiny_heads"])
            total = int(sizes.sum()) + len(data["tiny_bodies"]) + len(data["pure"])
            records.append({"dimension": dimension, "total": total})
    return records


if __name__ == "__main__":
    print(json.dumps(results(), indent=2))
