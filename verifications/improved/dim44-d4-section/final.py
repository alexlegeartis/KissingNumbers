"""Read the 44D population ledger. The coordinate and moment replay is verify.py."""
from pathlib import Path

import numpy as np


def results():
    with np.load(Path(__file__).with_name("families.npz"), allow_pickle=False) as data:
        sizes = data["populations"]
        assert sizes.shape == (17,) and np.all(sizes > 0)
        return [{"dimension": 44, "total": int(sizes.sum())}]


if __name__ == "__main__":
    print(results()[0])
