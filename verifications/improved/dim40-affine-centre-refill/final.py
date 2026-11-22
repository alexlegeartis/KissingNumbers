"""Count the witness described by the affine-centre certificate."""

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def results():
    data = json.loads((HERE / "certificate.json").read_text())
    with np.load(HERE / "auxiliary.npz", allow_pickle=False) as source:
        medium = len(source["norm36"])
        tiny = len(source["norm30"])
    total = 1084590 + medium + tiny + len(data["refill"]) + len(data["norm36_refill"])
    assert total == data["total"]
    return [{"dimension": 40, "total": total}]


if __name__ == "__main__":
    for row in results():
        print(json.dumps(row))
