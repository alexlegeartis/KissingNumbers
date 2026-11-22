"""Count the P48q affine witness from its ledger and explicit short layers."""

import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def results():
    data = json.loads((HERE / "certificate.json").read_text())
    counts = data["counts"]
    starting = {"54": 1_084_590, "36": 40_590, "30": 1134, "6": 0}
    for norm in ("36", "30", "6"):
        assert counts[norm] == starting[norm] + data["refill_counts"][norm] - data["removal_counts"][norm]
    assert counts["54"] == starting["54"]
    with np.load(HERE / "points.npz", allow_pickle=False) as points:
        assert all(len(points[f"norm{norm}"]) == counts[str(norm)] for norm in (36, 30, 6))
    total = sum(counts.values())
    assert total == data["total"] and data["dimension"] == 40
    return [{"dimension": 40, "total": total}]


if __name__ == "__main__":
    print(json.dumps(results()[0]))
