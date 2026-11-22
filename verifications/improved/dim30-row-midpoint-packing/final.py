"""Read the complete coordinate counts for the repository registry."""

import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent


def results():
    records = json.loads((ROOT / "manifest.json").read_text())["records"]
    rows = []
    for record in records:
        with np.load(ROOT / record["coordinates"], allow_pickle=False) as archive:
            count, dimension = archive["points"].shape
        if (count, dimension) != (record["count"], record["dimension"]):
            raise ValueError("Certificate count or dimension differs from its manifest")
        rows.append({"dimension": dimension, "total": count})
    return rows


if __name__ == "__main__":
    for row in results():
        print(json.dumps(row))
