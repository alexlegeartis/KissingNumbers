"""Recompute the 36D bound from the pinned predecessor and phase delta."""

from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def results():
    cert = json.loads((HERE / "certificate.json").read_bytes())
    source = (HERE.parent / "dim36-ct-nonlocal" / "certificate.json").read_bytes()
    assert sha256(source).hexdigest() == cert["source_sha256"]
    parent = json.loads(source)
    assert cert["dimension"] == parent["dimension"] == 36
    added = len({tuple(row) for row in cert["added_points"]})
    removed = len({tuple(row) for row in cert["removed_points"]})
    assert added == len(cert["added_points"]) and removed == len(cert["removed_points"])
    total = parent["points"] + added - removed
    assert total == cert["points"] > parent["points"]
    return [{"dimension": 36, "total": total}]


if __name__ == "__main__":
    for row in results():
        print(json.dumps(row))
