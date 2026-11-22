"""Pinned inputs for the discrete 27D searches."""
import hashlib
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OCTAD = np.array([0, 2, 4, 5, 6, 10, 11, 23])
COMPLEMENT = np.setdiff1d(np.arange(24), OCTAD)
RESIDUAL = np.r_[OCTAD, [24, 25, 26]]

SOURCES = {
    "d27_coordinates.npz": "ad94d28223a03a19907c5c2f15cea02ae2f499d48e7c4f180c7196ba6ceca708",
    "d27_certificate.npz": "8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94",
}
sys.path.insert(0, str(ROOT / "common"))
from leech import build as leech_shell


def checked_archive(name):
    path = HERE / "data" / name
    if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCES[name]:
        raise ValueError(f"Unexpected input digest for {path}")
    return np.load(path, allow_pickle=False)


def keys(rows):
    rows = np.ascontiguousarray(rows, dtype="<i2")
    return rows.view(np.dtype((np.void, 2 * rows.shape[1]))).ravel()


def load_points():
    with checked_archive("d27_coordinates.npz") as archive:
        points = archive["points"]
        if points.dtype.kind != "i" or points.shape != (210108, 27):
            raise ValueError("Invalid 27D coordinate matrix")
        if int(archive["squared_norm"]) != 128 or int(archive["inner_product_limit"]) != 64:
            raise ValueError("Invalid source normalization")
    return points


def verify_recipe(points, shell):
    """Check sources and all midpoint endpoint identities in repo ordering."""
    ordering = np.lexsort(shell.T[::-1])
    with checked_archive("d27_certificate.npz") as archive:
        additions = archive["added_integer"]
        rows = archive["midpoint_rows"]
        pairs = ordering[archive["midpoint_endpoints"]]
        deleted = ordering[archive["deleted_leech_indices"]]
    if len(rows) != 26496 or len(deleted) != 13248:
        raise ValueError("Unexpected midpoint population")
    if not np.array_equal(shell[pairs].sum(axis=1), 2 * additions[rows, :24]):
        raise ValueError("Incorrect midpoint endpoint identity")
    if not np.all(np.sum(shell[pairs[:, 0]] * shell[pairs[:, 1]], axis=1) == 64):
        raise ValueError("Incorrect Leech contact")
    if not np.array_equal(np.unique(pairs), np.sort(deleted)):
        raise ValueError("Incorrect declared Leech deletion set")
    retained = np.delete(shell, deleted, axis=0)
    retained = np.pad(retained, ((0, 0), (0, 3)))
    if not np.array_equal(np.sort(keys(np.r_[retained, additions])), np.sort(keys(points))):
        raise ValueError("The two source archives do not describe the same code")
    return {"leech_shell": len(shell), "source_size": len(points), "midpoint_rows": len(rows)}


if __name__ == "__main__":
    import json

    source = load_points()
    print(json.dumps(verify_recipe(source, 2 * leech_shell().astype(np.int16)), indent=2))
