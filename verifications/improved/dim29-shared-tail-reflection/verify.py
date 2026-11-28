#!/usr/bin/env python3
"""Verify the 232,242-point 29D tail-antipode splice from its pinned witness."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE_SHA256 = "81d5f3043438cbf3f8f3521033ecaa5de44a50e774b4610caee8ce54d2b234c6"
RESULT_SHA256 = "ac8fb373b0321a42ce71ffbeef0dde75b9e74a3e4c8ebb807d24ba20d910d848"

def point_hash(rows):
    order = np.lexsort(rows.T[::-1])
    return hashlib.sha256(rows[order].astype("<i2").tobytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path, help="232236-point NPZ written by the pinned Kempe verifier")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with np.load(args.baseline, allow_pickle=False) as f:
        raw = f["points"]
        assert raw.dtype.kind in "iu" and raw.shape == (232236, 29)
        assert np.all(np.abs(raw.astype(np.int64)) <= 128)
        base = raw.astype(np.int64)
    assert point_hash(base) == BASE_SHA256
    edits = json.loads((HERE / "edits.json").read_text())
    add = np.asarray(edits["add"], dtype=np.int64)
    delete = np.asarray(edits["remove"], dtype=np.int64)
    assert add.shape == (10, 29) and delete.shape == (4, 29)
    mask = np.ones(len(base), dtype=bool)
    for row in delete:
        ix = np.flatnonzero((base == row).all(axis=1))
        assert len(ix) == 1 and mask[ix[0]]
        mask[ix[0]] = False
    for row in add:
        mate = row.copy()
        mate[24:] *= -1
        assert (base == mate).all(axis=1).sum() == 1
    keep = base[mask]
    assert len(keep) == 232232
    assert np.all(np.einsum("ij,ij->i", add, add) == 128)
    g_old = add @ keep.T
    g_new = add @ add.T
    np.fill_diagonal(g_new, -128)
    assert g_old.max() <= 64 and g_new.max() <= 64
    final = np.concatenate((keep, add))
    assert final.shape == (232242, 29)
    assert point_hash(final) == RESULT_SHA256
    print(json.dumps({
        "points": len(final), "inherited_rows": len(keep), "removed": len(delete),
        "added": len(add), "new_vs_retained_max": int(g_old.max()),
        "new_vs_new_max": int(g_new.max()), "sorted_int16_rows_sha256": RESULT_SHA256,
    }, indent=2))
    if args.output:
        order = np.lexsort(final.T[::-1])
        np.savez_compressed(args.output, points=final[order].astype("<i2"),
                            squared_norm=np.int64(128), inner_product_limit=np.int64(64))

if __name__ == "__main__":
    main()
