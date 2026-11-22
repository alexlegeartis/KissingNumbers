"""Enumerate low-cost lifted contact midpoints on a chosen native Leech frame."""

import argparse
import json
from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "common"))
from leech import build

OCTAD = np.array([0, 2, 4, 5, 6, 10, 11, 23])
COMPLEMENT = np.setdiff1d(np.arange(24), OCTAD)


def keys(rows):
    return np.ascontiguousarray(rows, dtype="<i2").view(f"V{2 * rows.shape[1]}").ravel()


def tail_roots():
    rows = []
    for pair in combinations(range(6), 2):
        for signs in product((-4, 4), repeat=2):
            row = np.zeros(6, dtype=np.int16)
            row[list(pair)] = signs
            rows.append(row)
    return np.stack(rows)


def frames(shell):
    # Each of 4,320 complement projections determines a signed octad frame.
    native = shell[(shell[:, OCTAD] ** 2).sum(axis=1) == 64]
    fibre = {}
    for root in native:
        fibre.setdefault(tuple(root[COMPLEMENT]), []).append(tuple(root[OCTAD]))
    classes = {}
    for complement, octad in fibre.items():
        classes.setdefault(tuple(sorted(octad)), []).append(complement)
    return [(np.array(a, dtype=np.int16), np.array(sorted(w), dtype=np.int16))
            for a, w in sorted(classes.items())]


def enumerate_family(fixed, moving, fixed_axes, moving_axes, lifted, equator, tails, cap):
    tail_types, type_ids = np.unique(lifted[:, 24:], axis=0, return_inverse=True)
    tail_products = tails @ tail_types.T
    maximum_tail = tail_products.max(axis=0)[type_ids]
    frame_products = moving @ lifted[:, moving_axes].T
    left, right = np.where(np.triu(moving @ moving.T == 0, 1))
    result = []
    for anchor in fixed:
        constant = anchor @ lifted[:, fixed_axes].T
        endpoints = np.zeros((len(moving), 24), dtype=np.int16)
        endpoints[:, fixed_axes] = anchor
        endpoints[:, moving_axes] = moving
        retained = np.isin(keys(endpoints), equator)
        for i, j in zip(left, right):
            # Only the two endpoints can block the midpoint in the Leech equator.
            root_cost = int(retained[i]) + int(retained[j])
            if root_cost > cap:
                continue
            body_dot = constant + (frame_products[i] + frame_products[j]) // 2
            visible = body_dot + maximum_tail > 64
            violations = tail_products[:, type_ids[visible]] + body_dot[visible] > 64
            keep = root_cost + np.count_nonzero(violations, axis=1) <= cap
            if np.any(keep):
                points = np.zeros((np.count_nonzero(keep), 30), dtype=np.int16)
                points[:, fixed_axes] = anchor
                points[:, moving_axes] = (moving[i] + moving[j]) // 2
                points[:, 24:] = tails[keep]
                result.append(points)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frame", type=int, choices=range(135), metavar="FRAME")
    parser.add_argument("output", type=Path)
    parser.add_argument("--family", choices=("row", "column", "both"), default="both")
    parser.add_argument("--max-blockers", type=int, default=6)
    args = parser.parse_args()
    record, = json.loads((HERE / "manifest.json").read_text())["records"]
    with np.load(REPO / record["source"]["path"], allow_pickle=False) as archive:
        source = archive["points"].astype(np.int16)
    shell = 2 * build().astype(np.int16)
    has_tail = np.any(source[:, 24:] != 0, axis=1)
    equator = keys(source[~has_tail, :24])
    if not np.isin(equator, keys(shell)).all():
        raise ValueError("Source equator contains points outside the Leech shell")
    octad, complement = frames(shell)[args.frame]
    families = []
    if args.family in ("row", "both"):
        families.append((octad, complement, OCTAD, COMPLEMENT))
    if args.family in ("column", "both"):
        families.append((complement, octad, COMPLEMENT, OCTAD))
    kept = []
    tails = tail_roots()
    lifted = source[has_tail]
    for family in families:
        kept.extend(enumerate_family(*family, lifted, equator, tails, args.max_blockers))
    candidates = np.concatenate(kept) if kept else np.empty((0, 30), dtype=np.int16)
    candidates = candidates[~np.isin(keys(candidates), keys(source))]
    np.savez_compressed(args.output, candidates=candidates)
    print(json.dumps(dict(frame=args.frame, family=args.family, candidates=len(candidates))))


if __name__ == "__main__":
    main()
