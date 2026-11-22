"""Enumerate the pure Barnes–Wall midpoint refills with an optional new fibre class."""
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np

from source import (COMPLEMENT, HERE, OCTAD, RESIDUAL, checked_archive,
                    keys, leech_shell, load_points, verify_recipe)

TAILS = np.array([t for t in product((-4, 0, 4), repeat=3)
                  if sum(v != 0 for v in t) == 2 and next(v for v in t if v) > 0],
                 dtype=np.int16)


def overlaps(left, right):
    """Identify points on the right that violate a constraint from the left."""
    left = left.astype(np.float32)
    result = []
    for start in range(0, len(right), 2048):
        result.append((left @ right[start:start + 2048].astype(np.float32).T > 64).any(0))
    return np.concatenate(result)


def tilted_section(shell, use_anchor, head_frames=None):
    with checked_archive("d27_certificate.npz") as archive:
        source_additions = archive["added_integer"]
        midpoint_rows = archive["midpoint_rows"]
        ordering = np.lexsort(shell.T[::-1])
        contact_endpoints = ordering[archive["midpoint_endpoints"]]

    pure = np.zeros(len(source_additions), dtype=bool)
    pure[midpoint_rows] = ~np.any(source_additions[midpoint_rows][:, OCTAD], axis=1)
    frozen = source_additions[~pure]
    midpoints = source_additions[midpoint_rows]
    pure_parents = contact_endpoints[pure[midpoint_rows]]

    if use_anchor:
        with np.load(HERE / "data" / "refill_anchor.npz", allow_pickle=False) as archive:
            projections, residuals = archive["projections"], archive["residuals"]
        new = np.empty((len(projections), 2 * len(residuals), 27), dtype=np.int16)
        new[:, :, COMPLEMENT] = projections[:, None, :]
        reflected = residuals.copy()
        reflected[:, 8:] *= -1
        new[:, :, RESIDUAL] = np.r_[residuals, reflected][None]
        new = new.reshape(-1, 27)
        assert new.shape == (256, 27)
        assert np.all(np.sum(new * new, axis=1) == 128) and len(np.unique(keys(new))) == len(new)
        pairwise = new @ new.T
        np.fill_diagonal(pairwise, 0)
        assert pairwise.max() <= 64
        bad = overlaps(new, frozen)
        # The added class may displace row midpoints. All other non-BW
        # source additions have to remain compatible.
        assert np.all(np.sum(frozen[bad][:, OCTAD] ** 2, axis=1) == 64)
        assert np.all(np.sum(frozen[bad][:, COMPLEMENT] ** 2, axis=1) == 32)
        assert np.all(np.sum(frozen[bad][:, 24:] ** 2, axis=1) == 32)
        frozen = np.r_[frozen[~bad], new]
        new_deletions = np.flatnonzero(overlaps(new[:, :24], shell))
        assert len(new_deletions) == 128
    else:
        new_deletions = np.empty(0, dtype=np.int64)
    if head_frames is not None:
        from fibres import native_fibres
        from head_frames import make_heads

        frames, _ = native_fibres(shell)
        old_heads = ((frozen[:, COMPLEMENT] ** 2).sum(1) == 64) & (
            (frozen[:, 24:] ** 2).sum(1) == 64
        )
        assert old_heads.sum() == 192
        other = frozen[~old_heads]
        new_heads = np.concatenate([make_heads(frames[c], a) for a, c in enumerate(head_frames)])
        gram = new_heads @ new_heads.T
        np.fill_diagonal(gram, 0)
        if gram.max() > 64 or overlaps(new_heads, other).any() or overlaps(new_heads[:, :24], shell).any():
            raise ValueError("The requested heads conflict with the fixed construction")
        frozen = np.concatenate([other, new_heads])
    still_present = np.isin(keys(midpoints), keys(frozen))
    deleted = np.unique(np.r_[contact_endpoints[still_present].ravel(), new_deletions])
    assert not np.isin(pure_parents, deleted).any()
    surviving_bw = source_additions[pure]
    if use_anchor or head_frames is not None:
        surviving_bw = surviving_bw[~overlaps(frozen, surviving_bw)]
    return frozen, deleted, surviving_bw


def contact_blocks(shell):
    bw_ids = np.flatnonzero(~np.any(shell[:, OCTAD], axis=1))
    assert len(bw_ids) == 4320
    bw = shell[bw_ids][:, COMPLEMENT].astype(np.float32)
    for i in range(0, len(bw), 128):
        first, second = np.nonzero(bw[i:i + 128] @ bw[i:].T == 64)
        forward = first < second
        yield np.c_[bw_ids[i + first[forward]], bw_ids[i + second[forward]]]


def refill_options(shell, frozen):
    tail_values, indices = np.unique(frozen[:, 24:], axis=0, return_inverse=True)
    groups = []
    for k in np.argsort(-np.bincount(indices)):
        bodies = np.unique(frozen[indices == k][:, COMPLEMENT], axis=0).astype(np.float32).T
        limit = 64 - np.abs(TAILS @ tail_values[k])
        groups.append((bodies, limit))
    total_contacts = 0
    result = []
    for edges in contact_blocks(shell):
        total_contacts += len(edges)
        bodies = (shell[edges[:, 0]][:, COMPLEMENT] + shell[edges[:, 1]][:, COMPLEMENT]) // 2
        allowed = np.ones((len(bodies), len(TAILS)), dtype=bool)
        for constraint, threshold in groups:
            test = np.flatnonzero(allowed.any(1))
            if not len(test):
                break
            dots = (bodies[test].astype(np.float32) @ constraint).max(1)
            allowed[test] &= dots[:, None] <= threshold
        idx, tail_index = np.where(allowed)
        result.append((bodies[idx], edges[idx], tail_index))
    assert total_contacts == 604800
    return tuple(map(np.concatenate, zip(*result)))


def run(use_anchor, output=None, head_frames=None):
    shell = 2 * leech_shell().astype(np.int16)
    verify_recipe(load_points(), shell)
    fixed, deleted, survivors = tilted_section(shell, use_anchor, head_frames)
    bodies, parents, colors = refill_options(shell, fixed)
    allowed = np.c_[bodies, TAILS[colors]]
    incumbent = np.isin(keys(allowed), keys(survivors[:, np.r_[COMPLEMENT, [24, 25, 26]]]))
    assert 2 * incumbent.sum() == len(survivors)
    base = len(shell) - len(deleted) + len(fixed)
    score = 2 * int(incumbent.sum()) - len(np.unique(parents[incumbent]))
    if output:
        np.savez_compressed(output, bodies=bodies, parents=parents, colors=colors, tails=TAILS,
                            incumbent=incumbent, fixed=fixed, deleted=deleted, shell=shell)
    return {
        "extra_fibre_class": use_anchor,
        "contact_pairs": 604800,
        "admissible_tail_pairs": len(bodies),
        "eligible_endpoints": len(np.unique(parents)),
        "retained_source_tail_pairs": int(incumbent.sum()),
        "fixed_points_before_refill": base,
        "construction_with_retained_source_refill": base + score,
        "refill_gain_needed_to_beat_210108": 210109 - base,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extra-class", action="store_true")
    parser.add_argument("--head-frames", nargs=3, type=int, choices=range(135),
                        metavar=("FIRST", "SECOND", "THIRD"))
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.extra_class, args.out, args.head_frames), indent=2))
