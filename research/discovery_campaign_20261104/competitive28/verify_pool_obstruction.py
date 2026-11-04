"""Independent direct-integer Hall witness for an explicit candidate pool.

Only one forbidden baseline neighbour per candidate is needed: if these
neighbours are distinct, any subset of k added points forces at least k
deletions. Candidate/candidate conflicts can only make the family weaker.
No source PASS receipt, optimizer status, or full blocker table is trusted.
"""
from pathlib import Path
import hashlib
import json
import sys
from math import isqrt

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "verifications/improved/dim27-31-joint-refills/data/coordinates28_221988.npz"
SOURCE_SHA256 = "509b7894dd8534946a81d8bfd354c9ef947e4fe5c82482796aac76799169f129"


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def check(source, candidates, matching, norm):
    require(source.dtype.kind in 'iu' and candidates.dtype.kind in 'iu' and matching.dtype.kind in 'iu', "Integer witness arrays")
    require(source.ndim == candidates.ndim == 2 and source.shape[1] == candidates.shape[1], "Matching coordinate widths")
    require(matching.shape == (len(candidates),), "One witness owner per candidate")
    require(np.all(matching >= 0) and np.all(matching < len(source)), "Baseline owner index range")
    require(len(np.unique(matching)) == len(candidates), "Owner map must be injective")
    require(len(np.unique(candidates, axis=0)) == len(candidates), "Distinct candidate points")
    limit = isqrt((2**62 - 1) // source.shape[1])
    # Bound raw arrays before abs/casts: abs(minimum signed int64) overflows.
    require(np.all(source >= -limit) and np.all(source <= limit)
            and np.all(candidates >= -limit) and np.all(candidates <= limit), "Raw coordinate magnitude bound")
    bound = source.shape[1] * int(max(abs(source).max(initial=0), abs(candidates).max(initial=0)))**2
    require(bound < 2**62, "Proven int64 arithmetic bound")
    source_i, candidate_i = source.astype(np.int64), candidates.astype(np.int64)
    require(np.all(np.sum(source_i**2, axis=1) == norm), "Baseline norms")
    require(np.all(np.sum(candidate_i**2, axis=1) == norm), "Candidate norms")
    original = {row.tobytes() for row in source.astype(np.int16)}
    require(all(row.tobytes() not in original for row in candidates.astype(np.int16)), "Candidates are genuine additions")
    dots = np.sum(candidate_i * source_i[matching], axis=1)
    require(np.all(2 * dots > norm), "Every assigned owner must be an actual forbidden neighbour")
    return {"status": "PASS", "baseline_points": len(source), "distinct_candidates": len(candidates),
            "distinct_forbidden_owners": len(matching), "squared_norm": norm,
            "minimum_witness_dot": int(dots.min(initial=norm)),
            "maximum_possible_gain": 0,
            "proof": "Each distinct added point forces deletion of its assigned distinct incompatible baseline point. Every candidate subset costs at least its cardinality; all candidate/candidate constraints may be ignored.",
            "arithmetic": "direct int64 coordinate sums with proved overflow bound; no BLAS, matching optimizer, saved graph or saved PASS used",
            "scope": "all subsets of the explicit finite candidate pool with arbitrary baseline deletions; no bound on ungenerated candidates"}


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name('global_blocker_witness.npz')
    require(hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256, "Pinned baseline coordinate file hash")
    with np.load(SOURCE, allow_pickle=False) as archive:
        source = archive['points'].copy()
    with np.load(path, allow_pickle=False) as archive:
        candidates, matching = archive['points'].copy(), archive['matching'].copy()
    report = check(source, candidates, matching, 128)
    report.update(source_sha256=SOURCE_SHA256, witness_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
