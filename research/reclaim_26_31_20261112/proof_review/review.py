#!/usr/bin/env python3
"""Read-only independent review replay for the November sextet certificates.

Checks upstream source pins, replays the direct-int64 delta implementation,
and compares the complete geometry verifier with a separate brute-force
integer predicate on nine small fixtures. This is not a second large all-pairs
verification. It does not alter the certificate package or its manifest.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np


REPO = Path(__file__).resolve().parents[3]


def require(condition, label):
    if not condition:
        raise ValueError(label)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def brute_force_valid(points):
    """Direct integer predicate, independent of the verifier's line grouping."""
    rows = [tuple(int(c) for c in row) for row in points]
    if len(set(rows)) != len(rows):
        return False
    if any(sum(c*c for c in row) != 128 for row in rows):
        return False
    for i, row in enumerate(rows):
        for other in rows[i+1:]:
            if sum(a*b for a, b in zip(row, other)) > 64:
                return False
    return True


def fixtures():
    roots = []
    for i, j in itertools.combinations(range(5), 2):
        for a, b in itertools.product((-8, 8), repeat=2):
            row = np.zeros(28, dtype=np.int16)
            row[i], row[j] = a, b
            roots.append(row)
    roots = np.asarray(roots)
    positive = np.asarray([r for r in roots if r[np.flatnonzero(r)[0]] > 0])
    mixed = np.delete(roots, np.arange(0, len(roots), 3), axis=0)
    a = np.zeros(28, dtype=np.int16)
    a[:2] = 8
    b = np.zeros(28, dtype=np.int16)
    b[:5] = (8, 4, 4, 4, 4)
    wrong = a.copy()
    wrong[0] = 9
    return {
        "valid_all_paired": roots,
        "valid_all_unpaired": positive,
        "valid_mixed": mixed,
        "valid_negative_large_product": np.asarray([a, -b]),
        "wrong_norm": np.asarray([wrong, -a]),
        "duplicate": np.asarray([a, a]),
        "illegal_unpaired": np.asarray([a, b]),
        "illegal_paired_unpaired": np.asarray([a, -a, b]),
        "illegal_paired_paired": np.asarray([a, -a, b, -b]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--package", type=Path,
        default=REPO/"verifications/improved/dim28-29-november-sextets")
    parser.add_argument(
        "--source-dir", type=Path,
        default=REPO.parent/"PackingStar/25D-31D/new_bounds_20261111",
        help="Directory containing the pinned D28 and D29 PackingStar NPZ files")
    parser.add_argument("--output", type=Path,
                        help="Optional JSON output; no file is written if omitted")
    args = parser.parse_args()
    package = args.package.resolve()
    geometry = load_module("review_geometry", package/"verify_geometry.py")
    delta = load_module("review_delta", package/"verify_delta.py")
    provenance = json.loads((package/"data/provenance.json").read_text())
    require(provenance["packingstar_commit"] ==
            "2f3453b79268d55ff7cac3e2f3a3b9d14ef590c6", "Expected source commit")
    source_checks = []
    for record in provenance["records"]:
        dim, count = record["dimension"], record["packingstar"]
        path = args.source_dir/f"D{dim}_{count}_coordinates.npz"
        require(digest(path) == record["source_file_sha256"], "Upstream file hash")
        with np.load(path, allow_pickle=False) as z:
            require(int(z["coordinate_denominator"]) == 1, "Source denominator")
            require(int(z["squared_norm"]) == 32, "Source squared norm")
            require(int(z["inner_product_limit"]) == 16, "Source product limit")
            require(z["points"].dtype.kind in "iu", "Integer source coordinates")
            source = z["points"].astype(np.int64)*2
        require(source.shape == (count, dim), "Source shape")
        require(np.all(np.sum(source*source, axis=1) == 128), "Source norms")
        # This normalization is independent of the delta replay's Python set
        # reconstruction and sorted(tuple(...)) implementation.
        ordered = np.asarray(source[np.lexsort(source.T[::-1])], dtype="<i2")
        require(hashlib.sha256(ordered.tobytes()).hexdigest() ==
                record["source_canonical_sha256"], "Normalized source set hash")
        replay = delta.replay(record)
        require(replay == record["direct_integer_replay"], "Recorded delta replay")
        source_checks.append({
            "dimension": dim,
            "source_file_hash": "PASS",
            "normalized_source_set_hash": "PASS",
            "delta": replay,
        })

    truth = {}
    checks = 0
    precisions = (32, 64)
    block_sizes = (1, 2, 3, 7, 64)
    with tempfile.TemporaryDirectory(prefix="sextet-proof-review-") as temporary:
        for name, points in fixtures().items():
            expected = brute_force_valid(points)
            truth[name] = expected
            path = Path(temporary)/(name+".npz")
            np.savez(path, points=points, squared_norm=128,
                     inner_product_limit=64, coordinate_denominator=2)
            for precision, block_size in itertools.product(precisions, block_sizes):
                try:
                    result = geometry.verify(path, precision, block_size)
                    accepted = result["status"] == "PASS"
                except ValueError:
                    accepted = False
                require(accepted == expected,
                        f"Fixture {name}, binary{precision}, block {block_size}")
                checks += 1
    require(checks == 90, "Expected fixture comparison count")
    report = {
        "status": "PASS",
        "scope": "Independent source pins and direct-int64 delta replay; small-fixture coverage audit",
        "large_configuration_all_pairs_repeated": False,
        "reviewed_implementation_sha256": {
            name: digest(package/name)
            for name in ("verify_geometry.py", "verify_delta.py")
        },
        "packingstar_commit": provenance["packingstar_commit"],
        "sources": source_checks,
        "independent_fixture_comparisons": checks,
        "precisions": list(precisions),
        "block_sizes": list(block_sizes),
        "fixtures": truth,
    }
    output = json.dumps(report, indent=2)+"\n"
    if args.output:
        args.output.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
