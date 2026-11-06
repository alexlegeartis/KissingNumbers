#!/usr/bin/env python3
"""Independent exact sparse and dense/sparse audit of the 769436-point code.

This file deliberately imports no producer implementation. It consumes only
JSON certificates and source text data. All dot products are integer operations.
The separate dense/dense audit must establish dense orbit size and separation.
"""

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import time

import numpy as np


FROZEN_ORBIT_SHA256 = "e9f2dd59bf5169bf915be8e7d622a7ff14481471e2e0a96e0f9d90b46a5c6ad7"


def check(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def histogram(values):
    return {str(k): v for k, v in sorted(Counter(values).items())}


def signs(masks, width):
    return np.array([[1 - 2 * ((mask >> bit) & 1) for mask in masks]
                     for bit in range(width)], dtype=np.int64)


def squared_norms(vectors):
    return np.sum(vectors * vectors, axis=1, dtype=np.int64)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    start = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=True)
    orbit_path = args.source / "orbit_certificate.json"
    base_path = args.source / "base_certificate.json"
    check(sha256(orbit_path) == FROZEN_ORBIT_SHA256, "Frozen candidate SHA-256 mismatch")
    orbit = json.loads(orbit_path.read_text())
    base = json.loads(base_path.read_text())
    input_hashes = {str(orbit_path): sha256(orbit_path), str(base_path): sha256(base_path)}
    for filename, expected in base["source_sha256"].items():
        path = args.source / "source" / filename
        check(sha256(path) == expected, "Source SHA-256 mismatch: " + filename)
        input_hashes[str(path)] = sha256(path)

    n = orbit["dimension"]
    check(n == base["dimension"] == 39, "Dimension mismatch")
    lower = orbit["norm_lower"]
    upper = orbit["norm_upper"]
    check(lower == 99999996 and upper == 100000004, "Norm limits changed")
    reps_list = orbit["representatives"]
    check(len(reps_list) == 10256, "Representative count")
    check(all(len(row) == n and all(type(v) is int for v in row) for row in reps_list),
          "Representatives must be integer 39-tuples")
    reps = np.array(reps_list, dtype=np.int64)
    max_coordinate = max(abs(v) for row in reps_list for v in row)
    # Establish integer overflow safety before using NumPy arithmetic.
    check(n * max_coordinate ** 2 < 2 ** 63, "int64 squared norm overflow bound")
    check(8 * max_coordinate < 2 ** 63, "int64 middle product overflow bound")
    norms = squared_norms(reps)
    check(all(lower ** 2 <= int(q) <= upper ** 2 for q in norms), "Dense norm bounds")
    group = orbit["group_words"]
    check(len(group) == len(set(group)) == 32, "Sign-mask count/distinctness")
    check(all(type(g) is int and 0 <= g < 2 ** n for g in group), "Sign-mask range")
    check(0 in group, "Group identity missing")
    group_set = set(group)
    check(all(g ^ h in group_set for g in group for h in group), "XOR group closure")

    lines = (args.source / "source" / "supports_3424.txt").read_text().splitlines()
    check(len(lines) == 3424, "Source support count")
    check(all(len(line) == n and set(line) <= {"0", "1"} for line in lines),
          "Source bitstring encoding")
    original = [tuple(i for i, c in enumerate(line) if c == "1") for line in lines]
    check(all(len(s) == 8 for s in original), "Source support weight")
    check(len(set(original)) == 3424, "Duplicate source support")
    supports = list(original)
    replacement_records = []
    edited_indices = set()
    check(len(base["support_replacements"]) == 4, "Replacement count")
    for edit in base["support_replacements"]:
        i, old, new = edit["remove_index"], tuple(edit["remove"]), tuple(edit["add"])
        check(type(i) is int and 0 <= i < len(supports), "Replacement index")
        check(i not in edited_indices and supports[i] == old, "Replacement source mismatch")
        check(new == tuple(sorted(set(new))) and len(new) == 8
              and all(0 <= j < n for j in new), "Replacement support encoding")
        edited_indices.add(i)
        supports[i] = new
        replacement_records.append({"index": i, "source": list(old), "actual": list(new)})
    check(len(set(supports)) == 3424, "Duplicate actual support")
    check({i for i in range(3424) if supports[i] != original[i]} == edited_indices,
          "Undeclared support changes")

    even = tuple(mask for mask in range(256) if mask.bit_count() % 2 == 0)
    removed_by_support = {}
    check(len(orbit["partial_middle"]) == 0, "Partial-support record count")
    for item in orbit["partial_middle"]:
        i = item["support_index"]
        removed = item["removed_patterns"]
        check(i not in removed_by_support and type(i) is int and 0 <= i < 3424,
              "Partial-support index")
        check(all(type(p) is int and p in even for p in removed), "Deleted pattern parity/range")
        check(len(removed) == len(set(removed)), "Duplicate deleted pattern")
        removed_by_support[i] = set(removed)
    check(not removed_by_support,
          "Partial deletion location/count")
    check(supports[1149] == (1, 3, 14, 17, 30, 32, 33, 37), "Partial support coordinates")

    records = []
    retained = []
    reachable = []
    for i, support in enumerate(supports):
        patterns = tuple(p for p in even if p not in removed_by_support.get(i, set()))
        projected = sorted({sum(((g >> coordinate) & 1) << j
                                for j, coordinate in enumerate(support)) for g in group})
        # A literal union, not a parity shortcut: <g*x, sparse(p)> =
        # <x, sparse(p XOR projection(g))>.
        actual = tuple(sorted({p ^ q for p in patterns for q in projected}))
        retained.append(patterns)
        reachable.append(actual)
        records.append({"support_index": i, "coordinates": list(support),
                        "retained_masks": list(patterns), "projected_group_masks": projected,
                        "reachable_masks": list(actual),
                        "retained_parity_counts": histogram(p.bit_count() % 2 for p in patterns),
                        "reachable_parity_counts": histogram(p.bit_count() % 2 for p in actual)})
    check(len(retained[1149]) == 128 and len(reachable[1149]) == 256
          and len(records[1149]["projected_group_masks"]) == 32, "Partial literal-mask counts")
    middle_count = sum(map(len, retained))
    check(middle_count == 3424 * 128 == 438272, "Actual middle count")

    # Exact sparse geometry, starting with every pair of distinct supports.
    support_bits = [sum(1 << j for j in support) for support in supports]
    intersections = Counter()
    for i, mask in enumerate(support_bits):
        for other in support_bits[:i]:
            intersections[(mask & other).bit_count()] += 1
    check(max(intersections) <= 4, "Middle supports intersect in more than four coordinates")
    minimum_same_support_hamming = 8
    maximum_same_support_dot = -32
    # Enumerate every actual pattern class, rather than assuming parity.
    pattern_classes = sorted(set(retained))
    for patterns in pattern_classes:
        for j, p in enumerate(patterns):
            for q in patterns[:j]:
                distance = (p ^ q).bit_count()
                minimum_same_support_hamming = min(minimum_same_support_hamming, distance)
                maximum_same_support_dot = max(maximum_same_support_dot, 32 - 8 * distance)
    check(maximum_same_support_dot <= 16, "Middle same-support separation")

    star_specs = base["heavy_stars"]
    check(len(star_specs) == 4, "Heavy-star count")
    deleted_pairs = set()
    stars = []
    for spec in star_specs:
        center, leaves = spec["center"], spec["leaves"]
        check(type(center) is int and len(leaves) == len(set(leaves)) == 4
              and center not in leaves and all(type(j) is int and 0 <= j < n
                                              for j in [center] + leaves), "Star coordinates")
        for leaf in leaves:
            deleted_pairs.add(tuple(sorted((center, leaf))))
        for center_sign in (-1, 1):
            for mask in range(16):
                if mask.bit_count() % 2:
                    continue
                row = [0] * n
                row[center] = 4 * center_sign
                for j, leaf in enumerate(leaves):
                    row[leaf] = 2 * (1 - 2 * ((mask >> j) & 1))
                stars.append(row)
    check(len(deleted_pairs) == 14, "Deleted root-pair count")
    check(len(stars) == len({tuple(v) for v in stars}) == 64, "Star distinctness/count")
    stars = np.array(stars, dtype=np.int64)
    root_pairs = [pair for pair in itertools.combinations(range(n), 2) if pair not in deleted_pairs]
    roots = []
    for a, b in root_pairs:
        for sa, sb in itertools.product((-1, 1), repeat=2):
            row = [0] * n
            row[a], row[b] = 4 * sa, 4 * sb
            roots.append(row)
    check(len(roots) == len({tuple(v) for v in roots}) == 2908, "Root distinctness/count")
    roots = np.array(roots, dtype=np.int64)
    check(np.all(squared_norms(roots) == 32) and np.all(squared_norms(stars) == 32),
          "Sparse vector norm")
    # The actual magnitude profiles (eight 2s; two 4s; one 4 and four 2s)
    # differ, so different sparse layers cannot contain duplicate vectors.
    check(39 * 4 * 4 < 2 ** 15, "Sparse arithmetic overflow bound")

    def off_diagonal_max(vectors):
        result = -33
        # int64 throughout, excluding exactly each vector's self product.
        for lo in range(0, len(vectors), 128):
            product = vectors[lo:lo + 128] @ vectors.T
            for j in range(len(product)):
                product[j, lo + j] = -33
            result = max(result, int(product.max()))
        return result

    root_root_max = off_diagonal_max(roots)
    star_star_max = off_diagonal_max(stars)
    root_star_max = int((roots @ stars.T).max())
    check(max(root_root_max, star_star_max, root_star_max) <= 16, "Root/star sparse separation")

    middle_star_max = -33
    middle_star_witness = None
    middle_root_max = -33
    for i, support in enumerate(supports):
        p = signs(retained[i], 8).T * 2
        product = p @ stars[:, support].T
        local = int(product.max())
        if local > middle_star_max:
            middle_star_max = local
            j, k = np.unravel_index(int(product.argmax()), product.shape)
            middle_star_witness = {"support_index": i, "pattern": retained[i][j], "star_index": int(k)}
        # Each retained root pair includes all four signs, hence its maximum
        # against ANY middle sign pattern is exactly 8 times its overlap.
        support_set = set(support)
        middle_root_max = max(middle_root_max,
                              max(8 * ((a in support_set) + (b in support_set)) for a, b in root_pairs))
    check(max(middle_star_max, middle_root_max) <= 16, "Middle/root/star separation")

    # Dense/middle: enumerate every reachable actual 8-coordinate signed sum
    # for all 10256 representatives. No producer parity routine is used.
    # Each term has magnitude <= max_coordinate, so int64 exactness follows
    # from the separately checked 8*max_coordinate bound. Squaring below
    # is done in Python arbitrary-precision integers, without tolerance.
    middle_max = -1
    middle_witness = None
    middle_signed_sums_checked = 0
    per_support_dense_maxima = []
    for i, support in enumerate(supports):
        products = reps[:, support] @ signs(reachable[i], 8)
        local = int(products.max())
        check(local ** 2 <= 2 * lower ** 2, f"Dense/middle failure at support {i}")
        per_support_dense_maxima.append(local)
        middle_signed_sums_checked += products.size
        if local > middle_max:
            middle_max = local
            r, p = np.unravel_index(int(products.argmax()), products.shape)
            middle_witness = {"representative": int(r), "support_index": i,
                              "reachable_mask": reachable[i][p]}
        if (i + 1) % 500 == 0:
            print(f"Checked dense/middle supports {i + 1}/3424", flush=True)

    # Dense/root: complete signs on each retained pair make the root layer
    # invariant under every coordinate sign flip, not just this group.
    # Maximum divided by 4 is exactly |x_a|+|x_b|.
    absolute_reps = np.abs(reps)
    root_max = -1
    root_witness = None
    for a, b in root_pairs:
        values = absolute_reps[:, a] + absolute_reps[:, b]
        local = int(values.max())
        if local > root_max:
            root_max, root_witness = local, {"representative": int(values.argmax()), "pair": [a, b]}
    check(2 * root_max ** 2 <= lower ** 2, "Dense/root failure")

    # Dense/star: enumerate the literal union of all 64 stars transformed by
    # all 32 group masks. Divide sparse vectors by 2 to compare M^2<=2L^2.
    group_signs = np.array([[1 - 2 * ((g >> j) & 1) for j in range(n)] for g in group], dtype=np.int64)
    transformed_stars = sorted({tuple(int(v) for v in row)
                                for star in stars // 2 for row in group_signs * star})
    transformed_stars_array = np.array(transformed_stars, dtype=np.int64)
    check(all(sum(abs(v) for v in row) == 6 for row in transformed_stars), "Star absolute coefficient sum")
    check(6 * max_coordinate < 2 ** 63, "Star product overflow bound")
    star_products = reps @ transformed_stars_array.T
    star_max = int(star_products.max())
    r, p = np.unravel_index(int(star_products.argmax()), star_products.shape)
    star_witness = {"representative": int(r), "transformed_star_half": list(transformed_stars[p])}
    check(star_max ** 2 <= 2 * lower ** 2, "Dense/star failure")

    sparse_count = middle_count + len(roots) + len(stars)
    check(sparse_count == 441244, "Sparse count")
    # All integer dense representatives have nonzero coordinates. Coordinate
    # sign flips and radial normalization preserve that, while sparse vectors
    # have at most eight nonzero coordinates. Thus dense/sparse duplicates
    # are excluded independently of dense/dense distinctness.
    check(np.all(reps != 0), "Dense/sparse support distinction")
    implied_dense_count = len(reps) * len(group)
    check(implied_dense_count == 328192 and sparse_count + implied_dense_count == orbit["claimed_size"] == 769436,
          "Claimed total count reconciliation")

    patterns_path = args.output / "actual_sparse_patterns39.json"
    patterns_document = {
        "candidate_sha256": FROZEN_ORBIT_SHA256,
        "coordinate_convention": "Source bitstring character j is coordinate j. A mask bit j denotes a negative sign at sorted support coordinate j.",
        "partial_deletion_semantics": "No partial deletions. Every actual support retains all 128 even masks. Reachable masks are used only to reduce dense/sparse checking and are not additional sparse points.",
        "middle_supports": records,
        "deleted_root_pairs": [list(p) for p in sorted(deleted_pairs)],
        "root_vectors": roots.tolist(),
        "star_vectors": stars.tolist(),
        "group_transformed_star_vectors_divided_by_two": [list(v) for v in transformed_stars],
        "per_support_maximum_dense_signed_sum": per_support_dense_maxima,
    }
    patterns_path.write_text(json.dumps(patterns_document, separators=(",", ":")) + "\n")
    report = {
        "status": "PASS",
        "audit_scope": "Independent sparse/sparse geometry, actual support and sign-pattern reconstruction, deletion semantics, and every dense/sparse pair via exact sign-group reduction. Dense/dense separation and dense orbit distinctness are deliberately outside this file's scope.",
        "implementation": "Python 3 and NumPy integer matrix multiplication only; no producer Python implementation imported or executed.",
        "input_sha256": input_hashes,
        "implementation_sha256": sha256(Path(__file__)),
        "patterns_artifact": str(patterns_path),
        "patterns_artifact_sha256": sha256(patterns_path),
        "integer_arithmetic": {"type": "int64 dot products; arbitrary-precision Python squared comparisons",
                               "maximum_absolute_dense_coordinate": max_coordinate,
                               "maximum_possible_dense_squared_norm_sum": n * max_coordinate ** 2,
                               "maximum_possible_middle_signed_sum_absolute": 8 * max_coordinate,
                               "maximum_possible_star_half_dot_absolute": 6 * max_coordinate,
                               "float_arithmetic_used_for_acceptance": False,
                               "numerical_tolerance_used": False},
        "dense_norms": {"lower": lower, "upper": upper,
                        "minimum_squared_norm": int(norms.min()), "maximum_squared_norm": int(norms.max())},
        "support_replacements": replacement_records,
        "middle_patterns": {"support_count": len(supports),
                            "retained_pattern_count_distribution": histogram(map(len, retained)),
                            "reachable_pattern_count_distribution": histogram(map(len, reachable)),
                            "projected_group_size_distribution": histogram(len(r["projected_group_masks"]) for r in records),
                            "reachable_parity_distribution_by_support": dict(Counter(str(r["reachable_parity_counts"]) for r in records)),
                            "special_support_1149": records[1149],
                            "deleted_mask_count": sum(map(len, removed_by_support.values())),
                            "middle_points": middle_count},
        "sparse_geometry": {"all_squared_norms": 32,
                            "maximum_allowed_inner_product": 16,
                            "distinct_support_pairs_checked": sum(intersections.values()),
                            "distinct_support_intersection_histogram": dict(sorted(intersections.items())),
                            "distinct_support_middle_dot_upper_bound": 4 * max(intersections),
                            "same_support_pattern_classes_enumerated": len(pattern_classes),
                            "same_support_minimum_sign_hamming_distance": minimum_same_support_hamming,
                            "same_support_maximum_middle_dot": maximum_same_support_dot,
                            "root_root_maximum_dot": root_root_max,
                            "star_star_maximum_dot": star_star_max,
                            "root_star_maximum_dot": root_star_max,
                            "middle_root_maximum_dot": middle_root_max,
                            "middle_star_maximum_dot": middle_star_max,
                            "middle_star_witness": middle_star_witness,
                            "sparse_distinctness": "Unique supports, unique retained masks, unique roots/stars, and disjoint absolute coordinate profiles."},
        "dense_sparse": {
            "normalization_argument": "Normalize every dense representative and its sign orbit to norm sqrt(32). A sparse vector y has norm sqrt(32), so <sqrt(32)*x/||x||,y><=16 is equivalent to <x,y><=sqrt(8)*||x||. This audit uses the stronger uniform lower norm L.",
            "middle": {"literal_signed_sums_checked": middle_signed_sums_checked,
                       "maximum_signed_sum": middle_max, "witness": middle_witness,
                       "exact_acceptance_inequality": "M^2 <= 2 L^2",
                       "lhs": middle_max ** 2, "rhs": 2 * lower ** 2,
                       "nonnegative_integer_slack": 2 * lower ** 2 - middle_max ** 2},
            "root": {"representative_pair_checks": len(reps) * len(root_pairs),
                     "maximum_absolute_coordinate_pair_sum": root_max, "witness": root_witness,
                     "exact_acceptance_inequality": "2 M^2 <= L^2",
                     "lhs": 2 * root_max ** 2, "rhs": lower ** 2,
                     "nonnegative_integer_slack": lower ** 2 - 2 * root_max ** 2},
            "star": {"unique_group_transformed_star_vectors": len(transformed_stars),
                     "signed_products_checked": int(star_products.size),
                     "maximum_half_star_dot": star_max, "witness": star_witness,
                     "exact_acceptance_inequality": "M^2 <= 2 L^2",
                     "lhs": star_max ** 2, "rhs": 2 * lower ** 2,
                     "nonnegative_integer_slack": 2 * lower ** 2 - star_max ** 2},
            "dense_sparse_distinctness": "Every dense coordinate is nonzero, preserved by group action and normalization; every sparse vector has at least 31 zero coordinates."},
        "counts": {"middle": middle_count, "roots": len(roots), "stars": len(stars),
                   "sparse_total": sparse_count, "dense_representatives": len(reps), "group_order": len(group),
                   "dense_total_conditional_on_independent_dense_distinctness": implied_dense_count,
                   "total_conditional_on_independent_dense_distinctness": sparse_count + implied_dense_count},
        "elapsed_seconds": time.monotonic() - start,
    }
    report_path = args.output / "sparse_dense_sparse_audit39.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "report": str(report_path),
                      "elapsed_seconds": report["elapsed_seconds"],
                      "counts": report["counts"], "dense_middle_max": middle_max,
                      "dense_root_max": root_max, "dense_star_max": star_max}, indent=2), flush=True)


if __name__ == "__main__":
    main()
