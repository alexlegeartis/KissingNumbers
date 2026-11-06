#!/usr/bin/env python3
"""Resolve reported extreme products back to actual original point pairs.

Uses pure Python integer arithmetic, independently of NumPy dot products.
Run with --size 769388 or --size 769436.
"""
import argparse
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--size", type=int, choices=(769388, 769436), default=769436)
parser.add_argument('--directory',type=Path,default=Path(__file__).resolve().parent)
args = parser.parse_args()
directory = args.directory
report = json.loads((directory / "sparse_dense_sparse_audit39.json").read_text())
patterns_path = directory / "actual_sparse_patterns39.json"
check = lambda condition: condition or (_ for _ in ()).throw(ValueError("Witness verification failed"))
check(hashlib.sha256(patterns_path.read_bytes()).hexdigest() == report["patterns_artifact_sha256"])
patterns = json.loads(patterns_path.read_text())
orbit_path = next(Path(p) for p in report["input_sha256"] if p.endswith("/orbit_certificate.json"))
check(hashlib.sha256(orbit_path.read_bytes()).hexdigest() == report["input_sha256"][str(orbit_path)])
orbit = json.loads(orbit_path.read_text())
rows = orbit["representatives"]
group = orbit["group_words"]


def apply_group(vector, word):
    return [v * (1 - 2 * ((word >> j) & 1)) for j, v in enumerate(vector)]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


middle_report = report["dense_sparse"]["middle"]
witness = middle_report["witness"]
record = patterns["middle_supports"][witness["support_index"]]
support = record["coordinates"]
middle_resolution = None
for word in group:
    projection = sum(((word >> coordinate) & 1) << j for j, coordinate in enumerate(support))
    retained_mask = witness["reachable_mask"] ^ projection
    if retained_mask in record["retained_masks"]:
        sparse = [0] * 39
        for j, coordinate in enumerate(support):
            sparse[coordinate] = 2 * (1 - 2 * ((retained_mask >> j) & 1))
        actual_dot = dot(apply_group(rows[witness["representative"]], word), sparse)
        check(actual_dot == 2 * middle_report["maximum_signed_sum"])
        middle_resolution = {"representative": witness["representative"], "group_word": word,
                             "support_index": witness["support_index"], "retained_mask": retained_mask,
                             "sparse_vector": sparse, "raw_integer_dot": actual_dot}
        break
check(middle_resolution is not None)

root_report = report["dense_sparse"]["root"]
witness = root_report["witness"]
row = rows[witness["representative"]]
root = [0] * 39
for coordinate in witness["pair"]:
    root[coordinate] = 4 if row[coordinate] >= 0 else -4
check(root in patterns["root_vectors"])
check(dot(row, root) == 4 * root_report["maximum_absolute_coordinate_pair_sum"])
root_resolution = {"representative": witness["representative"], "group_word": 0,
                   "sparse_vector": root, "raw_integer_dot": dot(row, root)}

star_report = report["dense_sparse"]["star"]
witness = star_report["witness"]
star_resolution = None
for star_index, star in enumerate(patterns["star_vectors"]):
    for word in group:
        if [v // 2 for v in apply_group(star, word)] == witness["transformed_star_half"]:
            actual_dot = dot(apply_group(rows[witness["representative"]], word), star)
            check(actual_dot == 2 * star_report["maximum_half_star_dot"])
            star_resolution = {"representative": witness["representative"], "group_word": word,
                               "star_index": star_index, "sparse_vector": star,
                               "raw_integer_dot": actual_dot}
            break
    if star_resolution is not None:
        break
check(star_resolution is not None)

result = {"status": "PASS", "arithmetic": "Pure Python arbitrary-precision integers",
          "scope": "Independent reconstruction of the reported maximum dense/middle, dense/root, and dense/star products as actual original sparse points against actual group-transformed representatives.",
          "middle": middle_resolution, "root": root_resolution, "star": star_resolution}
(directory / "actual_extreme_pair_witnesses39.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
