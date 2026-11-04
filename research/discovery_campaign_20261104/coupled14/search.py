"""Seeded exact two-copy exchanges outside the fixed E7 x E7 alphabet.

Run from repository root, e.g.:
OPENBLAS_NUM_THREADS=1 python research/discovery_campaign_20261104/coupled14/search.py --trials 360 --out research/discovery_campaign_20261104/coupled14/run360
No floating screen is used. Every tested graph has all exact cross pairs and
a checked equal-size matching/cover certificate. Independent full coordinate
verification is run on the baseline and every new best, including ties first.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import random
import sys
from time import monotonic

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from discovery.core import Artifact, RunContext, digest
from discovery.certification import verify_artifact
from discovery.engines.groups import _build14
from discovery.engines.selection import select_isometric_rows


def normal(support, coefficients):
    result = [0] * 14
    for index, value in zip(support, coefficients):
        result[index] = value
    return result


def candidate_parameters(seed, trials):
    """Finite seeded list; no claim to enumerate rational O(14)."""
    rng = random.Random(seed)
    candidates = []
    # Local new-alphabet reflectors and mixed-coordinate reflectors.
    patterns = [[1, 2], [1, 1, 1], [1, 1, 2], [1, 2, 2], [1, 1, 1, 1],
                [1, 1, 1, 2], [1, 1, 2, 2], [1, 1, 1, 1, 1], [1]*7]
    for pattern in patterns:
        for placement in range(6):
            width = len(pattern)
            if placement < 2:
                support = rng.sample(range(7*placement, 7*(placement+1)), width)
            else:
                support = rng.sample(range(14), width)
            candidates.append(("structured_reflection", {"normals": [normal(support, pattern)]}))
    # Plane rotations with exact rational cos/sin, including near/block mixing.
    for p, q in [(1, 2), (1, 3), (2, 3), (1, 4), (3, 4), (2, 5)]:
        denominator = p*p + q*q
        cosine, sine = q*q-p*p, 2*p*q
        for placement in range(6):
            axes = rng.sample(range(14), 2) if placement >= 2 else rng.sample(range(7*placement, 7*(placement+1)), 2)
            matrix = [[int(i == j) for j in range(14)] for i in range(14)]
            a, b = axes
            matrix[a][a] = matrix[b][b] = [cosine, denominator]
            matrix[a][b], matrix[b][a] = [sine, denominator], [-sine, denominator]
            candidates.append(("rational_plane_rotation", {"matrix": matrix}))
    while len(candidates) < trials:
        index = len(candidates) - 90
        mode = index % 6
        length = 1 if mode in (0, 1, 2) else (2 if mode in (3, 4) else 4)
        normals = []
        for step in range(length if mode != 5 else 2):
            if mode in (0, 3):
                side = rng.randrange(2)
                support = rng.sample(range(7*side, 7*(side+1)), rng.randint(2, 7))
            else:
                support = rng.sample(range(14), rng.randint(2, 14))
            normals.append(normal(support, [rng.choice([-2, -1, 1, 2]) for _ in support]))
        if mode == 5:
            # Rational commutator; both blocks can participate in each normal.
            normals = [normals[0], normals[1], normals[0], normals[1]]
        family = ["random_block_reflection", "random_mixed_reflection", "random_dense_reflection",
                  "block_reflection_word", "mixed_reflection_word", "mixed_commutator"][mode]
        candidates.append((family, {"normals": normals}))
    return candidates[:trials]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20261104)
    parser.add_argument("--trials", type=int, default=360)
    parser.add_argument("--seconds", type=float, default=600)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.trials <= 2000 or not 0 < args.seconds <= 3600:
        raise ValueError("Bound trials to 1..2000 and seconds to (0,3600]")
    args.out.mkdir(parents=True, exist_ok=False)
    context = RunContext(args.out)
    rows, frames = _build14()
    baseline = Artifact("spherical_code", {"coordinates": rows, "weights": [1]*14, "norm_squared": 8})
    baseline_report = verify_artifact(baseline, context)
    assert baseline_report.status == "exact"
    assert baseline_report.details["facts"]["cardinality"] == 1932
    (args.out/"baseline_verification.json").write_text(json.dumps(baseline_report.to_dict(), indent=2)+"\n")
    np.savez_compressed(args.out/"baseline.npz", coordinates=np.asarray(rows), weights=np.ones(14, dtype=int), norm_squared=8)
    candidates = candidate_parameters(args.seed, args.trials)
    start = monotonic()
    best, best_index = 1931, None
    family_counts, size_counts = Counter(), Counter()
    edge_min, edge_max, maximum_new_points = None, None, 0
    attempted, stop = 0, "finite_list_exhausted"
    with (args.out/"trials.jsonl").open("w") as summaries, gzip.open(args.out/"certificates.jsonl.gz", "wt") as certificates:
        for index, (family, params) in enumerate(candidates):
            if monotonic()-start >= args.seconds:
                stop = "time_limit_reached"
                break
            tick = monotonic()
            data, certificate = select_isometric_rows(rows, [1]*14, 8, params, context)
            count = len(data["coordinates"])
            full_verified = False
            if count > best:
                result = Artifact("spherical_code", data)
                report = verify_artifact(result, context)
                assert report.status == "exact" and report.details["facts"]["cardinality"] == count
                full_verified = True
                best, best_index = count, index
                (args.out/"best_artifact.json").write_text(json.dumps({"kind": "spherical_code", "data": data,
                    "facts": report.details["facts"]}, separators=(",", ":"))+"\n")
                (args.out/"best_verification.json").write_text(json.dumps(report.to_dict(), indent=2)+"\n")
                (args.out/"best_exchange.json").write_text(json.dumps({"params": params, "certificate": certificate}, separators=(",", ":"))+"\n")
            summary = {"index": index, "family": family, "params": params, "cardinality": count,
                       "gain_over_1932": count-1932, "overlap_points": certificate["overlap_points"],
                       "image_new_points": certificate["image_new_points"], "conflict_edges": certificate["conflict_edges"],
                       "matching_size": len(certificate["matching"]), "cross_pairs_checked": certificate["cross_pairs_checked"],
                       "certificate_sha256": digest(certificate), "full_output_independently_verified": full_verified,
                       "seconds": monotonic()-tick}
            summaries.write(json.dumps(summary, separators=(",", ":"))+"\n")
            summaries.flush()
            certificates.write(json.dumps({"index": index, "certificate": certificate}, separators=(",", ":"))+"\n")
            family_counts[family] += 1
            size_counts[count] += 1
            edge_min = certificate["conflict_edges"] if edge_min is None else min(edge_min, certificate["conflict_edges"])
            edge_max = certificate["conflict_edges"] if edge_max is None else max(edge_max, certificate["conflict_edges"])
            maximum_new_points = max(maximum_new_points, certificate["image_new_points"])
            attempted += 1
            if attempted % 10 == 0 or count > 1932:
                print(json.dumps({"attempted": attempted, "family": family, "count": count, "best": best,
                                  "elapsed_seconds": monotonic()-start}), flush=True)
    source_paths = [Path(__file__), Path("discovery/engines/selection.py"), Path("discovery/certification.py"), Path("discovery/engines/groups.py")]
    report = {"seed": args.seed, "planned_trials": args.trials, "attempted_trials": attempted,
              "finite_list_exhausted": attempted == args.trials, "stop_reason": stop,
              "best_cardinality": best, "best_trial_index": best_index, "baseline_cardinality": 1932,
              "new_sota_found": best > 1932, "family_counts": dict(family_counts),
              "cardinality_histogram": dict(size_counts), "conflict_edges_range": [edge_min, edge_max],
              "maximum_new_image_points": maximum_new_points, "elapsed_seconds": monotonic()-start,
              "scope": "These finite two-isometric-copy unions only. No global bound or exhaustion of rational maps.",
              "baseline_proof": "full independent integer norm/distinctness/pair/rank verification",
              "per_trial_proof": "exact isometry, full cross graph, matching and equal-size vertex cover checked against every edge",
              "output_proof": "full independent integer verification of the first best and every strict improvement"}
    report["source_files_sha256_at_report_time"] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths}
    (args.out/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
