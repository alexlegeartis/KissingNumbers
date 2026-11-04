"""Independent, solver-free replay of no-gain two-copy matching witnesses.

Each right-only point is paired with a different original point at normalized
inner product >1/2.  These disjoint conflicting pairs, plus the remaining
original singletons, partition the union into 1932 cliques. Thus its independent
set number is at most 1932. The independently verified baseline attains 1932.
Only matched conflict edges are needed for this upper bound, so replay does
not have to rebuild the optimizer's entire conflict graph or trust its cover.
"""
import argparse
import gzip
import json
from pathlib import Path
import sys
from time import monotonic

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from discovery.core import Artifact, RunContext, digest
from discovery.certification import verify_artifact
from discovery.engines.groups import _build14


def replay(directory):
    context = RunContext(directory)
    rows, _ = _build14()
    baseline = Artifact("spherical_code", {"coordinates": rows, "weights": [1]*14, "norm_squared": 8})
    baseline_report = verify_artifact(baseline, context)
    if baseline_report.status != "exact" or baseline_report.details["cardinality"] != 1932:
        raise ValueError("Baseline did not independently verify at 1932")
    start = monotonic()
    summaries = [json.loads(line) for line in (directory/"trials.jsonl").read_text().splitlines()]
    count = matched_edges = positive_gains = 0
    with gzip.open(directory/"certificates.jsonl.gz", "rt") as stream:
        for line in stream:
            raw = json.loads(line)
            index, certificate = raw["index"], raw["certificate"]
            if index != count or summaries[index]["certificate_sha256"] != digest(certificate):
                raise ValueError("Certificate ordering or content digest mismatch")
            numerator, denominator = certificate["matrix_numerator"], certificate["matrix_denominator"]
            if type(denominator) is not int or denominator <= 0 or len(numerator) != 14 or any(len(row) != 14 for row in numerator):
                raise ValueError("Malformed integer rational map")
            if any(type(value) is not int for row in numerator for value in row):
                raise ValueError("Map entries must be integers")
            # Independent integer isometry proof in the unweighted 14D metric.
            for i in range(14):
                for j in range(i, 14):
                    dot = sum(numerator[i][k]*numerator[j][k] for k in range(14))
                    if dot != (denominator**2 if i == j else 0):
                        raise ValueError("Map is not an exact rational isometry")
            bound = 2 * max(sum(abs(numerator[k][j]) for k in range(14)) for j in range(14))
            dtype = np.int64 if bound < 2**62 and 2*denominator < 2**62 else object
            image = np.asarray(rows, dtype=dtype) @ np.asarray(numerator, dtype=dtype)
            left = [[denominator*x for x in row] for row in rows]
            existing = set(map(tuple, left))
            right_indices = [i for i, row in enumerate(image) if tuple(row) not in existing]
            right = [image[i] for i in right_indices]
            if certificate["right_source_indices"] != right_indices:
                raise ValueError("Overlap deduplication differs from the witness")
            pairs = certificate["matching"]
            if len(pairs) != len(right):
                # Positive trials require full selected-witness replay rather
                # than a no-gain conclusion from an unsaturated matching.
                positive_gains += 1
                raise ValueError("Right side is not saturated; replay the positive coordinate witness separately")
            used_left, used_right = set(), set()
            for u, v in pairs:
                if type(u) is not int or type(v) is not int or not 0 <= u < len(left) or not 0 <= v < len(right):
                    raise ValueError("Invalid matching endpoint")
                if u in used_left or v in used_right:
                    raise ValueError("Matching endpoints are not distinct")
                used_left.add(u)
                used_right.add(v)
                dot = sum(a*int(b) for a, b in zip(left[u], right[v]))
                if 2*dot <= 8*denominator**2:
                    raise ValueError("Matching pair is not an actual conflicting geometric pair")
            if used_right != set(range(len(right))):
                raise ValueError("A transformed new point is unmatched")
            if summaries[index]["cardinality"] != 1932:
                raise ValueError("Reported optimum disagrees with the exact clique-partition bound")
            count += 1
            matched_edges += len(pairs)
    if count != len(summaries):
        raise ValueError("Witness count does not match summary count")
    best = Artifact.from_dict(json.loads((directory/"best_artifact.json").read_text()))
    best_report = verify_artifact(best, context)
    return {"status": "exact", "trials_replayed": count, "matched_conflicts_verified": matched_edges,
            "all_tested_union_optima": 1932, "baseline_pairs_checked": baseline_report.details["pairs_checked"],
            "best_pairs_checked": best_report.details["pairs_checked"], "best_cardinality": best_report.details["cardinality"],
            "elapsed_seconds": monotonic()-start,
            "method": "fresh full baseline and best geometry; exact integer rational isometries, exact overlaps and saturating disjoint conflict pairs",
            "scope": "Only the enumerated two-copy unions; no global upper bound or novelty claim"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    result = replay(args.directory)
    (args.directory/"replay_report.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
