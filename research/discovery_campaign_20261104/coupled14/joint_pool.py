"""Bounded general-graph trials with heterogeneous exact product masses.

All candidate pairs are checked with exact integer coordinates before a
numerical binary MILP. Only independently replayed final coordinates are
certificates; solver bounds, timeouts and infeasibility are not exact proofs.
"""
import argparse
from collections import Counter
import json
from math import gcd, isfinite, lcm
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, csr_matrix, vstack

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from discovery.core import Artifact, RunContext, digest
from discovery.certification import verify_artifact
from discovery.engines.groups import _build14
from discovery.engines.selection import rational_isometry, transform_rows


def build_pool(mode, context):
    baseline, _ = _build14()
    if mode == "isometries":
        first, second = [0]*14, [0]*14
        first[0], first[7] = 1, 2
        second[1], second[7] = 2, 1
        words = [[first], [second], [first, second]]
        maps = [rational_isometry(14, [1]*14, {"normals": word}, context) for word in words]
        denominator = lcm(*(d for _, d in maps))
        rows = [[denominator*x for x in row] for row in baseline]
        groups = [0]*len(rows)
        seen = set(map(tuple, rows))
        for group, (numerator, d) in enumerate(maps, 1):
            _, images = transform_rows(baseline, numerator, d, context)
            for row in images:
                image = [x*(denominator//d) for x in row]
                if tuple(image) not in seen:
                    seen.add(tuple(image))
                    rows.append(image)
                    groups.append(group)
        return rows, groups, 8*denominator**2, [{"normals": word} for word in words]
    weights = [(7, 17, 13), (17, 7, 13), (359, 439, 401), (439, 359, 401)]
    denominator = lcm(*(d for _, _, d in weights))
    rows = [[denominator*x for x in row] for row in baseline]
    groups = [0] * len(rows)
    seen = set(map(tuple, rows))
    for group, (a, b, d) in enumerate(weights, 1):
        if a*a + b*b != 2*d*d:
            raise ValueError("Product masses do not preserve the exact norm")
        for row in baseline[:1764]:
            image = [a*(denominator//d)*x for x in row[:7]] + [b*(denominator//d)*x for x in row[7:]]
            if tuple(image) not in seen:
                seen.add(tuple(image))
                rows.append(image)
                groups.append(group)
    return rows, groups, 8*denominator**2, weights


def full_graph(rows, norm, context):
    count, width = len(rows), len(rows[0])
    context.limit("max_materialized_points", count)
    context.limit("max_dimension", width)
    context.limit("max_pairs", count*(count-1)//2)
    bound = max(abs(x) for row in rows for x in row)**2 * width
    if bound >= 2**62:
        raise ValueError("This bounded pool driver requires a proved int64 dot bound")
    coordinates = np.asarray(rows, dtype=np.int64)
    if np.any(np.sum(coordinates*coordinates, axis=1) != norm):
        raise ValueError("Candidate pool has unequal exact integer norms")
    left_chunks, right_chunks, edges = [], [], 0
    for start in range(0, count, 64):
        dots = coordinates[start:start+64] @ coordinates.T
        for offset, row in enumerate(dots):
            i = start+offset
            neighbors = np.flatnonzero(row[i+1:] > norm//2) + i+1
            if len(neighbors):
                edges += len(neighbors)
                context.limit("max_gram_entries", edges)
                left_chunks.append(np.full(len(neighbors), i, dtype=np.int32))
                right_chunks.append(neighbors.astype(np.int32))
    left = np.concatenate(left_chunks) if left_chunks else np.array([], dtype=np.int32)
    right = np.concatenate(right_chunks) if right_chunks else np.array([], dtype=np.int32)
    return left, right


def simplify(rows, norm):
    common = 0
    for row in rows:
        for value in row:
            common = gcd(common, value)
    if common > 1:
        rows = [[x//common for x in row] for row in rows]
        norm //= common**2
    return {"coordinates": rows, "weights": [1]*14, "norm_squared": norm}


def optimize_pool(rows, groups, norm, left, right, seconds, deleted_limit, force_mixture):
    baseline_count, count = 1932, len(rows)
    if np.any(right < baseline_count):
        raise ValueError("The freshly verified baseline has an internal conflict")
    edge_count = len(left)
    cross = left < baseline_count
    # y=deleted baseline point, x=selected candidate. Every conflict is retained:
    # cross edges impose x <= y; candidate edges impose x_u+x_v <= 1.
    row_ids = np.repeat(np.arange(edge_count, dtype=np.int32), 2)
    columns = np.column_stack((left, right)).ravel()
    values = np.ones(2*edge_count, dtype=np.float64)
    values[2*np.flatnonzero(cross)] = -1
    matrix = coo_matrix((values, (row_ids, columns)), shape=(edge_count, count)).tocsr()
    lower, upper = np.full(edge_count, -np.inf), np.where(cross, 0., 1.)
    deletion_row = csr_matrix((np.ones(baseline_count), (np.zeros(baseline_count, dtype=int), np.arange(baseline_count))), shape=(1, count))
    extra_rows, extra_lower, extra_upper = [deletion_row], [-np.inf], [deleted_limit]
    if force_mixture:
        for group in sorted(set(groups)-{0}):
            indices = np.flatnonzero(np.asarray(groups) == group)
            extra_rows.append(csr_matrix((np.ones(len(indices)), (np.zeros(len(indices), dtype=int), indices)), shape=(1, count)))
            extra_lower.append(1)
            extra_upper.append(np.inf)
    matrix = vstack([matrix, *extra_rows], format="csc")
    lower = np.r_[lower, extra_lower]
    upper = np.r_[upper, extra_upper]
    objective = np.r_[np.ones(baseline_count), -np.ones(count-baseline_count)]
    tick = monotonic()
    result = milp(objective, integrality=np.ones(count), bounds=Bounds(0, 1),
                  constraints=LinearConstraint(matrix, lower, upper),
                  options={"time_limit": seconds, "mip_rel_gap": 0., "presolve": True})
    def finite_diagnostic(value):
        return float(value) if value is not None and isfinite(float(value)) else None
    report = {"force_all_candidate_cohorts": force_mixture, "deleted_baseline_limit": deleted_limit,
              "solver_seconds_limit": seconds, "seconds": monotonic()-tick,
              "solver_status": int(result.status), "solver_message": str(result.message),
              "numerical_mip_gap": finite_diagnostic(getattr(result, "mip_gap", None)),
              "numerical_gain_upper_bound": finite_diagnostic(-result.mip_dual_bound) if getattr(result, "mip_dual_bound", None) is not None else None,
              "proof_scope": "Numerical search diagnostics only. Final coordinates require independent exact verification."}
    if result.x is None:
        report["incumbent_found"] = False
        return None, report
    integer_solution = np.rint(result.x).astype(np.int64)
    if np.max(np.abs(result.x-integer_solution)) > 1e-5:
        raise ValueError("MILP incumbent is not binary to tolerance")
    deleted = integer_solution[:baseline_count].astype(bool)
    chosen = np.r_[~deleted, integer_solution[baseline_count:].astype(bool)]
    if np.any(chosen[left] & chosen[right]):
        raise ValueError("Rounded MILP incumbent has an exact conflict edge")
    if deleted.sum() > deleted_limit:
        raise ValueError("Deletion budget violated")
    group_counts = Counter(groups[i] for i in np.flatnonzero(chosen))
    if force_mixture and any(group_counts[group] == 0 for group in set(groups)-{0}):
        raise ValueError("Mixture cohort constraint violated")
    selected_rows = [rows[i] for i in np.flatnonzero(chosen)]
    report.update(incumbent_found=True, cardinality=len(selected_rows), gain_over_1932=len(selected_rows)-1932,
                  retained_group_counts=dict(group_counts), deleted_baseline_points=int(deleted.sum()),
                  deleted_accessory_points=int(deleted[1764:].sum()),
                  chosen_indices=np.flatnonzero(chosen).tolist())
    return simplify(selected_rows, norm), report


def greedy_mixture(rows, groups, norm, left, right):
    """Concrete four-cohort fallback; all edges still checked independently."""
    n, count = 1932, len(rows)
    blockers = [set() for _ in rows]
    candidate_neighbors = [set() for _ in rows]
    for u, v in zip(left, right):
        if u < n:
            blockers[v].add(int(u))
        else:
            candidate_neighbors[u].add(int(v))
            candidate_neighbors[v].add(int(u))
    deleted, selected, forbidden = set(), set(), set()
    for group in sorted(set(groups)-{0}):
        available = [i for i in range(n, count) if groups[i] == group and i not in forbidden]
        chosen = min(available, key=lambda i: (len(blockers[i]-deleted), i))
        selected.add(chosen)
        deleted.update(blockers[chosen])
        forbidden.update(candidate_neighbors[chosen])
    # Refill every candidate whose blockers were already paid for. This is a
    # point-level mixture; no uniform mass layer is forced to be retained.
    for i in sorted(range(n, count), key=lambda i: (len(blockers[i]), i)):
        if i not in selected and i not in forbidden and blockers[i] <= deleted:
            selected.add(i)
            forbidden.update(candidate_neighbors[i])
    indices = sorted(set(range(n))-deleted | selected)
    chosen = np.zeros(count, dtype=bool)
    chosen[indices] = True
    if np.any(chosen[left] & chosen[right]):
        raise ValueError("Greedy mixture contains a graph conflict")
    return simplify([rows[i] for i in indices], norm), {"method": "shared deletion greedy mixture",
        "cardinality": len(indices), "gain_over_1932": len(indices)-1932,
        "retained_group_counts": dict(Counter(groups[i] for i in indices)),
        "deleted_baseline_points": len(deleted), "deleted_accessory_points": sum(i >= 1764 for i in deleted),
        "chosen_indices": indices}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--mode", choices=("weights", "isometries"), default="weights")
    parser.add_argument("--seconds", type=float, default=30)
    parser.add_argument("--deleted-limit", type=int, default=128)
    args = parser.parse_args()
    if not 0 < args.seconds <= 120 or not 1 <= args.deleted_limit <= 1932:
        raise ValueError("This experiment is limited to 120 seconds per solve and 1932 deletions")
    args.out.mkdir(parents=True, exist_ok=False)
    context = RunContext(args.out, {"max_materialized_points": 10000, "max_pairs": 50_000_000, "max_gram_entries": 3_000_000})
    start = monotonic()
    baseline, _ = _build14()
    base_report = verify_artifact(Artifact("spherical_code", {"coordinates": baseline, "norm_squared": 8}), context)
    rows, groups, norm, weights = build_pool(args.mode, context)
    left, right = full_graph(rows, norm, context)
    count = len(rows)
    cross = left < 1932
    degrees = np.bincount(right[cross]-1932, minlength=count-1932)
    inventory = {"mode": args.mode, "pool_cardinality": count, "baseline_cardinality": 1932, "cohort_descriptions": weights,
                 "cohort_counts": dict(Counter(groups)), "pairs_checked": count*(count-1)//2,
                 "conflict_edges": len(left), "baseline_candidate_edges": int(cross.sum()),
                 "candidate_candidate_edges": int((~cross).sum()), "blocker_degree_histogram": dict(Counter(map(int, degrees))),
                 "norm_squared": norm, "pool_data_sha256": digest(rows),
                 "baseline_verification": base_report.to_dict(), "all_accessory_points_jointly_removable": True}
    (args.out/"inventory.json").write_text(json.dumps(inventory, indent=2)+"\n")
    np.savez_compressed(args.out/"pool_and_graph.npz", coordinates=np.asarray(rows, dtype=np.int64),
                        norm_squared=norm, groups=np.asarray(groups), edge_left=left, edge_right=right)
    print(json.dumps({key: value for key, value in inventory.items() if key not in ("baseline_verification", "blocker_degree_histogram")}), flush=True)
    results = []
    fallback_data, fallback_report = greedy_mixture(rows, groups, norm, left, right)
    fallback_proof = verify_artifact(Artifact("spherical_code", fallback_data), context)
    fallback_report["independent_verification"] = fallback_proof.to_dict()
    (args.out/"greedy_mixture_artifact.json").write_text(json.dumps({"kind": "spherical_code", "data": fallback_data}, separators=(",", ":"))+"\n")
    results.append(fallback_report)
    for forced in (False, True):
        data, report = optimize_pool(rows, groups, norm, left, right, args.seconds, args.deleted_limit, forced)
        if data is not None:
            proof = verify_artifact(Artifact("spherical_code", data), context)
            report["independent_verification"] = proof.to_dict()
            name = "forced_mixture" if forced else "unrestricted"
            (args.out/(name+"_artifact.json")).write_text(json.dumps({"kind": "spherical_code", "data": data}, separators=(",", ":"))+"\n")
        results.append(report)
        (args.out/"results.json").write_text(json.dumps({"results": results, "elapsed_seconds": monotonic()-start,
            "new_sota_found": any(row.get("cardinality", 0) > 1932 for row in results),
            "scope": "Bounded heterogeneous-mass pool and deletion budget only. Numerical solver bounds do not prove exclusions."}, indent=2)+"\n")
        print(json.dumps({key: value for key, value in report.items() if key not in ("chosen_indices", "independent_verification")}), flush=True)


if __name__ == "__main__":
    main()
