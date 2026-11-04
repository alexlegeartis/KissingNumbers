"""Bounded exact two-copy exchanges beyond the fixed-normal-octad sextets.

This is a research experiment, not a registered construction operator. Every
source point is retained initially; candidates are exact rational isometric
images. A complete source geometry replay plus direct-integer checks of every
inserted point independently certifies any accepted delta.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from collections import defaultdict
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys
import tempfile
import time

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

ROOT = Path(__file__).resolve().parents[3]
PACKAGE = ROOT / "verifications/improved/dim27-31-joint-refills"
NORMAL = {0, 2, 4, 5, 6, 10, 11, 23}
SOURCE = PACKAGE / "data/coordinates28_221988.npz"
SOURCE_SHA256 = "509b7894dd8534946a81d8bfd354c9ef947e4fe5c82482796aac76799169f129"


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def all_sextets():
    generator = np.array([1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1], dtype=np.uint8)
    G = np.zeros((12, 23), np.uint8)
    for i in range(12):
        G[i, i:i + 12] = generator
    labels = ((np.arange(4096)[:, None] >> np.arange(12)) & 1).astype(np.uint8)
    code = np.zeros((4096, 24), np.uint8)
    code[:, :23] = labels @ G % 2
    code[:, 23] = code[:, :23].sum(axis=1) % 2
    octads = [tuple(map(int, np.flatnonzero(word))) for word in code[code.sum(axis=1) == 8]]
    require(len(octads) == 759, "Golay octad count")
    containing = defaultdict(list)
    for octad in octads:
        mask = sum(1 << index for index in octad)
        for tetrad in combinations(octad, 4):
            containing[sum(1 << index for index in tetrad)].append(mask)
    frames = set()
    for tetrad, octad_masks in containing.items():
        blocks = [tuple(i for i in range(24) if tetrad >> i & 1)]
        blocks += [tuple(i for i in range(24) if (octad ^ tetrad) >> i & 1) for octad in octad_masks]
        require(len(blocks) == 6 and sorted(sum((list(row) for row in blocks), [])) == list(range(24)), "Golay sextet partition")
        frames.add(tuple(sorted(blocks)))
    require(len(frames) == 1771, "Full Golay sextet count")
    ordered = sorted(frames)
    old = [frame for frame in ordered if any(set(block) <= NORMAL for block in frame)]
    outside = [frame for frame in ordered if frame not in set(old)]
    require(len(old) == 35 and len(outside) == 1736, "Fixed-octad subfamily count")
    return ordered, outside


def transformation(blocks, dimension):
    matrix = 2 * np.eye(dimension, dtype=np.int64)
    for block in blocks:
        for i in block:
            for j in block:
                matrix[i, j] = 1 - 2 * int(i == j)
    matrix[list(blocks[0])] *= -1
    require(np.array_equal(matrix @ matrix.T, 4 * np.eye(dimension, dtype=np.int64)), "Exact sextet isometry")
    return matrix


def image(points, frames, word):
    dimension = points.shape[1]
    matrix = np.eye(dimension, dtype=np.int64)
    denominator = 1
    for index in word:
        matrix = transformation(frames[index], dimension) @ matrix
        denominator *= 2
    require(np.array_equal(matrix @ matrix.T, denominator**2 * np.eye(dimension, dtype=np.int64)), "Composite isometry identity")
    numerator = points.astype(np.int64) @ matrix.T
    common = gcd(denominator, int(np.gcd.reduce(abs(numerator).ravel())))
    factor = denominator // common
    original = points.astype(np.int64) * factor
    transformed = numerator // common
    require(np.max(abs(original)) < 32768 and np.max(abs(transformed)) < 32768, "Coordinate storage bound")
    norm = 128 * factor**2
    require(np.all(np.sum(transformed**2, axis=1) == norm), "Image norms")
    return original.astype(np.int16), transformed.astype(np.int16), norm, factor, matrix, denominator


def core_flags(points, factor, lookup):
    flags = np.zeros(len(points), bool)
    possible = np.flatnonzero(~np.any(points[:, 24:], axis=1) & ~np.any(points[:, :24] % factor, axis=1))
    for i in possible:
        flags[i] = (points[i, :24] // factor).astype(np.int16).tobytes() in lookup
    return flags


def exchange(original, transformed, norm, *, blocker_limit, max_checks, max_edges, deadline, core_lookup=None, factor=1):
    started = time.monotonic()
    source = {row.tobytes(): i for i, row in enumerate(original)}
    lookup = {row.tobytes(): i for i, row in enumerate(transformed)}
    require(len(source) == len(original) and len(lookup) == len(transformed), "Source/image injection")
    common = source.keys() & lookup.keys()
    ai = np.array([i for key, i in source.items() if key not in common], np.int32)
    bi = np.array([i for key, i in lookup.items() if key not in common], np.int32)
    record = {"common": len(common), "exclusive_source": len(ai), "exclusive_candidates": len(bi),
              "blocker_limit": blocker_limit, "dot_products": 0}
    if not len(bi):
        return original, {**record, "state": "complete", "gain": 0, "removed": 0, "added": 0}, None
    a, b = original[ai], transformed[bi]
    bound = original.shape[1] * int(max(abs(a).max(initial=0), abs(b).max(initial=0)))**2
    require(bound < 2**24, "Exact binary32 integer-product bound")
    af, bf = a.astype(np.float32), b.astype(np.float32)
    active = np.arange(len(bi))
    blockers = [[] for _ in bi]
    rng = np.random.default_rng(281060)
    core_a = np.zeros(len(a), bool) if core_lookup is None else core_flags(a, factor, core_lookup)
    core_b = np.zeros(len(b), bool) if core_lookup is None else core_flags(b, factor, core_lookup)
    groups = [rng.permutation(np.flatnonzero(~core_a)), rng.permutation(np.flatnonzero(core_a))]
    batches = [(group[at:at + 64], is_core) for is_core, group in enumerate(groups) for at in range(0, len(group), 64)]
    checks, skipped, stored_edges = 0, 0, 0
    for batch, (indices, is_core) in enumerate(batches):
        positions = np.flatnonzero(~core_b[active]) if is_core else np.arange(len(active))
        evaluated = active[positions]
        required = len(evaluated) * len(indices)
        if checks + required > max_checks or time.monotonic() >= deadline:
            return original, {**record, "state": "incomplete", "reason": "cross_work_or_wall_cutoff", "dot_products": checks,
                              "unresolved_candidates": len(active), "source_blocks_scanned": batch,
                              "certified_core_pairs_skipped": skipped}, None
        bad = bf[evaluated] @ af[indices].T > norm // 2
        checks += required
        skipped += (len(active) - len(evaluated)) * len(indices)
        counts = np.fromiter((len(blockers[i]) for i in evaluated), dtype=np.int32, count=len(evaluated)) + bad.sum(axis=1)
        good = counts <= blocker_limit
        retained = evaluated[good]
        for discarded in evaluated[~good]:
            stored_edges -= len(blockers[int(discarded)])
            blockers[int(discarded)] = []
        rows, columns = np.nonzero(bad[good])
        if stored_edges + len(rows) > max_edges:
            return original, {**record, "state": "incomplete", "reason": "edge_cutoff", "dot_products": checks,
                              "certified_core_pairs_skipped": skipped}, None
        for row, column in zip(rows, columns):
            blockers[int(retained[row])].append(int(indices[column]))
        stored_edges += len(rows)
        keep = np.ones(len(active), bool)
        keep[positions[~good]] = False
        active = active[keep]
        if not len(active):
            break
    selected_blockers = [blockers[int(i)] for i in active]
    edges = sum(map(len, selected_blockers))
    if edges > max_edges:
        return original, {**record, "state": "incomplete", "reason": "edge_cutoff", "dot_products": checks}, None
    indptr = np.r_[0, np.cumsum([len(row) for row in selected_blockers])]
    graph = csr_matrix((np.ones(edges, np.uint8), np.array(sum(selected_blockers, []), np.int32), indptr), shape=(len(active), len(a)))
    matching = maximum_bipartite_matching(graph, perm_type="column")
    size = int(np.count_nonzero(matching >= 0))
    require(len(set(map(int, matching[matching >= 0]))) == size, "Matching endpoint injection")
    for i, j in enumerate(matching):
        require(j < 0 or int(j) in selected_blockers[i], "Matching edge membership")
    reverse = np.full(len(a), -1, np.int64)
    for i, j in enumerate(matching):
        if j >= 0:
            reverse[j] = i
    left = set(map(int, np.flatnonzero(matching < 0)))
    right, queue = set(), list(left)
    for i in queue:
        for j in selected_blockers[i]:
            if matching[i] == j or j in right:
                continue
            right.add(j)
            partner = int(reverse[j])
            if partner >= 0 and partner not in left:
                left.add(partner)
                queue.append(partner)
    # Independent combinatorial optimality witness: the complement of the
    # alternating left set plus the reached right set is a minimum vertex cover.
    cover_left = set(range(len(active))) - left
    require(len(cover_left) + len(right) == size, "Matching/vertex-cover cardinality equality")
    require(all(i in cover_left or j in right for i, row in enumerate(selected_blockers) for j in row), "Vertex cover edge coverage")
    gain = len(active) - size
    keep = np.ones(len(original), bool)
    removed = ai[sorted(right)]
    added = bi[active[sorted(left)]]
    keep[removed] = False
    result = np.concatenate((original[keep], transformed[added]))
    require(len(result) == len(original) + gain, "Exact exchange count")
    record.update(state="complete", restricted_candidates=len(active), pruned_candidates=len(bi) - len(active),
                  dot_products=checks, edges=edges, matching_size=size, vertex_cover_size=size,
                  certified_core_pairs_skipped=skipped,
                  gain=gain, removed=len(removed), added=len(added), total=len(result),
                  arithmetic="exact bounded integer binary32 products; no tolerance", seconds=time.monotonic() - started,
                  scope="maximum compatible subset of this blocker-restricted two-copy union only")
    witness = {"removed_indices": removed, "added_image_indices": added,
               "restricted_image_indices": bi[active], "matching": matching,
               "cover_left": np.array(sorted(cover_left), np.int32), "cover_right": np.array(sorted(right), np.int32),
               "restricted_image_points": transformed[bi[active]], "blocker_indptr": indptr,
               "blocker_baseline_indices": ai[np.array(sum(selected_blockers, []), np.int32)]}
    return result, record, witness


def shared_blocker_search(baseline, pools, max_tests=1000000):
    """Point-level multi-image search: two additions/one deletion, then 3/2."""
    common_factor = 1
    from math import lcm
    for pool in pools:
        common_factor = lcm(common_factor, pool["factor"])
    norm = 128 * common_factor**2
    points, blockers, origins, seen = [], [], [], {}
    for pool_index, pool in enumerate(pools):
        scale = common_factor // pool["factor"]
        for i, row in enumerate(pool["points"]):
            start, stop = pool["indptr"][i:i + 2]
            blocked = tuple(sorted(map(int, pool["blockers"][start:stop])))
            if not 1 <= len(blocked) <= 2:
                continue
            point = (row.astype(np.int64) * scale).astype(np.int16)
            key = point.tobytes()
            if key in seen:
                require(blockers[seen[key]] == blocked, "Repeated image point has inconsistent exact blockers")
                origins[seen[key]].add(pool_index)
                continue
            seen[key] = len(points)
            points.append(point)
            blockers.append(blocked)
            origins.append({pool_index})
    buckets = defaultdict(list)
    for index, blocked in enumerate(blockers):
        buckets[blocked].append(index)
    tests, cache = 0, {}

    def compatible(i, j):
        nonlocal tests
        key = (min(i, j), max(i, j))
        if key not in cache:
            tests += 1
            if tests > max_tests:
                raise TimeoutError("Multi-image candidate-pair budget")
            cache[key] = 2 * int(points[i].astype(np.int64) @ points[j].astype(np.int64)) <= norm
        return cache[key]

    def found(indices, removed):
        original = (baseline.astype(np.int64) * common_factor).astype(np.int16)
        added = np.array([points[i] for i in indices], np.int16)
        keep = np.ones(len(baseline), bool)
        keep[list(removed)] = False
        final = np.concatenate((original[keep], added))
        proof = direct_delta(original, final, added)
        return {"state": "gain", "gain": len(indices) - len(removed), "removed_indices": list(removed),
                "added": added, "points": final, "factor": common_factor, "proof": proof,
                "pair_tests": tests, "unique_low_blocker_candidates": len(points),
                "source_pools": [sorted(origins[i]) for i in indices]}

    try:
        for blocked, indices in sorted(buckets.items()):
            if len(blocked) != 1:
                continue
            for i, j in combinations(indices, 2):
                if compatible(i, j):
                    return found([i, j], blocked)
        for blocked, middle in sorted(buckets.items()):
            if len(blocked) != 2:
                continue
            indices = sorted(set(middle + buckets.get((blocked[0],), []) + buckets.get((blocked[1],), [])))
            for i, j in combinations(indices, 2):
                if not compatible(i, j):
                    continue
                for k in indices:
                    if k > j and compatible(i, k) and compatible(j, k):
                        return found([i, j, k], blocked)
    except TimeoutError:
        return {"state": "incomplete", "reason": "candidate_pair_budget", "pair_tests": tests,
                "unique_low_blocker_candidates": len(points)}
    return {"state": "complete_no_gain_in_two_three_move_pool", "pair_tests": tests,
            "unique_low_blocker_candidates": len(points), "buckets": len(buckets),
            "scope": "two compatible additions sharing one deletion; three compatible additions with two deletions"}


def direct_delta(original, final, added):
    """Independent int64 products: all newly introduced pairs, no BLAS."""
    require(len(np.unique(final, axis=0)) == len(final), "Distinct output points")
    norm = int(np.sum(original[0].astype(np.int64)**2))
    require(np.all(np.sum(final.astype(np.int64)**2, axis=1) == norm), "Output norms")
    lookup = {row.tobytes(): i for i, row in enumerate(final)}
    maximum, count = -norm, 0
    values = final.astype(np.int64)
    for source in added:
        row = source.astype(np.int64)
        dots = np.zeros(len(final), np.int64)
        for j in np.flatnonzero(row):
            dots += values[:, j] * int(row[j])
        dots[lookup[source.tobytes()]] = -norm
        require(int(dots.max(initial=-norm)) * 2 <= norm, "Direct integer inserted-point incompatibility")
        maximum = max(maximum, int(dots.max(initial=-norm)))
        count += len(final) - 1
    old = {row.tobytes() for row in original}
    new = {row.tobytes() for row in added}
    require(all(row.tobytes() in old or row.tobytes() in new for row in final), "Output is retained source plus checked insertions")
    return {"status": "PASS", "direct_integer_pairs": count, "maximum_inserted_product": maximum,
            "norm_squared": norm, "cardinality": len(final), "ambient_dimension": final.shape[1],
            "proof_scope": "full freshly checked source geometry plus every inserted/final integer pair; unchanged pairs inherited"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--frames", type=int, default=8)
    parser.add_argument("--blocker-limit", type=int, default=8)
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--max-checks", type=int, default=1000000000)
    parser.add_argument("--max-edges", type=int, default=1000000)
    parser.add_argument("--composites", action="store_true")
    parser.add_argument("--outside-indices", help="explicit comma-separated outside-family indices")
    parser.add_argument("--core-skip", action="store_true", help="freshly verify full generated shell, then omit its internal cross pairs")
    parser.add_argument("--multi-image", action="store_true")
    args = parser.parse_args()
    require(args.frames > 0 and args.blocker_limit >= 0 and args.seconds > 0 and args.max_checks > 0 and args.max_edges > 0, "Search limits")
    args.out.mkdir(parents=True, exist_ok=False)
    require(hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256, "Pinned baseline file hash")
    spec = importlib.util.spec_from_file_location("_competitive28_source_verifier", PACKAGE / "verify_geometry.py")
    source_verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source_verifier)
    baseline_report = source_verifier.verify(SOURCE)
    require(baseline_report["points"] == 221988, "Baseline count")
    with np.load(SOURCE, allow_pickle=False) as archive:
        baseline = archive["points"].astype(np.int16)
    frames, outside = all_sextets()
    chosen = list(map(int, args.outside_indices.split(','))) if args.outside_indices else list(range(args.frames))
    require(all(0 <= i < len(outside) for i in chosen) and len(set(chosen)) == len(chosen), "Outside-frame indices")
    frame_indices = [frames.index(outside[i]) for i in chosen]
    words = [[index] for index in frame_indices]
    if args.composites and len(frame_indices) > 1:
        words += [[frame_indices[i], frame_indices[(i + 1) % len(frame_indices)]] for i in range(len(frame_indices))]
    # Do not charge identical matrix words as new search breadth.
    unique_words, matrix_keys = [], set()
    for word in words:
        matrix = np.eye(28, dtype=np.int64)
        denominator = 1
        for i in word:
            matrix = transformation(frames[i], 28) @ matrix
            denominator *= 2
        common = gcd(denominator, int(np.gcd.reduce(abs(matrix).ravel())))
        key = (denominator // common, (matrix // common).tobytes())
        if key not in matrix_keys:
            matrix_keys.add(key)
            unique_words.append(word)
    words = unique_words
    core_lookup, core_report = None, None
    if args.core_skip:
        sys.path.insert(0, str(ROOT / "common"))
        import leech
        shell = np.pad(2 * leech.build().astype(np.int16), ((0, 0), (0, 4)))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'shell.npz'
            np.savez_compressed(path, points=shell, squared_norm=128, inner_product_limit=64, coordinate_denominator=2)
            core_report = source_verifier.verify(path)
        core_lookup = {row[:24].tobytes() for row in shell}
        require(len(core_lookup) == 196560, "Regenerated shell membership set")
    started, records, best, pools = time.monotonic(), [], len(baseline), []
    deadline = started + args.seconds
    for ordinal, word in enumerate(words):
        if time.monotonic() >= deadline:
            break
        original, transformed, norm, factor, matrix, denominator = image(baseline, frames, word)
        result, record, witness = exchange(original, transformed, norm, blocker_limit=args.blocker_limit,
                                           max_checks=args.max_checks, max_edges=args.max_edges, deadline=deadline,
                                           core_lookup=core_lookup, factor=factor)
        record.update(word=word, word_frames=[frames[i] for i in word], common_integer_scale=factor,
                      original_relative_fractional_image_points=int(np.any(transformed % factor != 0, axis=1).sum()) if factor > 1 else 0,
                      norm_squared=norm)
        if record["state"] == "complete" and witness is not None:
            pool = {"points": witness["restricted_image_points"], "indptr": witness["blocker_indptr"],
                    "blockers": witness["blocker_baseline_indices"], "factor": factor}
            pools.append(pool)
            with (args.out / f"pool-{ordinal:03d}.npz").open("xb") as stream:
                np.savez_compressed(stream, **pool)
        if record.get("gain", 0) > 0:
            added = transformed[witness["added_image_indices"]]
            record["independent_delta"] = direct_delta(original, result, added)
            path = args.out / f"candidate-{ordinal:03d}.npz"
            with path.open("xb") as stream:
                np.savez_compressed(stream, points=result, squared_norm=norm, inner_product_limit=norm // 2,
                                    coordinate_denominator=2 * factor)
            with (args.out / f"witness-{ordinal:03d}.npz").open("xb") as stream:
                np.savez_compressed(stream, **witness, transform_numerator=matrix, transform_denominator=denominator)
            record["candidate"] = path.name
            record["candidate_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            best = max(best, len(result))
        records.append(record)
        print(json.dumps(record), flush=True)
    multi = shared_blocker_search(baseline, pools) if args.multi_image else {"state": "not_requested"}
    if multi.get("state") == "gain":
        with (args.out / "multi-image-gain.npz").open("xb") as stream:
            np.savez_compressed(stream, points=multi.pop("points"), added=multi.pop("added"),
                                removed_indices=multi["removed_indices"], squared_norm=128 * multi["factor"]**2,
                                inner_product_limit=64 * multi["factor"]**2, coordinate_denominator=2 * multi["factor"])
        best = max(best, len(baseline) + multi["gain"])
    report = {"source": str(SOURCE.relative_to(ROOT)), "source_sha256": SOURCE_SHA256,
              "source_fresh_full_verification": baseline_report,
              "fresh_shared_core_verification": core_report, "multi_image": multi,
              "full_sextet_partitions": len(frames), "outside_normal_octad_subfamily": len(outside),
              "planned_words": words, "attempted_words": len(records), "records": records,
              "best_cardinality": best, "seconds": time.monotonic() - started,
              "limits": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
              "scope": "bounded restricted two-copy exchanges from the current competitive baseline; no unrestricted nonexistence or global priority claim",
              "status": "complete_in_declared_restricted_pool" if len(records) == len(words) and all(row["state"] == "complete" for row in records) and multi.get("state") != "incomplete" else "incomplete"}
    with (args.out / "report.json").open("x") as stream:
        stream.write(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"report": str(args.out / "report.json"), "best": best, "status": report["status"]}), flush=True)


if __name__ == "__main__":
    main()
