"""Take exact tail hyperplane sections of the 223,528-point 28D witness."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations, permutations, product
import json
from multiprocessing import get_context

import numpy as np

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, maximum_bipartite_matching

from fusion import cross_conflicts
from source import ROOT, keys, leech_shell, load_points


def complete_normals(tails):
    # A triple determines the span of any rank-three section.
    normals = set()
    for triple in combinations(tails.astype(np.int64), 3):
        matrix = np.stack(triple)
        normal = np.empty(4, dtype=np.int64)
        for i in range(4):
            a, b, c = np.delete(matrix, i, axis=1)
            normal[i] = (-1) ** i * (np.cross(a, b) @ c)
        if np.any(normal):
            normal //= np.gcd.reduce(np.abs(normal))
            normals.add(tuple(normal * np.sign(normal[np.flatnonzero(normal)[0]])))
    return sorted(normals)


def root_blockers(tails, index, deleted, auxiliary):
    blocked = np.zeros((len(tails), len(deleted)), dtype=bool)
    root_matrix = deleted.astype(np.float32).T
    for t in range(len(tails)):
        bodies = auxiliary[index == t, :24].astype(np.float32)
        for start in range(0, len(bodies), 128):
            blocked[t] |= np.any(bodies[start:start + 128] @ root_matrix > 64, axis=0)
    return blocked


def aligned_sections():
    # An integral isometry is available for each coordinate or diagonal normal.
    for axis in range(4):
        yield np.delete(2 * np.eye(4, dtype=np.int16), axis, axis=1)
    frame = np.array([[1, 1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, 1]], dtype=np.int16)
    for signs in product((-1, 1), repeat=3):
        yield np.array([1, *signs], dtype=np.int16)[:, None] * frame


def check_alignment(job):
    section, order, signs = job
    layer = sections27[section].copy()
    layer[:, 24:] = layer[:, 24:][:, order] * signs
    candidates = np.r_[layer, missing_roots]
    graph = cross_conflicts(candidates, affected)
    match = maximum_bipartite_matching(graph, perm_type="column")
    matched = match >= 0
    assert np.all(np.sum(candidates[matched] * affected[match[matched]], axis=1) > 64)
    if matched.all():
        return 0, 0

    # Outside these alternating paths, each candidate has its own blocker
    # outside the neighborhood of the reached candidates.
    inverse = np.full(len(affected), -1)
    inverse[match[matched]] = np.flatnonzero(matched)
    reached = ~matched
    blockers = np.zeros(len(affected), dtype=bool)
    queue = np.flatnonzero(reached).tolist()
    while queue:
        i = queue.pop()
        neighbors = graph.indices[graph.indptr[i]:graph.indptr[i + 1]]
        new = neighbors[~blockers[neighbors]]
        blockers[new] = True
        for j in inverse[new]:
            assert j >= 0
            if not reached[j]:
                reached[j] = True
                queue.append(j)

    selected = candidates[reached]
    old = affected[blockers]
    edges = graph[reached][:, blockers].tocoo()
    n = len(selected)
    subgraph = csr_matrix((edges.data, (edges.row, n + edges.col)), shape=(n + len(old),) * 2)
    count, labels = connected_components(subgraph, directed=False)
    largest_gain = 0
    for component in range(count):
        new = selected[labels[:n] == component]
        base = old[labels[n:] == component]
        # The largest component on the named input has twelve candidates.
        if len(new) > 20:
            raise ValueError("Reached component exceeds the exhaustive subset limit")
        subsets = (np.arange(1 << len(new))[:, None] & (1 << np.arange(len(new)))) > 0
        conflicts = new @ new.T > 64
        np.fill_diagonal(conflicts, False)
        valid = subsets[~np.any((subsets @ conflicts) & subsets, axis=1)]
        cuts = valid @ (new @ base.T > 64)
        profits = valid.sum(axis=1) - cuts.sum(axis=1)
        best = int(np.argmax(profits))
        if profits[best] > largest_gain:
            largest_gain = int(profits[best])
            keep = ~np.isin(keys(source), keys(base[cuts[best]]))
            witness = np.r_[source[keep], new[valid[best]]]
            np.savez_compressed(f"section_win_{section}_{order}_{signs}.npz", points=witness,
                                squared_norm=128, inner_product_limit=64)
    return largest_gain, count


def main():
    global source, missing_roots, affected, sections27
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--alignments", type=int, default=576, choices=range(577), metavar="0..576")
    args = parser.parse_args()
    package = ROOT / "verifications/improved/dim28-takhanov-root-exchange/data"
    with np.load(package / "coordinates28_223528.npz") as archive:
        p = archive["points"].astype(np.int16)
    shell = 2 * leech_shell().astype(np.int16)
    auxiliary = p[np.any(p[:, 24:] != 0, axis=1)]
    old_roots = p[np.all(p[:, 24:] == 0, axis=1), :24]
    deleted = shell[~np.isin(keys(shell), keys(old_roots), assume_unique=True)]
    tails, index = np.unique(auxiliary[:, 24:], return_inverse=True, axis=0)
    blocked = root_blockers(tails, index, deleted, auxiliary)
    counts = Counter()
    for normal in complete_normals(tails):
        keep = tails @ normal == 0
        size = len(shell) + np.sum(keep[index]) - np.any(blocked[keep], axis=0).sum()
        counts[int(size)] += 1
    print(json.dumps(dict(section_size_histogram=counts)), flush=True)

    source = load_points()
    source_equator = source[np.all(source[:, 24:] == 0, axis=1), :24]
    missing_roots = shell[~np.isin(keys(shell), keys(source_equator), assume_unique=True)]
    missing_roots = np.pad(missing_roots, ((0, 0), (0, 3)))
    affected_roots = deleted[np.isin(keys(deleted), keys(source_equator), assume_unique=True)]
    affected = np.r_[source[np.any(source[:, 24:] != 0, axis=1)],
                     np.pad(affected_roots, ((0, 0), (0, 3)))]
    sections27 = []
    for frame in aligned_sections():
        projection = tails @ frame
        keep = np.sum(projection ** 2, axis=1) == 4 * np.sum(tails ** 2, axis=1)
        assert np.all(projection[keep] % 2 == 0)
        layer = np.c_[auxiliary[keep[index], :24], projection[index[keep[index]]] // 2]
        sections27.append(layer)
    jobs = list(product(range(12), permutations(range(3)), product((-1, 1), repeat=3)))
    with ProcessPoolExecutor(args.workers, mp_context=get_context("fork")) as pool:
        results = list(pool.map(check_alignment, jobs[:args.alignments]))
    print(json.dumps(dict(alignments=len(results), matching_certificates=sum(c == 0 for _, c in results),
                          components_checked=sum(c for _, c in results),
                          largest_component_gain=max((g for g, _ in results), default=0))))


if __name__ == "__main__":
    main()
