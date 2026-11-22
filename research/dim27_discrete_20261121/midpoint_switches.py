"""Classify all paired D3 tail assignments on the 13,248 incumbent midpoint bodies."""

import argparse
from collections import deque
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from bw_refill import TAILS
from head_frames import enumerate_frames
from source import HERE, checked_archive, keys, load_points


def source_layers():
    with checked_archive("d27_certificate.npz") as data:
        mid = data["added_integer"][data["midpoint_rows"]]
    mid = mid[mid[np.arange(len(mid)), 24 + np.argmax(mid[:, 24:] != 0, axis=1)] > 0]
    index = {tuple(t): i for i, t in enumerate(TAILS)}
    labels = np.array([index[tuple(t)] for t in mid[:, 24:]])
    order = np.argsort(labels)
    source = load_points()
    heads = source[(source[:, 24:] ** 2).sum(1) == 64]
    return mid[order], labels[order], heads


def graph_components(n, edges):
    matrix = csr_matrix((np.ones(len(edges), dtype=np.int8), edges.T), shape=(n, n))
    return connected_components(matrix, directed=False)[1]


def force_planes(domains, neighbors):
    """Eliminate choices using signed parity constraints between already fixed planes."""
    domains = domains.copy()
    size = len(domains)
    while True:
        block = np.full(size, -1)
        phase = np.zeros(size, dtype=np.int8)
        for start in range(size):
            if domains[start].bit_count() != 1 or block[start] >= 0:
                continue
            block[start] = start
            todo = deque([start])
            while todo:
                u = todo.popleft()
                for v, forbidden in neighbors[u]:
                    if domains[v] != domains[u]:
                        continue
                    assert forbidden != 3
                    opposite = int(forbidden == 1)
                    if block[v] < 0:
                        block[v] = start
                        phase[v] = phase[u] ^ opposite
                        todo.append(v)
                    else:
                        assert phase[v] ^ phase[u] == opposite

        updates = domains.copy()
        for u, domain in enumerate(domains):
            if domain.bit_count() == 1:
                continue
            required = {}
            for v, forbidden in neighbors[u]:
                if block[v] < 0 or not (domain & domains[v]):
                    continue
                group = block[v]
                value = int(phase[v]) ^ int(forbidden == 1)
                if forbidden == 3 or required.get(group, value) != value:
                    updates[u] &= ~domains[v]
                required[group] = value
        assert np.all(updates)
        if np.array_equal(domains, updates):
            return domains
        domains = updates


def switch_groups(mid, labels, heads):
    body = mid[:, :24].astype(np.float32)
    plane = np.argmin(np.abs(TAILS), axis=1)
    half = np.array([sum(plane[:c] == p) for c, p in enumerate(plane)])
    strong, weak = [], []
    for k in range(0, len(mid), 192):
        dot = body[k:k + 192] @ body.T
        i, j = np.where(dot > 32)
        keep = i + k < j
        i, j = i[keep], j[keep]
        val = dot[i, j]
        assert np.all(val <= 64)
        for result, use in [(strong, val > 48), (weak, val <= 48)]:
            result.append(np.column_stack([i[use] + k, j[use]]).astype(np.int32))
    strong, weak = np.vstack(strong), np.vstack(weak)
    p, b = plane[labels], half[labels]
    i, j = strong.T
    assert np.all(p[i] == p[j]) and np.all(b[i] != b[j])
    component = graph_components(len(mid), strong)
    count = int(component.max()) + 1
    incidence = body @ heads[:, :24].astype(np.float32).T
    domains = np.zeros(count, dtype=np.int64)
    for axis in range(3):
        eligible = np.all(incidence[:, heads[:, 24 + axis] == 0] <= 32, axis=1)
        # Heads on the omitted axis pay zero tail cost. All body products must still be <= 64.
        assert np.all(incidence <= 64)
        common = np.ones(count, dtype=bool)
        np.logical_and.at(common, component, eligible)
        domains |= common.astype(np.int64) << axis

    i, j = weak.T
    u, v = np.minimum(component[i], component[j]), np.maximum(component[i], component[j])
    code, inv = np.unique(u.astype(np.int64) * count + v, return_inverse=True)
    forbidden = np.zeros(len(code), dtype=np.int8)
    np.bitwise_or.at(forbidden, inv, 1 << (b[i] ^ b[j]))
    neighbors = [[] for _ in range(count)]
    for pair, mask in zip(code, forbidden):
        u, v = divmod(int(pair), count)
        if u != v:
            neighbors[u].append((v, int(mask)))
            neighbors[v].append((u, int(mask)))
        else:
            assert mask == 2
    forced = force_planes(domains, neighbors)
    assert np.array_equal(forced[component], 1 << p)
    groups = graph_components(len(mid), np.vstack([strong, weak[p[i] == p[j]]]))
    return groups, p, {
        "midpoint_bodies": len(mid), "strong_edges": len(strong), "weak_edges": len(weak),
        "strong_components": count, "initially_plane_pinned": int(sum(x.bit_count() == 1 for x in domains)),
        "forced_components": int(sum(x.bit_count() == 1 for x in forced)), "flip_groups": int(groups.max() + 1),
    }


def contact_spectrum(mid, labels, groups, plane):
    body = mid[:, :24].astype(np.float32)
    count = int(groups.max()) + 1
    changes = np.zeros((count, count), dtype=np.int64)
    for k in range(0, len(mid), 192):
        i, j = np.where(body[k:k + 192] @ body.T == 32)
        i += k
        keep = (i < j) & (plane[i] == plane[j]) & (groups[i] != groups[j])
        i, j = i[keep], j[keep]
        a, b = np.minimum(groups[i], groups[j]), np.maximum(groups[i], groups[j])
        np.add.at(changes, (a, b), np.where(labels[i] == labels[j], -2, 2))

    spectrum, best = {0}, np.zeros(count, dtype=bool)
    for axis in range(3):
        group = np.unique(groups[plane == axis])
        # Flipping all bits in a plane leaves its contact count invariant.
        masks = np.arange(1 << (len(group) - 1))
        bits = (masks[:, None] >> np.arange(len(group))) & 1
        weights = changes[np.ix_(group, group)]
        total = np.zeros(len(masks), dtype=np.int64)
        for i, j in np.argwhere(weights):
            total += weights[i, j] * (bits[:, i] ^ bits[:, j])
        best[group] = bits[np.argmin(total)]
        spectrum = {int(x + y) for x in spectrum for y in np.unique(total)}
    return spectrum, best


def count_contacts(points):
    """Count unordered touching pairs directly in blocks with exact binary32 products."""
    assert points.dtype.kind == "i"
    assert np.all((points.astype(np.int64) ** 2).sum(1) == 128)
    rows = points.astype(np.float32)
    total = 0
    for start in range(0, len(rows), 192):
        block = rows[start:start + 192]
        products = block @ rows[start:].T
        assert np.all(np.triu(products[:, :len(block)], 1) <= 64)
        assert np.all(products[:, len(block):] <= 64)
        total += int(np.count_nonzero(np.triu(products[:, :len(block)] == 64, 1)))
        total += int(np.count_nonzero(products[:, len(block):] == 64))
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--groups", type=Path, help="write phase groups for the contact scan")
    parser.add_argument("--count-contacts", type=Path, help="independently count touching pairs in a coordinate archive")
    args = parser.parse_args()
    if args.count_contacts:
        with np.load(args.count_contacts, allow_pickle=False) as archive:
            print(json.dumps({"contacts": count_contacts(archive["points"])}))
        return
    mid, labels, heads = source_layers()
    groups, plane, stats = switch_groups(mid, labels, heads)
    spectrum, best = contact_spectrum(mid, labels, groups, plane)
    lines = np.argmin(np.abs(TAILS), axis=1)
    by_plane = np.array([np.where(lines == axis)[0] for axis in range(3)])
    if args.groups:
        group_planes = np.empty(int(groups.max()) + 1, dtype=np.int32)
        group_planes[groups] = plane
        np.concatenate([by_plane.ravel(), group_planes, groups]).astype("<i4").tofile(args.groups)
    frames = enumerate_frames()
    head_deltas = {int(x) - frames["source_contacts_incident_to_heads"] for x in frames["contact_histogram"]}
    combined = {x + y for x in spectrum for y in head_deltas}
    if args.out:
        with np.load(HERE / "data/d27_head_frame_tie.npz") as rows:
            base = rows["points"]
        flip = np.empty(6, dtype=np.int32)
        for left, right in by_plane:
            flip[left], flip[right] = right, left
        changed = mid.copy()
        changed[:, 24:] = TAILS[np.where(best[groups], flip[labels], labels)]
        other = changed.copy()
        other[:, 24:] *= -1
        fixed = base[(base[:, 24:] ** 2).sum(1) != 32]
        points = np.vstack([fixed, changed, other])
        assert len(points) == len(base) and len(np.unique(keys(points))) == len(points)
        np.savez_compressed(args.out, points=points, squared_norm=128, inner_product_limit=64,
                            flipped_groups=np.flatnonzero(best))
    stats.update(
        midpoint_assignments=2 ** stats["flip_groups"], frame_choices=frames["frame_triples"],
        midpoint_contact_totals=len(spectrum), combined_contact_totals=len(combined),
        minimum_contact_delta=min(combined), maximum_contact_delta=max(combined),
    )
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
