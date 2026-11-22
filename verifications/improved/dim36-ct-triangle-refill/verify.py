"""Verify the 36D Coxeter–Todd triangle switches and phase refill."""

from collections import Counter, defaultdict
from dataclasses import dataclass
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "dim36-ct-nonlocal"


@dataclass
class Geometry:
    gram: np.ndarray
    cosets: list
    minima: np.ndarray
    distance: np.ndarray
    negative: np.ndarray
    axial: np.ndarray
    roots: np.ndarray
    cells: np.ndarray
    original: set
    final: set
    inserted: set
    switches: dict
    points: int


def read_geometry():
    source = (PARENT / "certificate.json").read_bytes()
    cert = json.loads((HERE / "certificate.json").read_bytes())
    assert sha256(source).hexdigest() == cert["source_sha256"]
    parent = json.loads(source)
    assert parent["dimension"] == cert["dimension"] == 36
    assert parent["points"] == 502821 and parent["squared_norm"] == 12
    gram = np.array(parent["gram_12"], dtype=np.int64)
    similar = np.array(parent["sublattice_matrix"], dtype=np.int64)
    assert gram.shape == similar.shape == (12, 12)
    assert np.array_equal(gram, gram.T)
    phi = np.eye(12, dtype=np.int64) - similar
    assert np.array_equal(phi @ gram @ phi.T, gram)
    assert np.array_equal(phi @ phi @ phi, np.eye(12, dtype=np.int64))

    spec = importlib.util.spec_from_file_location("ct_parent_geometry", PARENT / "verify.py")
    library = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(library)
    assert abs(library.determinant(similar.tolist())) == 729
    assert np.array_equal(similar @ gram @ similar.T, 3 * gram)
    first, second = parent["label_first"], parent["label_second"]
    assert len(first) == len(second) == 6 and sorted(first + second) == list(range(12))
    assert np.all((similar[:, first] + similar[:, second]) % 3 == 0)
    powers = 3 ** np.arange(6)
    short = list(library.enumerate_short(gram.tolist()))
    rows = np.array([row for row, _ in short] + [(0,) * 12], dtype=np.int64)
    norms = np.array([norm for _, norm in short] + [0], dtype=np.int64)
    assert np.array_equal((rows @ gram * rows).sum(axis=1), norms)
    assert Counter(norms) == {0: 1, 4: 756, 6: 4032, 8: 20412}
    labels = ((rows[:, first] + rows[:, second]) % 3) @ powers
    minima = np.full(729, 100, dtype=np.int64)
    np.minimum.at(minima, labels, norms)
    cosets = [rows[(labels == i) & (norms == minima[i])] for i in range(729)]
    assert Counter(zip(minima, map(len, cosets))) == {(0, 1): 1, (4, 3): 252, (6, 18): 224, (8, 81): 252}
    digits = np.arange(729)[:, None] // powers % 3
    distance = minima[((digits[:, None, :] - digits[None, :, :]) % 3) @ powers]
    negative = (-digits % 3) @ powers
    roots = np.zeros((729, 3, 12), dtype=np.int64)
    for i in np.flatnonzero(minima == 4):
        v = cosets[i][0]
        roots[i] = [v, v @ phi, v @ phi @ phi]
        assert set(map(tuple, roots[i])) == set(map(tuple, cosets[i]))
    axial = rows[norms == 4] @ similar
    assert np.all((axial @ gram * axial).sum(axis=1) == 12)

    cells = np.array(parent["outer_cosets"], dtype=np.int64)
    assert cells.shape == (5100, 3) and np.array_equal(cells[0], [0, 0, 0])
    assert np.all((cells >= 0) & (cells < 729))
    assert len({tuple(x) for x in cells}) == len(cells)
    assert np.all(minima[cells[1:]].sum(axis=1) == 12)
    original = {(*r, a, b, c) for r in cells[np.all(minima[cells] == 4, axis=1)]
                for a in range(3) for b in range(3) for c in range(3)}
    removed, added = map(lambda key: {tuple(row) for row in cert[key]}, ("removed_points", "added_points"))
    assert len(removed) == len(cert["removed_points"]) and len(added) == len(cert["added_points"])
    assert removed <= original and not (added & original)
    assert all(len(row) == 6 and all(type(n) is int for n in row) for row in added | removed)
    assert all(all(0 <= i < 729 and minima[i] == 4 for i in row[:3]) and
               all(0 <= i < 3 for i in row[3:]) for row in added)
    final = (original - removed) | added
    sizes = np.array(list(map(len, cosets)))
    frozen = cells[1:][(cells[1:] == 0).any(axis=1)]
    base_count = 3 * len(axial) + int(sizes[frozen].prod(axis=1).sum())
    assert len(original) == 99117 and base_count == 403704
    assert base_count + len(final) == cert["points"] > parent["points"]

    switches = defaultdict(set)
    for index, row in cert["switched_triangles"]:
        assert type(index) is int and 0 < index < len(cells)
        assert len(row) == 12 and all(type(k) is int for k in row)
        assert sorted(minima[cells[index]]) == [0, 4, 8]
        large_label = int(cells[index][np.argmax(minima[cells[index]])])
        value = tuple(row)
        assert value in set(map(tuple, cosets[large_label]))
        assert value not in switches[index]
        switches[index].add(value)
    return Geometry(gram, cosets, minima, distance, negative,
                    axial, roots, cells, original, final, added, switches, cert["points"])


def pairing(cell, g):
    small = int(np.flatnonzero(g.minima[cell] == 4)[0])
    large = int(np.flatnonzero(g.minima[cell] == 8)[0])
    return small, large, 3 - small - large


def verify_base(g):
    indices = np.flatnonzero((g.cells == 0).any(axis=1))
    entries = [g.cells[i].copy() for i in indices]
    owners = list(map(int, indices))
    for i in sorted(g.switches):
        old = g.cells[i]
        a, b, c = pairing(old, g)
        alt = np.zeros(3, dtype=np.int64)
        alt[b] = old[b]
        alt[c] = g.negative[old[a]]
        entries.append(alt)
        owners.append(i)
    e = np.array(entries)
    assert np.all(g.minima[e[1:]].sum(axis=1) == 12)
    m = sum(g.distance[e[:, k, None], e[:, k]] for k in range(3))
    different = np.not_equal.outer(owners, owners)
    assert np.all(m[different] >= 12)
    # Within an owner, two different large rows differ by a nonzero
    # kernel vector. Their distance squared is at least 12, so their
    # product is at most 2. The small parts of the two placements are
    # orthogonal. One placement is chosen for each large row.
    return (int(np.count_nonzero(different)) // 2, len(entries))


def verify_cross(g, roots, indices):
    v = roots @ g.gram
    # This also covers relocated small triangles against the axial layer.
    assert int((v @ g.axial.T).max()) <= 6
    tested = 0
    for j, cell in enumerate(g.cells[1:], start=1):
        kinds = sorted(g.minima[cell])
        if kinds == [4, 4, 4]:
            continue
        if kinds == [0, 6, 6]:
            a, b = np.flatnonzero(cell)
            first = (v @ g.cosets[cell[a]].T).max(axis=1)
            second = (v @ g.cosets[cell[b]].T).max(axis=1)
            assert np.all(first[indices[:, a]] + second[indices[:, b]] <= 6)
            tested += len(indices)
            continue
        assert kinds == [0, 4, 8]
        a, b, c = pairing(cell, g)
        large = g.cosets[cell[b]]
        moved = g.switches.get(j, set())
        chosen = np.array([tuple(row) in moved for row in large])
        for place, tag, keep in ((a, cell[a], ~chosen), (c, g.negative[cell[a]], chosen)):
            if keep.any():
                first = (v @ g.roots[tag].T).max(axis=1)
                second = (v @ large[keep].T).max(axis=1)
                assert np.all(first[indices[:, place]] + second[indices[:, b]] <= 6)
                tested += len(indices)
    return tested


def verify_root_angles(g, roots, indices):
    table = roots @ g.gram @ roots.T
    assert np.array_equal(np.diag(table), np.full(len(roots), 4))
    off_diagonal = table.copy()
    np.fill_diagonal(off_diagonal, -4)
    assert int(off_diagonal.max()) <= 2
    checked = 0
    for block in range(3):
        for root in np.unique(indices[:, block]):
            subset = indices[indices[:, block] == root]
            a, b = np.triu_indices(len(subset), k=1)
            if block:
                mask = np.all(subset[a, :block] != subset[b, :block], axis=1)
                a, b = a[mask], b[mask]
            assert np.all(table[subset[a], subset[b]].sum(axis=1) <= 6)
            checked += len(a)
    # If the points share no root, all three block products are at most 2.
    return checked


def main():
    assert __debug__, "run without -O"
    g = read_geometry()
    points = np.array(sorted(g.final), dtype=np.int64)
    labels = np.flatnonzero(g.minima == 4)
    roots = g.roots[labels].reshape(-1, 12)
    root_number = np.full(729, -1, dtype=np.int64)
    root_number[labels] = np.arange(len(labels))
    indices = 3 * root_number[points[:, :3]] + points[:, 3:]
    base_pairs, base_cosets = verify_base(g)
    crosses = verify_cross(g, roots, indices)
    root_pairs = verify_root_angles(g, roots, indices)
    print(json.dumps({
        "status": "PASS",
        "dimension": 36,
        "points": g.points,
        "added_root_points": len(g.inserted),
        "removed_root_points": len(g.original - g.final),
        "changed_triangles": sum(map(len, g.switches.values())),
        "affected_base_owners": len(g.switches),
        "base_different_owner_pairs": base_pairs,
        "base_coset_placements": base_cosets,
        "root_point_against_base_checks": crosses,
        "root_point_pairs_with_shared_root": root_pairs,
    }, indent=2))


if __name__ == "__main__":
    main()
