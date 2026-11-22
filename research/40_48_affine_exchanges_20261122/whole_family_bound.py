"""Exact ceiling for the 40D complete-family P48 gluing model."""

from collections import Counter
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "verifications/improved/dim41-46-p48-heads/data"
BASE = 1_126_314


def clique_prefix(mask, adjacency, weights):
    """A proper colouring gives an upper bound for every search prefix."""
    result = []
    earlier = 0
    while mask:
        available = mask
        colour_weight = 0
        while available:
            v = available & -available
            i = v.bit_length() - 1
            colour_weight = max(colour_weight, weights[i])
            result.append((i, v, earlier + colour_weight))
            mask ^= v
            available &= ~(v | adjacency[i])
        earlier += colour_weight
    return earlier, result


def heavy_bound(digits, powers, mu, tiny):
    heavy = [(int(c), 16, 116_235 - 40 * int(tiny[c]))
             for c in np.flatnonzero(mu % 6 == 2)]
    medium = [(int(c), 14, 9180) for c in np.flatnonzero(mu % 6 == 4)]
    medium.extend((int(c), 12, 495) for c in np.flatnonzero(mu == 6))
    heavy.sort(key=lambda family: -family[2])
    families = heavy + medium
    classes, radii, weights = np.asarray(families).T
    weights = weights.tolist()

    adjacency = []
    for i, (label, radius, _) in enumerate(families):
        delta = ((digits[label] - digits[classes]) % 3) @ powers
        dot2 = radius + radii - 18 + mu[delta]
        compatible = (dot2 <= 0) | (dot2**2 <= radius * radii)
        compatible[i] = False
        adjacency.append(int.from_bytes(np.packbits(compatible, bitorder="little").tobytes(), "little"))

    heavy_mask = (1 << len(heavy)) - 1
    medium_mask = ((1 << len(families)) - 1) ^ heavy_mask
    target = BASE - int(tiny.sum())
    nodes = 0

    def cover_medium(mask, current, order):
        nonlocal nodes
        nodes += 1
        assert current < target
        for i, bit, cap in reversed(order):
            if current + cap < target:
                return
            remainder = mask & adjacency[i]
            _, prefix = clique_prefix(remainder, adjacency, weights)
            cover_medium(remainder, current + weights[i], prefix)
            mask ^= bit

    def cover_heavy(large, small, current):
        nonlocal nodes
        nodes += 1
        assert current < target
        large_bound, large_order = clique_prefix(large, adjacency, weights)
        small_bound, small_order = clique_prefix(small, adjacency, weights)
        if current + large_bound + small_bound < target:
            return
        if current + small_bound >= target:
            cover_medium(small, current, small_order)
        for i, bit, cap in reversed(large_order):
            if current + cap + small_bound < target:
                return
            cover_heavy(large & adjacency[i], small & adjacency[i], current + weights[i])
            large ^= bit

    # Anchor each nonempty selection at its least heavy index.
    for i in range(len(heavy)):
        later = heavy_mask ^ ((1 << (i + 1)) - 1)
        cover_heavy(adjacency[i] & later, adjacency[i] & medium_mask, weights[i])
    return {"anchors": len(heavy), "search_nodes": nodes, "ceiling": BASE - 1}


def light_bound(digits, powers, mu, tiny):
    def plus(a, b):
        return ((digits[a] + digits[b]) % 3) @ powers

    def minus(a, b):
        return ((digits[a] - digits[b]) % 3) @ powers

    medium = np.flatnonzero(mu == 6)
    large = np.flatnonzero(mu % 6 == 4)
    subspace = np.array([0])
    for generator in (11, 85, 247, 2192):
        subspace = np.unique(np.r_[subspace, plus(subspace, generator), minus(subspace, generator)])
    assert len(subspace) == 81 and np.isin(mu[subspace], [0, 6]).all()
    unseen = (1 << len(mu)) - 1
    groups = []
    while unseen:
        label = (unseen & -unseen).bit_length() - 1
        group = plus(label, subspace)
        assert np.isin(mu[minus(group[:, None], group)], [0, 6]).all()
        for item in group:
            bit = 1 << int(item)
            assert unseen & bit
            unseen ^= bit
        groups.append(group)

    assert len(groups) == 81
    assert np.count_nonzero(mu[groups[0]] == 6) == 80
    assert np.count_nonzero(mu[groups[0]] % 6 == 4) == 0
    assert all(np.count_nonzero(mu[g] == 6) == 27 for g in groups[1:])
    assert len(large) == 2160

    remaining = medium.tolist()
    disjoint_pairs = 0
    while len(remaining) > 1:
        label = remaining.pop()
        incompatible = np.flatnonzero(mu[minus(label, remaining)] == 8)
        if len(incompatible):
            remaining.pop(int(incompatible[0]))
            disjoint_pairs += 1

    # A radius-14 family displaces all radius-12 families in its coset.
    # Two radius-14 families in the same coset also fail the angle test.
    # Independently, at most one family from each matched pair survives.
    loss = 27 * 495 - 9180
    cap = max(min((len(medium) - disjoint_pairs) * 495 + k * 9180,
                  len(medium) * 495 - k * loss) for k in range(81))
    return {"cosets": len(groups), "disjoint_conflicts": disjoint_pairs,
            "ceiling": cap + int(tiny.sum())}


def core_bound(digits, powers, mu, tiny):
    counts = Counter()
    roots = np.flatnonzero(mu == 2)
    # With a radius-18 family, the gluing inequality excludes all
    # radius-14 and radius-16 families by their residue norms. It permits
    # radius-10 and radius-12 families only at differences of norm 0 or 2.
    for c in np.flatnonzero(np.isin(mu, [0, 6])):
        population = 1_092_000 if c == 0 else 1_084_590
        neighbours = np.r_[c, ((digits[c] + digits[roots]) % 3) @ powers]
        medium = neighbours[mu[neighbours] == 6]
        shortest = neighbours[tiny[neighbours] > 0]
        classes = np.r_[c, medium, shortest]
        radii = np.repeat([18, 12, 10], [1, len(medium), len(shortest)])
        delta = ((digits[classes, None] - digits[classes]) % 3) @ powers
        dot2 = radii[:, None] + radii - 18 + mu[delta]
        assert np.all((dot2 <= 0) | (dot2**2 <= radii[:, None] * radii))
        weight = population + 495 * len(medium) + int(tiny[shortest].sum())
        counts[weight] += 1
    return dict(counts)


def main():
    assert __debug__, "Run without Python -O"
    with np.load(DATA / "e8_geometry.npz", allow_pickle=False) as geometry:
        digits = geometry["digits"].astype(np.int64)
        powers = geometry["powers"].astype(np.int64)
        mu = (geometry["mu"] // 4).astype(np.int64)
    with np.load(DATA / "population_table.npz", allow_pickle=False) as table:
        tiny = table["tiny"].astype(np.int64)
    assert sorted(np.unique(mu).tolist()) == [0, 2, 4, 6, 8] and mu[0] == 0
    assert np.count_nonzero(mu == 0) == 1
    assert tiny.sum() == 30_240 and np.all(mu[tiny > 0] % 6 == 2)
    heavy = heavy_bound(digits, powers, mu, tiny)
    light = light_bound(digits, powers, mu, tiny)
    core = core_bound(digits, powers, mu, tiny)
    assert max(core) == BASE > max(light["ceiling"], heavy["ceiling"])
    print(json.dumps({"complete_family_maximum": BASE, "with_radius_18": core,
                      "with_radius_16_no_radius_18": heavy,
                      "radii_10_12_14_only": light}, indent=2))


if __name__ == "__main__":
    main()
