"""Bound all 64 exceptional fibres after releasing their midpoint choices together."""

import json
from collections import Counter

import numpy as np

import fibre_trades as ft
from bw_refill import TAILS
from source import OCTAD


def packed(mask):
    return int.from_bytes(np.packbits(mask, bitorder="little").tobytes(), "little")


def maximal(pools):
    result = []
    for pool in sorted(set(pools), key=int.bit_count, reverse=True):
        if all((pool & other) != pool for other in result):
            result.append(pool)
    return result


def domains(candidates, gates):
    pools = [(1 << len(candidates)) - 1]
    for j in range(32):
        options = maximal([packed((gates[:, state] & (1 << j)) == 0) for state in (0, 1)])
        pools = maximal([p & option for p in pools for option in options])
    return pools


def has_clique(pool, need, edges):
    if need == 0:
        return True
    if pool.bit_count() < need:
        return False
    # Each greedy color class is independent in the compatibility graph.
    remaining, order, upper = pool, [], []
    color = 0
    while remaining:
        color += 1
        part = remaining
        while part:
            bit = part & -part
            i = bit.bit_length() - 1
            order.append(i)
            upper.append(color)
            remaining ^= bit
            part = (part ^ bit) & ~edges[i]
    for i, bound in zip(reversed(order), reversed(upper)):
        if bound < need:
            return False
        if has_clique(pool & edges[i], need - 1, edges):
            return True
        pool ^= 1 << i
    return False


def reaches(points, target):
    gram = points @ points.T
    lookup = {tuple(v): i for i, v in enumerate(points)}
    antipodes = [i for i, v in enumerate(points) if lookup.get(tuple(-v), -1) > i]
    lines = [packed(gram[i, antipodes] == 0) for i in antipodes]
    orthogonal = [packed(gram[i] == 0) for i in antipodes]
    compatible = [packed(row <= 0) for row in gram]
    needed_pairs = max(0, 2 * target - 33)

    def visit(available, remaining, points_left):
        if remaining == 0:
            return has_clique(points_left, target - 2 * needed_pairs, compatible)
        while available.bit_count() >= remaining:
            bit = available & -available
            i = bit.bit_length() - 1
            available ^= bit
            if visit(available & lines[i], remaining - 1, points_left & orthogonal[i]):
                return True
        return False

    return visit((1 << len(lines)) - 1, needed_pairs, (1 << len(points)) - 1)


def run():
    ft.initialize()
    fibres = [32 * c + i for c in (24, 134) for i in range(32)]
    keep = ~np.isin(ft.owner, fibres)
    body = ft.mids[keep, :24].astype(np.float32).T
    old_color, group = ft.labels[keep], ft.groups[keep]
    histogram = Counter()
    for f in fibres:
        r = ft.roots_by_fibre[f]
        a, b = np.where(np.triu(ft.shell[r] @ ft.shell[r].T == 64, 1))
        mid = (ft.shell[r[a]] + ft.shell[r[b]]) // 2
        products = mid.astype(np.float32) @ body
        candidates, gates = [], []
        for p, dot in zip(mid, products):
            if dot.max() > 64:
                continue
            for c, tail in enumerate(TAILS):
                same = ft.plane[old_color] == ft.plane[c]
                if np.any((~same) & (dot > 48)):
                    continue
                close = same & (dot > 32)
                ones = sum(1 << int(j) for j in np.unique(group[close & (old_color == c)]))
                zeros = sum(1 << int(j) for j in np.unique(group[close & (old_color != c)]))
                if zeros & ones:
                    continue
                for sign in (1, -1):
                    candidates.append(np.r_[p[OCTAD], sign * tail])
                    gates.append((ones, zeros))
        for root in r:
            if any(keep[i] for i in ft.incident[root]):
                continue
            candidates.append(np.r_[ft.shell[root, OCTAD], [0, 0, 0]])
            gates.append((0, 0))
        candidates = np.array(candidates, dtype=np.int16)
        count = int(np.sum(~ft.deleted[r]) + 2 * np.sum(ft.owner == f))
        pools = domains(candidates, np.array(gates, dtype=np.int64))
        for pool in pools:
            selected = [i for i in range(len(candidates)) if (pool >> i) & 1]
            if reaches(candidates[selected], count + 1):
                raise ValueError(f"Fibre {f} exceeds its incumbent capacity in the relaxation")
        histogram[(count, len(pools))] += 1
    return {"jointly_released_fibres": len(fibres), "maximal_phase_domains": sum(n * v for (_, n), v in histogram.items()),
            "capacity_by_fibre": {str(k): sum(v for (size, _), v in histogram.items() if size == k)
                                  for k in sorted({size for size, _ in histogram})},
            "total_capacity": sum(k * v for (k, _), v in histogram.items())}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
