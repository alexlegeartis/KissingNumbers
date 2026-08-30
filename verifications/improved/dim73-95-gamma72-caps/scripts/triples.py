#!/usr/bin/env python3
"""Zero-sum triple partitions of the root systems A_1, A_2, A_3, D_4, D_5, E_6, E_7, E_8.

In the cap construction the gain is  2 * c * (|W| - M),  where W is a kissing configuration
of R^k used for the cap directions and M is the number of parts in a partition of W into
sets that are internally at cosine <= -1/2.  n points pairwise at cosine <= -1/2 satisfy
n <= 1 + 1/(1/2) = 3, so the best possible part is a ZERO-SUM TRIPLE  u + v + w = 0  (three
unit vectors pairwise at 120 degrees), and M >= ceil(|W|/3).

This script constructs each root system, finds a perfect partition into zero-sum triples,
and verifies it.  A perfect partition attains M = |W|/3 and hence the factor 2|W|/3.

Cross-check: with the Leech class of 248 lines these factors reproduce Cohn's table for
dimensions 26, 28, 29, 30, 31 exactly, and give dimension 27 = 200540 (better than the
tabulated 200044, which uses 5 parts where 4 suffice).
"""
import itertools
import numpy as np


def roots(name):
    """root systems normalised to squared norm 2; pairwise inner products in {0,+-1,+-2}."""
    if name == "A_1":
        R = [[1.0], [-1.0]]
        return np.array(R) * np.sqrt(2)
    if name == "A_2":
        a = np.arange(6) * np.pi / 3
        return np.sqrt(2) * np.stack([np.cos(a), np.sin(a)], 1)
    if name in ("A_3", "D_4", "D_5"):
        n = {"A_3": 3, "D_4": 4, "D_5": 5}[name]
        R = []
        for i, j in itertools.combinations(range(n), 2):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [0.0] * n
                    v[i] = si
                    v[j] = sj
                    R.append(v)
        return np.array(R)
    if name == "E_8":
        R = []
        for i, j in itertools.combinations(range(8), 2):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [0.0] * 8
                    v[i] = si
                    v[j] = sj
                    R.append(v)
        for m in range(256):
            s = [(1 if (m >> b) & 1 == 0 else -1) for b in range(8)]
            if s.count(-1) % 2 == 0:
                R.append([x * 0.5 for x in s])
        return np.array(R)
    if name in ("E_7", "E_6"):
        E8 = roots("E_8")
        u = E8[0]
        if name == "E_7":                     # E_8 roots orthogonal to one root
            return E8[np.abs(E8 @ u) < 1e-9]
        # E_6 is the orthogonal complement of an A_2, not of two orthogonal roots
        # (the latter gives D_6, 60 roots).  Pick v with <u,v> = 1 in norm-2 scaling.
        v = E8[np.abs(E8 @ u - 1.0) < 1e-9][0]
        return E8[(np.abs(E8 @ u) < 1e-9) & (np.abs(E8 @ v) < 1e-9)]
    raise ValueError(name)


def triple_partition(R, seed=0):
    """greedy + restart search for a perfect partition into zero-sum triples."""
    n = len(R)
    rng = np.random.default_rng(seed)

    def find(w):
        """index of the vector equal to w, by nearest neighbour with a tolerance."""
        d = np.abs(R - w).max(axis=1)
        i = int(d.argmin())
        return i if d[i] < 1e-9 else None
    for _ in range(300):
        used = np.zeros(n, bool)
        parts = []
        order = rng.permutation(n)
        ok = True
        for i in order:
            if used[i]:
                continue
            cand = [j for j in range(n) if not used[j] and j != i
                    and abs(R[i] @ R[j] + 0.5) < 1e-9]   # 120 degrees, unit vectors
            rng.shuffle(cand)
            placed = False
            for j in cand:
                k = find(-(R[i] + R[j]))
                if k is not None and not used[k] and k not in (i, j):
                    used[[i, j, k]] = True
                    parts.append((i, j, k))
                    placed = True
                    break
            if not placed:
                ok = False
                break
        if ok and used.all():
            return parts, [], []
    # |W| need not be divisible by 3 (D_5 has 40 roots).  Take as many triples as possible,
    # then pair the leftovers at cosine <= -1/2, then singletons.  The gain factor is
    # |W| - M = 2*(triples) + (pairs), so triples are worth twice a pair.
    # enumerate every zero-sum triple once, then pack them MOST-CONSTRAINED-FIRST:
    # repeatedly take the unused vertex lying in the fewest still-available triples.
    alltri = []
    for i in range(n):
        for j in range(i + 1, n):
            if abs(R[i] @ R[j] + 0.5) > 1e-9:
                continue
            k = find(-(R[i] + R[j]))
            if k is not None and k > j:
                alltri.append((i, j, k))
    inc = [[] for _ in range(n)]
    for t, (i, j, k) in enumerate(alltri):
        inc[i].append(t); inc[j].append(t); inc[k].append(t)
    best = None
    for _ in range(400):
        used = np.zeros(n, bool)
        dead = np.zeros(len(alltri), bool)
        tri = []
        while True:
            free = [v for v in range(n) if not used[v]]
            if not free:
                break
            opts = {v: [t for t in inc[v] if not dead[t]] for v in free}
            live = [v for v in free if opts[v]]
            if not live:
                break
            v = min(live, key=lambda x: (len(opts[x]), rng.random()))
            t = opts[v][rng.integers(len(opts[v]))]
            i, j, k = alltri[t]
            used[[i, j, k]] = True
            tri.append((i, j, k))
            for w in (i, j, k):
                for tt in inc[w]:
                    dead[tt] = True
        rest = [i for i in range(n) if not used[i]]
        pairs = []
        for i in list(rest):
            if used[i]:
                continue
            for j in rest:
                if j != i and not used[j] and not used[i] and R[i] @ R[j] <= -0.5 + 1e-9:
                    used[[i, j]] = True
                    pairs.append((i, j))
                    break
        singles = [i for i in range(n) if not used[i]]
        cand = (tri, pairs, singles)
        val = 2 * len(tri) + len(pairs)
        if best is None or val > best[0]:
            best = (val, cand)
    return best[1]


TAU = {"A_1": 2, "A_2": 6, "A_3": 12, "D_4": 24, "D_5": 40, "E_6": 72, "E_7": 126, "E_8": 240}
K = {"A_1": 1, "A_2": 2, "A_3": 3, "D_4": 4, "D_5": 5, "E_6": 6, "E_7": 7, "E_8": 8}

print(__doc__)
print("=" * 78)
print("  k  system  |W|   check      partition            M    factor |W|-M")
res = {}
for name in ("A_1", "A_2", "A_3", "D_4", "D_5", "E_6", "E_7", "E_8"):
    R = roots(name)
    R = R / np.linalg.norm(R, axis=1)[:, None]          # unit vectors
    G = R @ R.T
    np.fill_diagonal(G, -1.0)
    assert len(R) == TAU[name], (name, len(R))
    assert G.max() <= 0.5 + 1e-9, (name, G.max())        # is a kissing configuration
    tri, pairs, singles = triple_partition(R)
    for (i, j, l) in tri:                                # verify every triple exactly
        assert np.abs(R[i] + R[j] + R[l]).max() < 1e-9
        assert max(R[i] @ R[j], R[i] @ R[l], R[j] @ R[l]) <= -0.5 + 1e-9
    for (i, j) in pairs:                                 # verify every pair
        assert R[i] @ R[j] <= -0.5 + 1e-9
    assert 3 * len(tri) + 2 * len(pairs) + len(singles) == len(R)
    assert len(set([x for t in tri for x in t] + [x for q in pairs for x in q]
                   + singles)) == len(R)                 # a genuine partition
    M = len(tri) + len(pairs) + len(singles)
    note = "%d triples + %d pairs + %d singles" % (len(tri), len(pairs), len(singles))
    f = TAU[name] - M
    res[K[name]] = f
    print("  %d  %-6s %4d  kissing OK  %-22s %4d   %d"
          % (K[name], name, TAU[name], note, M, f))

print()
print("=" * 78)
print("CROSS-CHECK against Cohn's table, dimensions 25-31, with the Leech class c = 248")
print("=" * 78)
PUB = {25: 197056, 26: 198550, 27: 200044, 28: 204520, 29: 209496, 30: 220440, 31: 238350}
for k in range(1, 8):
    tot = 196560 + 2 * 248 * res[k] + TAU[{1: "A_1", 2: "A_2", 3: "A_3", 4: "D_4",
                                           5: "D_5", 6: "E_6", 7: "E_7"}[k]]
    p = PUB[24 + k]
    print("   dim %2d : %d   table %d   %s" % (24 + k, tot, p,
          "MATCH" if tot == p else "BEATS by %d" % (tot - p)))
print()
print("factors for k = 1..8 :", [res[k] for k in range(1, 9)])
print("previous (lambda(k) = tau(k)/2):", [TAU[n] // 2 for n in
      ("A_1", "A_2", "A_3", "D_4", "D_5", "E_6", "E_7", "E_8")])
