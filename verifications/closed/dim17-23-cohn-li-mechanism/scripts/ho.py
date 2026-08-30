#!/usr/bin/env python3
"""Dimension 19: pushing Ho's construction.

Boon Suan Ho (arXiv:2603.10425) proved tau(19) >= 11948 by the Cohn-Li odd-sign construction,
for which it suffices to exhibit a binary code of length 19 and minimum distance 5 inside the
5-punctured extended binary Golay code.  He constructed one of size 1280, and
tau(19) = 10668 + (code size), so 10668 + 1280 = 11948.  My own LP ceiling on that code is
1536, so up to +256 remains.

Since M_24 is 5-transitive every 5-subset of coordinates is equivalent, so puncturing the
first five is without loss of generality.  The punctured code is linear of dimension 12, and
"pairwise distance >= 5" is a difference condition, so this is exactly

    MAXIMUM INDEPENDENT SET in the Cayley graph on F_2^12 whose connection set S is the
    set of nonzero codewords of punctured weight <= 4.

|S| = 21 (one of weight 3, twenty of weight 4), so the graph is 21-regular on 4096 vertices.

Two things follow from the group structure and are exploited here:

  * the eigenvalues are the character sums sum_{s in S} (-1)^{<chi,s>}, computable exactly,
    giving the Hoffman ratio bound on the independence number;
  * any subgroup H with H n S = empty is an independent set, and a union of cosets of H is
    independent exactly when the corresponding set in the quotient F_2^12 / H is independent
    for the projected connection set.  So  alpha >= |H| * alpha(quotient), and Ho's
    1280 = 256 * 5 is of this shape.

    python ho.py
"""
import sys
import os
import numpy as np
import itertools

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'common'))
from capcheck import golay          # noqa: E402


def build():
    G = np.array(sorted(golay()), dtype=np.uint8)
    keep = list(range(5, 24))
    Q = G[:, keep]
    w = Q.sum(1)
    # index the 4096 codewords by a 12-bit label via a basis
    basis = []
    seen = {0}
    lab = {}
    pack = lambda v: int("".join(map(str, v)), 2)
    idx = {pack(Q[i]): i for i in range(len(Q))}
    return Q, w, idx, pack


Q, W, IDX, pack = build()
n = len(Q)
codes = np.array([pack(q) for q in Q], dtype=np.int64)
pos = {int(c): i for i, c in enumerate(codes)}
S = [int(codes[i]) for i in range(n) if 0 < W[i] <= 4]
print(__doc__)
print("codewords %d, connection set |S| = %d, weights %s"
      % (n, len(S), sorted(set(int(W[pos[s]]) for s in S))))

# ---------------------------------------------------------------- Hoffman ratio bound
# a basis for the code, so characters can be indexed by F_2^12
basis = []
span = {0}
for c in codes:
    c = int(c)
    if c not in span:
        basis.append(c)
        span |= {x ^ c for x in span}
    if len(basis) == 12:
        break
assert len(span) == 4096
# coordinates of each codeword in that basis
coord = {}
for m in range(4096):
    v = 0
    for b in range(12):
        if m >> b & 1:
            v ^= basis[b]
    coord[v] = m
Sc = [coord[s] for s in S]
lam = np.zeros(4096)
for chi in range(4096):
    lam[chi] = sum(-1.0 if bin(chi & s).count("1") % 2 else 1.0 for s in Sc)
d = len(S)
lmin = lam.min()
hoff = n * (-lmin) / (d - lmin)
print()
print("spectrum: degree %d, lambda_min = %.6f, distinct eigenvalues %s"
      % (d, lmin, sorted(set(np.round(lam, 6)))[:12])[:200])
print("Hoffman ratio bound: alpha <= %d * %.4f / (%d + %.4f) = %.2f -> %d"
      % (n, -lmin, d, -lmin, hoff, int(np.floor(hoff))))
print("   (Ho's construction gives 1280; my LP ceiling was 1536)")

# ---------------------------------------------------------------- coset structure
print()
print("subgroup + coset search:  alpha >= |H| * alpha(quotient)")
Sset = set(Sc)


def rand_subgroup(rng, dim):
    """random dim-dimensional subgroup of F_2^12 avoiding S, or None."""
    B = []
    span = {0}
    tries = 0
    while len(B) < dim and tries < 4000:
        tries += 1
        c = int(rng.integers(1, 4096))
        if c in span:
            continue
        new = {x ^ c for x in span}
        if new & Sset:
            continue
        B.append(c)
        span |= new
    return (B, span) if len(B) == dim else None


def quotient_alpha(span, B):
    """exact max independent set of the quotient graph (small), by brute force / greedy."""
    reps = {}
    for v in range(4096):
        key = min(v ^ x for x in span)
        reps.setdefault(key, []).append(v)
    keys = sorted(reps)
    m = len(keys)
    kid = {k: i for i, k in enumerate(keys)}
    adj = [[False] * m for _ in range(m)]
    for i, a in enumerate(keys):
        for j, b in enumerate(keys):
            if i < j and any((a ^ b ^ x) in Sset for x in span):
                adj[i][j] = adj[j][i] = True
    best = []
    if m <= 20:                                   # exact
        for r in range(m, 0, -1):
            found = None
            for comb in itertools.combinations(range(m), r):
                if all(not adj[a][b] for a, b in itertools.combinations(comb, 2)):
                    found = comb
                    break
            if found:
                best = list(found)
                break
    else:                                         # greedy with restarts
        rng2 = np.random.default_rng(0)
        for _ in range(3000):
            cur = []
            for v in rng2.permutation(m):
                if all(not adj[v][u] for u in cur):
                    cur.append(int(v))
            if len(cur) > len(best):
                best = cur
    return len(best), m


rng = np.random.default_rng(0)
overall = 0
for dim in (10, 9, 8, 7, 6):
    bestv = 0
    for _ in range(300):
        r = rand_subgroup(rng, dim)
        if r is None:
            continue
        B, span = r
        a, m = quotient_alpha(span, span)
        bestv = max(bestv, (1 << dim) * a)
    if bestv:
        print("   |H| = 2^%-2d  quotient size %4d   best  %d * a = %d"
              % (dim, 4096 >> dim, 1 << dim, bestv))
        overall = max(overall, bestv)
print()
print("best coset-structured independent set found: %d   (Ho 1280, ceiling %d)"
      % (overall, int(np.floor(hoff))))
