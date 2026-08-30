#!/usr/bin/env python3
"""Zero-sum triple partitions of the laminated lattices Lambda_9 ... Lambda_23.

For dimensions 81-95 the cap construction (KNOWLEDGE.md sections 50, 51) needs, in R^k, a
kissing configuration W partitioned into parts that are internally at cosine <= -1/2.  A part
has at most 3 elements, and a 3-element part is a ZERO-SUM TRIPLE a + b + c = 0.  The gain
factor is |W| - M where M is the number of parts, so a perfect triple partition gives
2|W|/3 against an antipodal pairing's |W|/2 -- a factor 4/3.

W is a cross-section of the Leech lattice by j = 24 - k minimal vectors forming an A_j root
system (pairwise inner product 16 at squared norm 32).  It lives in R^k and its members are
pairwise at cosine <= 1/2, which is all the construction needs of a direction set.

**IT IS NOT THE MINIMAL SHELL OF Lambda_k, and it is not meant to be.**  The A_j
cross-section coincides with Lambda_(24-j) only for j <= 3:

    j       1       2       3       4       5      ...     15
    A_j     93150   49896   27720   15540   8676   ...     112
    Lambda_ 93150   49896   27720   17400   10668  ...     272
            ^ equal ^               ^ smaller from here on

The `expected` column printed below is the Lambda_k size, shown for exactly that comparison;
from k = 20 down it legitimately differs and the script does not assert equality.  (The right
cross-section for Lambda_20 is by a different root system, not A_4.  See KNOWLEDGE.md
section 53, which corrects an earlier version of this that assumed otherwise.)

WHY A SMALLER W CAN STILL BE THE BETTER CHOICE.  The gain is 2 * sum_i |C_i| * (|Z_i| - 1),
so what matters is |W| - M, not |W|.  A perfect triple partition of a set of size w gives
2w/3; an antipodal pairing of the full Lambda_k gives tau_lat(k)/2.  At k = 19 that is
5708 from the A_5 cross-section against 5334 from Lambda_19 -- the smaller set wins, because
triples are worth 4 points per class line and pairs only 2.  The comparison is made
explicitly, and asserted, in final81-95.py.

A zero-sum triple among vectors of squared norm 32 means pairwise inner product -16, i.e.
120 degrees.  The packing is greedy with most-constrained-first restarts; each step is one
matrix-vector product against the shell, so it stays fast even at 93150 vectors.

    python lamtriples.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'common'))
from capcheck import leech            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXPECT = {23: 93150, 22: 49896, 21: 27720, 20: 17400, 19: 10668, 18: 7398, 17: 5346,
          16: 4320, 15: 2340, 14: 1422, 13: 906, 12: 648, 11: 438, 10: 336, 9: 272}


def a_chain(A, jmax):
    """j minimal vectors pairwise at inner product 16 (an A_j root system), greedily."""
    chain = [A[0]]
    S = A
    out = []
    while len(chain) < jmax:
        ip = S @ chain[-1]
        cand = S[ip == 16]
        for c in cand:
            if all(int(c @ x) == 16 for x in chain):
                chain.append(c)
                break
        else:
            break
        out.append(len(chain))
    return np.array(chain)


def cross(A, chain):
    M = A @ chain.T
    return A[(M == 0).all(axis=1)]


def triple_pack(V, rng, rounds=6):
    """greedy zero-sum triple packing; returns (triples, pairs, singles)."""
    n = len(V)
    key = {}
    for i, v in enumerate(V):
        key[v.tobytes()] = i
    best = None
    for _ in range(rounds):
        used = np.zeros(n, bool)
        tri = []
        for i in rng.permutation(n):
            if used[i]:
                continue
            ip = V @ V[i]
            cand = np.nonzero((ip == -16) & (~used))[0]
            rng.shuffle(cand)
            for j in cand[:400]:
                w = -(V[i] + V[j])
                k = key.get(np.ascontiguousarray(w).tobytes())
                if k is not None and not used[k] and k != i and k != j:
                    used[[i, j, k]] = True
                    tri.append((i, j, k))
                    break
        rest = np.nonzero(~used)[0]
        pairs = []
        for i in rest:
            if used[i]:
                continue
            ip = V[rest] @ V[i]
            for pos, j in enumerate(rest):
                if j != i and not used[j] and ip[pos] <= -16:
                    used[[i, j]] = True
                    pairs.append((i, j))
                    break
        singles = np.nonzero(~used)[0]
        val = 2 * len(tri) + len(pairs)
        if best is None or val > best[0]:
            best = (val, tri, pairs, list(singles))
    return best[1], best[2], best[3]


if __name__ == '__main__':
    rng = np.random.default_rng(3)
    A = leech()
    print(__doc__)
    print("Leech minimal vectors:", len(A))
    chain = a_chain(A, 15)
    print("A_j chain length built:", len(chain))
    G = chain @ chain.T
    assert (np.diag(G) == 32).all()
    off = G[~np.eye(len(chain), dtype=bool)]
    print("chain Gram: diagonal 32, off-diagonal all 16:", bool((off == 16).all()))
    print()
    print("   k   |W|    expected   triples  pairs  singles   M      factor |W|-M   vs |W|/2")
    res = {}
    for j in range(1, len(chain) + 1):
        k = 24 - j
        V = cross(A, chain[:j])
        exp = EXPECT.get(k)
        tri, pairs, singles = triple_pack(V, rng)
        for (a, b, c) in tri[:50]:
            assert (V[a] + V[b] + V[c] == 0).all()
        M = len(tri) + len(pairs) + len(singles)
        assert 3 * len(tri) + 2 * len(pairs) + len(singles) == len(V)
        f = len(V) - M
        res[k] = (len(V), f)
        print("  %2d %6d   %-8s   %6d %5d  %6d  %6d   %6d        %d"
              % (k, len(V), exp if exp else "?", len(tri), len(pairs), len(singles), M, f,
                 len(V) // 2), flush=True)
    np.save(os.path.join(HERE, '..', 'data', 'lam_factors.npy'),
            np.array([[k, res[k][0], res[k][1]] for k in sorted(res)]))
    print()
    print("saved ../data/lam_factors.npy  (k, |W|, |W| - M)")
