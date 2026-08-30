#!/usr/bin/env python3
"""The antipode construction on the LEECH lattice, for low dimensions.

Sun-Wang (arXiv:2607.20359) / Chen et al. (2025): for a lattice L of minimal norm mu and a
rank-k cross-section K with Gram K, put M = K^{-1} and choose S in M with maximal pairwise
squared distance delta < mu.  This gives a packing in dimension n - k whose kissing number is

    N_i = sum over j != i with |u_i - u_j|^2 = delta of #{ w in L : |w|^2 = mu, U(w) = u_j - u_i }
    kissing = max_i N_i.

For s = 1 this is the bare cross-section (the laminated lattice).  The point of the antipode
is that a rank-k simplex has up to k partners at the maximal distance, so the kissing number
becomes a SUM OF k FIBRES rather than the single zero fibre.

Taking L = Leech (mu = 4, 196560 minimal vectors) and k = 24 - n puts a lot of fibres in play
in low dimensions: k = 8 for n = 16, k = 11 for n = 13.  Every fibre is computed here by
direct enumeration of the Leech's minimal vectors, so nothing is estimated.

With the dual basis, U(w) = sum_l <u_l, w> m_l, so the fibre over the M-point with integer
coordinates c is exactly  #{w minimal : <u_l, w> = c_l for every l}.

    python antipode_leech.py [kmin kmax]
"""
import os
import sys
import itertools
import numpy as np
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from capcheck import leech          # noqa: E402

REC = {12: 841, 13: 1154, 14: 1932, 15: 2564, 16: 4320, 17: 5730, 18: 7654,
       19: 11948, 20: 19448, 21: 29768, 22: 49896, 23: 93150}


def inv_rational(G):
    n = len(G)
    A = [[F(G[i][j]) for j in range(n)] + [F(1 if i == j else 0) for j in range(n)]
         for i in range(n)]
    for i in range(n):
        p = next(r for r in range(i, n) if A[r][i] != 0)
        A[i], A[p] = A[p], A[i]
        pv = A[i][i]
        A[i] = [x / pv for x in A[i]]
        for r in range(n):
            if r != i and A[r][i] != 0:
                f = A[r][i]
                A[r] = [x - f * y for x, y in zip(A[r], A[i])]
    return [row[n:] for row in A]


def chain(A, k, ip=-1):
    """k minimal vectors with consecutive inner product `ip` and the rest 0 (an A_k-type chain)."""
    out = [A[0]]
    while len(out) < k:
        last = out[-1]
        cand = np.nonzero(A @ last == ip)[0]
        found = None
        for c in cand:
            v = A[c]
            if all(int(v @ w) == (ip if w is out[-1] else 0) for w in out):
                found = v
                break
        if found is None:
            return None
        out.append(found)
    return np.array(out)


def main():
    kmin, kmax = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (8, 12)
    A32 = leech().astype(np.int64)              # squared norm 32
    print(__doc__)
    print("Leech minimal vectors: %d" % len(A32))
    # rescale to minimal norm 4: inner products become 0, +-1, +-2, +-4
    assert (np.abs(A32 @ A32[0]) % 8 == 0).all()
    A = A32                                      # keep integers; divide inner products by 8
    for k in range(kmin, kmax + 1):
        n = 24 - k
        U = chain(A, k, ip=-8)                   # -8 at norm 32  ==  -1 at norm 4
        if U is None:
            print("k=%2d: no A_k chain found" % k)
            continue
        Gr = [[int(U[i] @ U[j]) // 8 for j in range(k)] for i in range(k)]   # norm-4 Gram
        Minv = inv_rational(Gr)
        # fibres over the dual basis vectors: <u_l, w> = delta_{l,i}
        IP = (A @ U.T) // 8                      # integer inner products at norm 4
        fib = []
        for i in range(k):
            tgt = np.zeros(k, dtype=np.int64)
            tgt[i] = 1
            fib.append(int((IP == tgt).all(1).sum()))
        bare = int((IP == 0).all(1).sum())
        d2 = [Minv[i][i] for i in range(k)]
        delta = max(d2 + [Minv[i][i] + Minv[j][j] - 2 * Minv[i][j]
                          for i in range(k) for j in range(i + 1, k)])
        # N_0 counts the dual-basis points at the maximal distance
        N0 = sum(fib[i] for i in range(k) if d2[i] == delta)
        print("k=%2d -> dim %2d : bare cross-section %6d ; delta = %s (%s mu=4) ; "
              "fibres %s ; N_0 = %d ; record %s"
              % (k, n, bare, delta, "<" if delta < 4 else ">=", fib[:6], N0, REC.get(n)),
              flush=True)


if __name__ == '__main__':
    main()
