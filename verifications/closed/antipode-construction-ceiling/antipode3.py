#!/usr/bin/env python3
"""The antipode construction with the OPTIMAL star, not the dual-basis simplex.

antipode2.py evaluated Chen et al.'s Theorem-2 choice S = {0, m_1, ..., m_k} (the dual basis).
That is one point of a large search space, and a bad one: its delta is close to mu, so the
fibres it sums are the sparse ones.  Re-deriving the constraint from scratch gives a much
better shape.

Let L be a lattice of minimal norm mu, K a rank-k cross-section with Gram K, and M = K^{-1}.
For a finite S in M with maximal pairwise squared distance delta < mu the antipode packing has

    N_i = sum over j with |u_i - u_j|^2 = delta  of  fibre(u_j - u_i),
    fibre(c) = #{ w in L : |w|^2 = mu, <u_l, w> = c_l for all l }.

Only DIAMETER pairs contribute.  Take S = {0, v_1, ..., v_m} with

    |v_i|^2 = delta  for every i        (so every pair (0, v_i) is a diameter pair)
    |v_i - v_j|^2 <= delta              (so delta really is the maximum)

The second condition, given the first, is exactly <v_i, v_j> >= delta/2 -- all the v_i lie in
a 30-degree cap.  Then

    N_0 = sum_i fibre(v_i)

is a sum of m fibres, and m is limited only by how many norm-delta points of M fit in a
30-degree cap.  This is a MAX-WEIGHT CLIQUE over the norm-delta shell of M with weights
fibre(v) -- a star, not a simplex, and it can collect far more than k fibres.

Since M is a lattice, centring the star at 0 is without loss of generality: only the
difference set matters.  Every clique found is a certificate, so the number printed is a
genuine lower bound for that cross-section.

    python antipode3.py [kmin kmax] [trials]
"""
import os
import sys
import itertools
import numpy as np
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from capcheck import leech          # noqa: E402

MU = 4                                            # Leech minimal norm, at the norm-4 scale
REC = {12: 841, 13: 1154, 14: 1932, 15: 2564, 16: 4320, 17: 5730, 18: 7654,
       19: 11948, 20: 19448, 21: 29768, 22: 49896, 23: 93150}
LAM = {12: 648, 13: 918, 14: 1422, 15: 2340, 16: 4320, 17: 5346, 18: 7398,
       19: 10668, 20: 17400, 21: 27720, 22: 49896, 23: 93150}


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


def build_chain(A, k, rng, pool=400):
    """k LINEARLY INDEPENDENT Leech minimal vectors with the largest orthogonal complement.

    The chain vectors are NOT required to be mutually orthogonal -- only the cross-section is
    orthogonal to all of them -- which is what recovers the laminated sections.  They MUST be
    independent, though: a greedy pick from the whole shell happily returns a vector already in
    the span, and the Gram then goes singular from that step on (this silently truncated every
    run past k = 5).  Independence is enforced with a running Gram-Schmidt basis.  All pool
    candidates are scored in one chunked matmul.
    """
    U, alive, basis = [], np.ones(len(A), bool), []

    def independent(v):
        r = v.astype(float).copy()
        for b in basis:
            r -= (r @ b) * b
        nr = np.linalg.norm(r)
        return (nr > 1e-6 * np.linalg.norm(v)), (r / nr if nr > 0 else r)

    for _ in range(k):
        cur = A[np.nonzero(alive)[0]]
        sample = rng.choice(len(A), min(len(A), pool), replace=False)
        C = A[sample]
        scored = []
        for s0 in range(0, len(C), 64):
            blk = C[s0:s0 + 64]
            z = (cur @ blk.T).__eq__(0).sum(0) if len(cur) else np.zeros(len(blk), np.int64)
            for j in range(len(blk)):
                scored.append((int(z[j]), int(sample[s0 + j])))
        scored.sort(key=lambda p: -p[0])
        for _, bi in scored:
            ok, r = independent(A[bi])
            if ok:
                basis.append(r)
                U.append(A[bi])
                alive = alive & (A @ A[bi] == 0)
                break
        else:
            break
    return np.array(U)


def maxweight_clique(nodes, adj, w, rng, restarts=400):
    """Randomised greedy + 1-swap local search.  Returns (weight, member indices)."""
    n = len(nodes)
    if n == 0:
        return 0, []
    order_w = np.argsort(-w)
    best, bestset = 0, []
    for r in range(restarts):
        if r == 0:
            order = order_w
        else:
            noise = rng.random(n) ** 2
            order = np.argsort(-(w * (0.35 + noise)))
        cur, mask = [], np.ones(n, bool)
        for i in order:
            i = int(i)
            if not mask[i]:
                continue
            cur.append(i)
            mask &= adj[i]
            mask[i] = False
        tot = int(w[cur].sum())
        if tot > best:
            best, bestset = tot, list(cur)
    return best, bestset


def main():
    kmin, kmax = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (1, 12)
    trials = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    A = leech().astype(np.int64)
    print(__doc__)
    print("Leech minimal vectors: %d" % len(A))
    best = {}
    for t in range(trials):
        rng = np.random.default_rng(1000 + t)
        U = build_chain(A, kmax, rng)
        print("chain %d built" % t, flush=True)
        for k in range(kmin, kmax + 1):
            Uk = U[:k]
            Gr = [[int(Uk[i] @ Uk[j]) // 8 for j in range(k)] for i in range(k)]
            try:
                Minv = inv_rational(Gr)
            except StopIteration:
                continue
            IP = (A @ Uk.T) // 8
            uniq, cnt = np.unique(IP, axis=0, return_counts=True)
            bare = int(cnt[(uniq == 0).all(1)].sum())
            Mn = np.array([[float(Minv[i][j]) for j in range(k)] for i in range(k)])
            Uf = uniq.astype(float)
            q = np.einsum('ai,ij,aj->a', Uf, Mn, Uf)
            for delta in sorted({round(float(x), 9) for x in q if 1e-9 < x < MU - 1e-9}):
                sel = np.nonzero(np.abs(q - delta) < 1e-7)[0]
                if len(sel) < 2:
                    continue
                V, w = Uf[sel], cnt[sel].astype(np.int64)
                if len(sel) > 700:
                    keep = np.argsort(-w)[:700]
                    V, w = V[keep], w[keep]
                adj = (V @ Mn @ V.T) >= delta / 2 - 1e-7
                np.fill_diagonal(adj, False)
                tot, mem = maxweight_clique(range(len(V)), adj, w, rng, restarts=250)
                rec = best.get(k)
                if rec is None or tot > rec[0]:
                    best[k] = (tot, delta, len(mem), bare, sorted(int(x) for x in w[mem])[-6:])
            print("  k=%2d done: best %s" % (k, best.get(k, (0,))[0]), flush=True)
    print()
    print("  k  dim   bare    delta   star   N_0        Lambda   record")
    for k in sorted(best):
        n = 24 - k
        tot, delta, m, bare, hw = best[k]
        flag = ""
        if tot > REC.get(n, 10 ** 9):
            flag = " *** BEATS THE RECORD ***"
        elif tot > LAM.get(n, 10 ** 9):
            flag = " (beats Lambda_%d)" % n
        print("  %2d  %3d  %6d  %7.4f  %5d   %8d   %6s   %6s%s"
              % (k, n, bare, delta, m, tot, LAM.get(n), REC.get(n), flag))
        print("        heaviest fibres used: %s" % (hw,))


if __name__ == '__main__':
    main()
