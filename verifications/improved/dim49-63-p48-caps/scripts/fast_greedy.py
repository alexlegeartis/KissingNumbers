#!/usr/bin/env python3
"""Block-greedy maximal independent set in the P_48 minimal-line conflict graph.

A CLASS is a set of minimal LINES of P_48 with pairwise |<u,u'>| <= 2.  P_48 has minimum 6,
so that is exactly |cos| <= 1/3, the condition the cap construction needs (KNOWLEDGE.md
section 50: gamma = 1/3 forces cap level t = 3/4 and antipodal PAIRS of cap directions).
Conflicts are exactly |<u,u'>| = 3.

Greedy over a uniformly random order of the pool, organised so the dominant work is BLAS
matrix-matrix products.

    python fast_greedy.py [nclasses] [seed] [poolfile]
Saves classes_seed<seed>.npz with keys c0, c1, ... (int8 coefficient vectors).
"""
import sys, os, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(HERE, 'p48_gram.npy')):
    raise SystemExit(
        "fast_greedy.py is a PIPELINE STAGE that BUILT the shipped data, not a "
        "verifier.  It\nreads a pool of P_48 minimal vectors and other intermediates that are far\ntoo large to ship.  What it produced IS shipped, and \nscripts/regenerate.py rebuilds and checks the class family from the small files\nthat are.")

THR = 2.5          # integer inner products; |ip| <= 2 accepted, |ip| = 3 rejected


def greedy_once(X, XG, order, blocksize=40000, cchunk=1024):
    chosen = []
    Cx = np.zeros((0, 48), dtype=np.float32)
    for a in range(0, len(order), blocksize):
        blk = order[a:a + blocksize]
        P = XG[blk].astype(np.float32)
        if len(Cx):
            ok = np.ones(len(blk), dtype=bool)
            for c0 in range(0, len(Cx), cchunk):
                idx = np.nonzero(ok)[0]
                if len(idx) == 0:
                    break
                ipc = P[idx] @ Cx[c0:c0 + cchunk].T
                good = (np.abs(ipc) <= THR).all(axis=1)
                ok[idx[~good]] = False
            blk = blk[ok]
            P = P[ok]
        newx = []
        while len(blk):
            i = blk[0]
            chosen.append(int(i))
            v = X[i].astype(np.float32)
            newx.append(v)
            ipv = P @ v
            keep = np.abs(ipv) <= THR
            keep[0] = False
            blk = blk[keep]
            P = P[keep]
        if newx:
            Cx = np.vstack([Cx, np.array(newx, dtype=np.float32)])
    return chosen


def verify(X, G, idx):
    A = X[idx].astype(np.int64)
    worst = 0
    step = 4000
    for a in range(0, len(A), step):
        M = A[a:a + step] @ G @ A.T
        for t in range(len(M)):
            M[t, a + t] = 0
        worst = max(worst, int(np.abs(M).max()))
    dia = np.einsum('ij,jk,ik->i', A, G, A)
    return worst, sorted(set(int(d) for d in dia))


def load_pool(poolfile):
    G = np.load(os.path.join(HERE, 'p48_gram.npy')).astype(np.int64)
    X = np.load(os.path.join(HERE, poolfile))
    N = len(X)
    XG = np.empty((N, 48), dtype=np.int8)
    Gi = G.astype(np.int32)
    for a in range(0, N, 400000):
        B = (X[a:a + 400000].astype(np.int32) @ Gi)
        assert int(np.abs(B).max()) < 127
        XG[a:a + 400000] = B.astype(np.int8)
    return G, X, XG


def main():
    nclasses = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    poolfile = sys.argv[3] if len(sys.argv) > 3 else 'p48_pool.npy'
    G, X, XG = load_pool(poolfile)
    N = len(X)
    print("pool %d lines from %s" % (N, poolfile))
    rng = np.random.default_rng(seed)
    used = np.zeros(N, dtype=bool)
    out = {}
    for c in range(nclasses):
        avail = np.nonzero(~used)[0]
        order = rng.permutation(avail)
        t0 = time.time()
        ch = greedy_once(X, XG, order)
        used[ch] = True
        w, dia = verify(X, G, ch)
        print("class %2d: %5d lines  max|off-diag ip| = %d  norms %s  (%.0fs)"
              % (c, len(ch), w, dia, time.time() - t0))
        assert w <= 2 and dia == [6]
        out['c%d' % c] = X[ch]
        np.savez_compressed(os.path.join(HERE, 'classes_seed%d.npz' % seed), **out)
    print("saved classes_seed%d.npz" % seed)


if __name__ == '__main__':
    main()
