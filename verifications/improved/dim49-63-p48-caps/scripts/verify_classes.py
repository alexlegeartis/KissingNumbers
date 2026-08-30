#!/usr/bin/env python3
"""Exact verification of the P_48 classes.

For each saved class C (an (m,48) integer array of coefficient vectors with respect to the
Catalogue-of-Lattices basis of P_48p, Gram matrix G):

  1. every row x satisfies x G x^T = 6  -- it is a minimal vector of P_48;
  2. for every pair x != y in C, |x G y^T| <= 2  -- the class condition (|cos| <= 1/3);
     in particular no pair has |x G y^T| = 6, so no two rows are +- each other and every
     row is a distinct LINE;
  3. the classes are pairwise DISJOINT as sets of lines.

ARITHMETIC.  All quantities are integers and every product below is formed in float64 or
float32 with a magnitude bound far under the exact-integer range of the format, so the BLAS
result is the exact integer:
   |x_i| <= 22 and |G_ij| <= 10, so (X G)_ij = <x_i, e_j> obeys |(XG)_ij| <= sqrt(6*G_jj)
   <= sqrt(60) < 8 (Cauchy-Schwarz, both vectors in the lattice), and the intermediate sums
   in X @ G are at most 48*22*10 = 10560 << 2^53;
   then (XG) X^T has terms at most 8*22 = 176 and sums at most 48*176 = 8448 << 2^24.
The float64/float32 assertions below check those bounds, so nothing is assumed.

DISJOINTNESS is decided by a 64-bit random linear hash of the canonicalised line: equal rows
necessarily have equal hashes, so all hashes distinct proves all rows distinct.

    python verify_classes.py file1.npy|file.npz [file2 ...]
"""
import sys, os, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def canon(V):
    nz = V != 0
    first = np.argmax(nz, axis=1)
    s = np.sign(V[np.arange(len(V)), first]).astype(V.dtype)
    return V * s[:, None]


def check_class(C, G64):
    A = C.astype(np.float64)
    AG = A @ G64                                    # exact: |terms| <= 10560
    assert np.abs(AG).max() < 2 ** 20
    nm = (AG * A).sum(axis=1)
    assert (nm == 6).all(), "row of norm != 6: %s" % sorted(set(nm.tolist()))
    Af = A.astype(np.float32)
    AGf = AG.astype(np.float32)
    assert (AGf == AG).all()
    worst = 0
    step = 8192
    for a in range(0, len(A), step):
        M = AGf[a:a + step] @ Af.T                  # exact: |sums| <= 8448 < 2^24
        for t in range(len(M)):
            M[t, a + t] = 0
        worst = max(worst, int(np.abs(M).max()))
    return worst


def main():
    G = np.load(os.path.join(HERE, '..', 'data', 'p48p_gram.npy')).astype(np.int64)
    assert (G == G.T).all() and (np.diag(G) % 2 == 0).all()
    assert int(np.abs(G).max()) <= 10
    G64 = G.astype(np.float64)
    rng = np.random.default_rng(12345)
    r = rng.integers(-(2 ** 62), 2 ** 62, size=48, dtype=np.int64)
    items = []
    for p in sys.argv[1:]:
        p = p if os.path.isabs(p) else os.path.join(HERE, p)
        if p.endswith('.npz'):
            z = np.load(p)
            ks = sorted(z.files, key=lambda s: int(s[1:]) if s[1:].isdigit() else 0)
            items += [(os.path.basename(p) + ':' + k, z[k]) for k in ks]
        else:
            items.append((os.path.basename(p), np.load(p)))
    t0 = time.time()
    total = 0
    worstall = 0
    hashes = []
    sizes = []
    for name, C in items:
        assert C.shape[1] == 48, (name, C.shape)
        if len(C) == 0:
            sizes.append(0)
            continue
        assert int(np.abs(C).max()) <= 30
        w = check_class(C, G64)
        assert w <= 2, (name, "class condition violated, max |ip| =", w)
        worstall = max(worstall, w)
        total += len(C)
        sizes.append(len(C))
        hashes.append(canon(C.astype(np.int64)) @ r)
        if len(items) <= 12:
            print("%-30s %6d lines   all norms 6   max |off-diagonal ip| = %d  OK"
                  % (name, len(C), w))
    if len(items) > 12:
        print("%d classes verified: every row has norm exactly 6, every off-diagonal "
              "inner product satisfies |ip| <= %d  (%.0fs)"
              % (len(items), worstall, time.time() - t0))
        print("   sizes min %d max %d mean %.1f" % (min(sizes), max(sizes),
                                                    sum(sizes) / len(sizes)))
    if len(hashes) > 1:
        h = np.concatenate(hashes)
        u = len(np.unique(h))
        print("disjointness: %d lines in total, %d distinct hashes -> %s"
              % (len(h), u, "DISJOINT" if u == len(h) else "*** OVERLAP ***"))
        assert u == len(h)
    print("sum of class sizes = %d   (%.0fs)" % (total, time.time() - t0))


if __name__ == '__main__':
    main()
