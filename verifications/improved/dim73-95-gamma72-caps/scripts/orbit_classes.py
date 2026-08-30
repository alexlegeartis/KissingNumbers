#!/usr/bin/env python3
"""Turn ONE good class of Gamma_72 into 120 pairwise DISJOINT classes of the same size.

If g is an automorphism of Gamma_72 (an integer matrix with g G g^T = G) and C is a class
-- a set of minimal lines with pairwise |<u,u'>| <= 2 -- then C g is again a class, because
g preserves both the norm and the inner product exactly.  Two random images C g_i, C g_j
overlap in about |C|^2 / 3109087800 lines, i.e. essentially never for |C| ~ 2000, so a
handful of random automorphisms already produces the tau(8)/2 = 120 disjoint classes the
cap construction needs.  Any accidental overlap is simply deleted from the later class, so
the output is disjoint by construction and is verified as such.

    python orbit_classes.py [classfile] [nclasses] [wordlen] [seed]
Writes disjoint_classes.npz with keys c0..c<n-1>.
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(HERE, 'gamma72_gram.npy')):
    raise SystemExit(
        "orbit_classes.py is a PIPELINE STAGE that BUILT the shipped data, not a "
        "verifier.  It\nreads a pool of Gamma_72 minimal vectors and other intermediates that are far\ntoo large to ship.  What it produced IS shipped, and \nscripts/verify_rebuild.py rebuilds and checks the 32000-class family from the small files\nthat are.")


def canon(V):
    nz = V != 0
    first = np.argmax(nz, axis=1)
    s = np.sign(V[np.arange(len(V)), first]).astype(V.dtype)
    return V * s[:, None]

def main():
    cf = sys.argv[1] if len(sys.argv) > 1 else 'class_improved.npy'
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    L = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    G = np.load(os.path.join(HERE, 'gamma72_gram.npy')).astype(np.int64)
    GE = np.load(os.path.join(HERE, 'gamma72_gens.npy')).astype(np.int64)
    C = np.load(os.path.join(HERE, cf)).astype(np.int64)
    print("base class: %d lines from %s" % (len(C), cf))
    A = C @ G @ C.T
    np.fill_diagonal(A, 0)
    assert (np.einsum('ij,jk,ik->i', C, G, C) == 8).all()
    assert int(np.abs(A).max()) <= 2
    print("base class verified: all norms 8, max |off-diagonal ip| =", int(np.abs(A).max()))

    rng = np.random.default_rng(seed)
    r = rng.integers(-(2**62), 2**62, size=72, dtype=np.int64)
    I = np.eye(72, dtype=np.int64)
    seen = np.zeros(0, dtype=np.int64)
    out = {}
    sizes = []
    overlaps = 0
    gmax = 0
    for i in range(n):
        g = I.copy()
        if i > 0:
            for _ in range(L):
                g = g @ GE[rng.integers(0, len(GE))]
        assert (g @ G @ g.T == G).all(), "not an isometry"
        gmax = max(gmax, int(np.abs(g).max()))
        Ci = C @ g
        assert (np.einsum('ij,jk,ik->i', Ci, G, Ci) == 8).all()
        M = Ci @ G @ Ci.T
        np.fill_diagonal(M, 0)
        assert int(np.abs(M).max()) <= 2, "image is not a class"
        K = canon(Ci)
        h = K @ r
        if len(seen):
            dup = np.isin(h, seen)
            if dup.any():
                overlaps += int(dup.sum())
                K, h = K[~dup], h[~dup]
        # also guard against a repeated line inside the image itself (cannot happen)
        _, k = np.unique(h, return_index=True)
        K, h = K[k], h[k]
        seen = np.concatenate([seen, h])
        out['c%d' % i] = K.astype(np.int16)
        sizes.append(len(K))
    S = np.vstack([out[k] for k in out])
    assert len(np.unique(S, axis=0)) == len(S), "classes are not disjoint"
    np.savez_compressed(os.path.join(HERE, 'disjoint_classes.npz'), **out)
    print("automorphism entries stay bounded: max |g_ij| = %d" % gmax)
    print("%d disjoint classes; sizes min %d max %d; total lines %d; overlaps removed %d"
          % (n, min(sizes), max(sizes), sum(sizes), overlaps))
    print("saved disjoint_classes.npz")

if __name__ == '__main__':
    main()
