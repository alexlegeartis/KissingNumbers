#!/usr/bin/env python3
"""Turn ONE good class of P_48 into many pairwise DISJOINT classes.

If g is an automorphism of P_48 (an integer matrix with g G g^T = G) and C is a class -- a
set of minimal lines with pairwise |<u,u'>| <= 2 -- then C g is again a class, because g
preserves the norm and the inner product exactly.  Random images barely overlap: two images
of a class of size c share about c^2 / 26 208 000 lines.  Any overlap is simply deleted from
the later class (a subset of a class is a class), so the output is disjoint by construction
and is verified as such.

For dimension 63 the cap construction needs lambda(15) = tau(15)/2 = 1282 disjoint
classes -- tau, not tau_lat: the record configuration of R^15 is antipodally symmetric.

    python orbit_classes.py [classfile] [nclasses] [wordlen] [seed] [outfile]
"""
import sys, os, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), 'data')
if not os.path.exists(os.path.join(DATA, 'p48p_gram.npy')):
    raise SystemExit(
        "orbit_classes.py is a PIPELINE STAGE that BUILT the shipped data, not a "
        "verifier.  It\nreads a pool of P_48 minimal vectors and other intermediates that are far\ntoo large to ship.  What it produced IS shipped, and \nscripts/regenerate.py rebuilds and checks the class family from the small files\nthat are.")



def canon(V):
    nz = V != 0
    first = np.argmax(nz, axis=1)
    s = np.sign(V[np.arange(len(V)), first]).astype(V.dtype)
    return V * s[:, None]


def main():
    cf = sys.argv[1] if len(sys.argv) > 1 else 'class_p48_7069.npy'
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 1170
    L = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    outf = sys.argv[5] if len(sys.argv) > 5 else 'disjoint_classes_p48.npz'
    G = np.load(os.path.join(DATA, 'p48p_gram.npy')).astype(np.int64)
    GE = np.load(os.path.join(DATA, 'p48p_gens.npy')).astype(np.int64)
    C = np.load(os.path.join(DATA, cf)).astype(np.int64)
    print("base class: %d lines from %s" % (len(C), cf))
    A = C @ G @ C.T
    np.fill_diagonal(A, 0)
    assert (np.einsum('ij,jk,ik->i', C, G, C) == 6).all()
    assert int(np.abs(A).max()) <= 2
    print("base class verified: all norms 6, max |off-diagonal ip| = %d" % int(np.abs(A).max()))
    del A

    rng = np.random.default_rng(seed)
    r = rng.integers(-(2 ** 62), 2 ** 62, size=48, dtype=np.int64)
    I = np.eye(48, dtype=np.int64)
    seen = np.zeros(0, dtype=np.int64)
    seen_sorted = np.zeros(0, dtype=np.int64)
    pending = np.zeros(0, dtype=np.int64)
    gr = rng.integers(-(2 ** 30), 2 ** 30, size=48 * 48, dtype=np.int64)
    gseen = set()
    out = {}
    sizes = []
    overlaps = 0
    gmax = 0
    t0 = time.time()
    for i in range(n):
        # draw a group element not used before: |G| = 145728, so 1170 random words
        # collide by the birthday bound about 5 times, and a repeat gives an empty image
        for _try in range(200):
            g = I.copy()
            if i > 0:
                for _ in range(L):
                    g = g @ GE[rng.integers(0, len(GE))]
            gh = int((g.reshape(-1) * gr).sum())
            if gh not in gseen:
                gseen.add(gh)
                break
        assert (g @ G @ g.T == G).all(), "not an isometry"
        gmax = max(gmax, int(np.abs(g).max()))
        Ci = C @ g
        assert (np.einsum('ij,jk,ik->i', Ci, G, Ci) == 6).all()
        K = canon(Ci)
        h = K @ r
        dup = np.zeros(len(h), dtype=bool)
        if len(seen_sorted):
            pos = np.clip(np.searchsorted(seen_sorted, h), 0, len(seen_sorted) - 1)
            dup |= seen_sorted[pos] == h
        if len(pending):
            dup |= np.isin(h, pending)
        if dup.any():
            overlaps += int(dup.sum())
            K, h = K[~dup], h[~dup]
        _, k = np.unique(h, return_index=True)
        K, h = K[k], h[k]
        pending = np.concatenate([pending, h])
        if len(pending) > 400000:                       # amortise the re-sort
            seen = np.concatenate([seen, pending])
            seen_sorted = np.sort(seen)
            pending = np.zeros(0, dtype=np.int64)
        out['c%d' % i] = K.astype(np.int8)
        sizes.append(len(K))
        if (i + 1) % 100 == 0:
            print("  %4d classes, sizes %d..%d, total %d lines, %d overlaps removed (%.0fs)"
                  % (i + 1, min(sizes), max(sizes), sum(sizes), overlaps, time.time() - t0))
    seen = np.concatenate([seen, pending])
    assert len(np.unique(seen)) == len(seen), "classes are not disjoint"
    np.savez_compressed(os.path.join(DATA, outf), **out)
    print("automorphism entries stay bounded: max |g_ij| = %d" % gmax)
    print("%d disjoint classes; sizes min %d max %d mean %.1f; total lines %d; overlaps removed %d"
          % (n, min(sizes), max(sizes), sum(sizes) / len(sizes), sum(sizes), overlaps))
    print("saved %s" % outf)


if __name__ == '__main__':
    main()
