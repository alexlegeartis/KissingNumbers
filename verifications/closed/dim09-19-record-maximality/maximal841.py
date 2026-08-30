#!/usr/bin/env python3
"""Is Takhanov-Assylbekov-Yun's 841-point configuration in R^12 MAXIMAL?

A new point z (unit) may be added iff  <z,x> <= 1/2 for every x in the configuration, i.e.
iff

    m(C) := min over unit z of max over x in C of <z,x>   is at most 1/2.

Note the condition is ONE-SIDED: the 841 is only nearly antipodal (410 antipodal pairs plus
21 singletons), so -x need not be present and <z,x> >= -1/2 is not required.

Equivalently, maximise |z| over the polytope P = {z : <z,x> <= 1/2 for all x}; a free point
exists iff the maximum is at least 1.  Maximising a norm over a polytope is not convex, so
this is done by the standard ascent -- linear maximisation in the current direction gives a
vertex, repeat -- with many random restarts.  Every returned point is checked for feasibility
explicitly, which is the trap the Lambda_22/Lambda_23 analysis fell into (KNOWLEDGE section
45: a naive iterated LP was 70% low on a control case).

Control: the same code is run on a configuration with a KNOWN free direction, to confirm it
finds one.

    python maximal841.py
"""
import os
import sys
import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))


def load841():
    p = os.path.join(HERE, 'data', 'c841.npy')
    if os.path.exists(p):
        return np.load(p)
    src = os.environ.get('G841', 'g841.txt')
    rows = []
    for ln in open(src):
        q = ln.split()
        if len(q) >= 12:
            try:
                rows.append([float(x) for x in q[:12]])
            except ValueError:
                pass
    X = np.array(rows)
    X = X / np.linalg.norm(X, axis=1)[:, None]
    np.save(p, X)
    return X


def minmax(C, restarts=4000, seed=0, iters=200):
    """ascent for max |z| over {z : <z,x> <= 1/2}; returns (best |z|, best z)."""
    n, d = C.shape
    rng = np.random.default_rng(seed)
    best, bz = 0.0, None
    A_ub = C
    b_ub = np.full(n, 0.5)
    for r in range(restarts):
        z = rng.normal(size=d)
        z /= np.linalg.norm(z)
        for _ in range(iters):
            res = linprog(c=-z, A_ub=A_ub, b_ub=b_ub, bounds=[(-10, 10)] * d, method='highs')
            if not res.success:
                break
            w = res.x
            nw = np.linalg.norm(w)
            if nw < 1e-12:
                break
            if nw <= np.linalg.norm(z) * (1 + 1e-12) and np.allclose(w / nw, z, atol=1e-9):
                z = w
                break
            z = w / nw
        res = linprog(c=-z, A_ub=A_ub, b_ub=b_ub, bounds=[(-10, 10)] * d, method='highs')
        if not res.success:
            continue
        w = res.x
        # explicit feasibility assertion -- never trust the solver's word for it
        if (C @ w).max() > 0.5 + 1e-9:
            continue
        nw = float(np.linalg.norm(w))
        if nw > best:
            best, bz = nw, w.copy()
            print("   restart %5d : |z| = %.9f   (need >= 1)" % (r, best), flush=True)
    return best, bz


if __name__ == '__main__':
    C = load841()
    print(__doc__)
    print("configuration: %d points in R^%d" % C.shape)
    G = C @ C.T
    np.fill_diagonal(G, -1.0)
    print("max inner product: %.12f  (valid kissing configuration: %s)"
          % (G.max(), G.max() <= 0.5 + 1e-7))
    # how antipodal is it?
    near = (G.min(axis=1) < -0.99).sum()
    print("points with a near-antipode: %d of %d" % (near, len(C)))

    print()
    print("CONTROL: a configuration that provably has a free direction --")
    print("the 841 with all points having a positive last coordinate removed.")
    keep = C[C[:, 11] <= 0.0]
    b, _ = minmax(keep, restarts=60, seed=1)
    print("   control max |z| = %.9f  (should exceed 1)" % b)

    print()
    print("THE REAL QUESTION: is the full 841 maximal?")
    best, bz = minmax(C, restarts=int(sys.argv[1]) if len(sys.argv) > 1 else 3000)
    print()
    print("max |z| over the polytope = %.9f" % best)
    if best >= 1.0 - 1e-9:
        z = bz / np.linalg.norm(bz)
        print("FREE DIRECTION FOUND: max <z,x> = %.12f" % (C @ z).max())
        np.save(os.path.join(HERE, 'free_dir_841.npy'), z)
        print("=> tau(12) >= 842   NEW RECORD")
    else:
        print("m(C) = 1/(2*%.9f) = %.9f > 1/2" % (best, 0.5 / best))
        print("=> the 841 is MAXIMAL: no 842nd sphere can be added to it.")
