#!/usr/bin/env python3
"""Explicit coordinate verification of the cap construction that is used for dimensions
   73-80 over Gamma_72.  We cannot write down Gamma_72's 6 218 175 600 minimal vectors, so
   we verify the construction on its exact analogue over the LEECH lattice, in real
   coordinates, for k = 1, 2 and 3.  Every step of the argument -- the deletion rule, the
   within-group condition, the cross-group condition and the poles -- is identical; only
   the numerical values of mu and of the class size differ.

   The construction (ambient R^n (+) R^k, all points unit vectors, all pairwise inner
   products required to be <= 1/2):

     equator   (v/|v|, 0)             for every minimal v of L not in the chosen classes
     caps      (s*sqrt(3)/2 * u/|u|, e/2 * z_i)   s,e in {+-1}, for each line u in class i
     poles     (0, w)                 for w in a kissing configuration of R^k

   where z_1..z_lambda are lines in R^k with pairwise |cos| <= 1/2 (lambda = tau(k)/2) and
   each class i is a set of lines of L with pairwise |<u,u'>| <= mu/3, the classes disjoint.

   Count:  |L-shell| - 2*sum|C_i| + 4*sum|C_i| + tau(k)  =  |L-shell| + 2*sum|C_i| + tau(k).
"""
import numpy as np
import itertools
from fractions import Fraction as F

# ------------------------------------------------------------------ Golay and Leech
def golay():
    g = [1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1]           # generator of the [23,12,7] QR code
    rows = []
    for i in range(12):
        r = [0] * 23
        for j, c in enumerate(g):
            r[(i + j) % 23] = c
        rows.append(r + [0])
    B = []
    for r in rows:
        r = r[:]
        r[23] = sum(r) % 2
        B.append(r)
    B = np.array(B, dtype=np.int8)
    code = set()
    for m in range(1 << 12):
        v = np.zeros(24, dtype=np.int8)
        mm = m
        i = 0
        while mm:
            if mm & 1:
                v ^= B[i]
            mm >>= 1
            i += 1
        code.add(tuple(int(x) for x in v))
    return code

def leech():
    """minimal vectors, squared norm 32; inner products in {0,+-8,+-16,+-32}."""
    C = golay()
    wt = {}
    for c in C:
        wt.setdefault(sum(c), []).append(c)
    assert sorted(wt) == [0, 8, 12, 16, 24] and len(wt[8]) == 759, sorted(
        (k, len(v)) for k, v in wt.items())
    V = []
    for i, j in itertools.combinations(range(24), 2):          # (+-4^2, 0^22)
        for si in (4, -4):
            for sj in (4, -4):
                v = [0] * 24
                v[i] = si
                v[j] = sj
                V.append(v)
    for oc in wt[8]:                                            # (+-2^8, 0^16) on an octad
        pos = [i for i in range(24) if oc[i]]
        for m in range(1 << 7):
            s = [1] * 8
            neg = 0
            for b in range(7):
                if m >> b & 1:
                    s[b] = -1
                    neg += 1
            s[7] = 1 if neg % 2 == 0 else -1                    # even number of minus signs
            v = [0] * 24
            for t, p in enumerate(pos):
                v[p] = 2 * s[t]
            V.append(v)
    for c in C:                                                 # (-3, +-1^23) type
        for k in range(24):
            v = [0] * 24
            for i in range(24):
                v[i] = -1 if c[i] else 1
            v = [x for x in v]
            if v[k] > 0:
                v[k] = -3 if True else v[k]
            else:
                v[k] = 3
            # sign convention: coordinate k carries +-3, the rest +-1 with the
            # negative positions forming the Golay codeword c
            vv = [0] * 24
            for i in range(24):
                vv[i] = -1 if c[i] else 1
            vv[k] = 3 * (-1 if c[k] else 1) * -1
            V.append(vv)
    A = np.array(V, dtype=np.int64)
    A = np.unique(A, axis=0)
    nrm = (A * A).sum(1)
    A = A[nrm == 32]
    return A

# ------------------------------------------------------------------ classes
def strong_classes(A, ip_thr, howmany, target, rng):
    """greedy disjoint sets of LINES with pairwise |<u,u'>| <= ip_thr."""
    n = len(A)
    canon = {}
    reps = []
    for i in range(n):
        t = tuple(int(x) for x in A[i])
        if tuple(-x for x in t) in canon:
            continue
        canon[t] = len(reps)
        reps.append(i)
    L = A[reps]
    used = np.zeros(len(L), dtype=bool)
    out = []
    for _ in range(howmany):
        order = rng.permutation(len(L))
        cls = []
        for i in order:
            if used[i]:
                continue
            if not cls:
                cls.append(i)
                continue
            ips = L[cls] @ L[i]
            if np.all(np.abs(ips) <= ip_thr):
                cls.append(i)
            if len(cls) >= target:
                break
        used[cls] = True
        out.append(L[cls].copy())
    return out

# ------------------------------------------------------------------ R^k pieces
def zlines(k):
    """lines in R^k with pairwise |cos| <= 1/2, and a kissing configuration of R^k."""
    if k == 1:
        Z = np.array([[1.0]])
        W = np.array([[1.0], [-1.0]])
    elif k == 2:
        Z = np.array([[np.cos(a), np.sin(a)] for a in (0, np.pi / 3, 2 * np.pi / 3)])
        W = np.array([[np.cos(a), np.sin(a)] for a in np.arange(6) * np.pi / 3])
    elif k == 3:
        p = (1 + 5 ** 0.5) / 2
        raw = [(0, 1, p), (0, 1, -p), (1, p, 0), (1, -p, 0), (p, 0, 1), (-p, 0, 1)]
        Z = np.array(raw, dtype=float)
        Z /= np.linalg.norm(Z, axis=1)[:, None]
        W = []                                              # D_3 = cuboctahedron, tau(3)=12
        for i, j in itertools.combinations(range(3), 2):
            for si in (1, -1):
                for sj in (1, -1):
                    w = [0.0] * 3
                    w[i] = si
                    w[j] = sj
                    W.append(w)
        W = np.array(W) / 2 ** 0.5
    else:
        raise ValueError(k)
    return Z, W

# ------------------------------------------------------------------ verify
def verify(k, classes, A, mu2=32.0):
    Z, W = zlines(k)
    lam = len(Z)
    assert len(classes) == lam
    zz = np.abs(Z @ Z.T)
    np.fill_diagonal(zz, 0.0)
    assert zz.max() <= 0.5 + 1e-12, ("z lines not at |cos| <= 1/2", zz.max())
    ww = W @ W.T
    np.fill_diagonal(ww, -1.0)
    assert ww.max() <= 0.5 + 1e-12, ("poles not a kissing configuration", ww.max())

    used = np.vstack(classes)
    usedset = set()
    for u in used:
        usedset.add(tuple(int(x) for x in u))
        usedset.add(tuple(-int(x) for x in u))
    keep = np.array([tuple(int(x) for x in v) not in usedset for v in A])
    E = A[keep] / np.sqrt(mu2)                                   # equator

    caps = []
    for i, Ci in enumerate(classes):
        U = Ci / np.sqrt(mu2)
        for s in (1.0, -1.0):
            for e in (1.0, -1.0):
                blk = np.hstack([s * np.sqrt(3) / 2 * U,
                                 np.tile(e * 0.5 * Z[i], (len(U), 1))])
                caps.append(blk)
    caps = np.vstack(caps)
    Efull = np.hstack([E, np.zeros((len(E), k))])
    poles = np.hstack([np.zeros((len(W), 24)), W])
    P = np.vstack([Efull, caps, poles])
    nr = np.linalg.norm(P, axis=1)
    assert np.abs(nr - 1).max() < 1e-9, ("not unit vectors", np.abs(nr - 1).max())

    worst = -2.0
    # equator-equator is the Leech property itself; check the rest exhaustively
    small = np.vstack([caps, poles])
    G = small @ small.T
    np.fill_diagonal(G, -1.0)
    worst = max(worst, G.max())
    step = 20000
    for a in range(0, len(Efull), step):
        blk = Efull[a:a + step] @ small.T
        worst = max(worst, blk.max())
    GE = None
    # equator-equator: verify via integer inner products on the deleted-shell
    ipE = None
    return len(P), worst, len(E), len(caps), len(poles)


if __name__ == '__main__':
    rng = np.random.default_rng(11)
    A = leech()
    print("Leech minimal vectors: %d  (expect 196560)" % len(A))
    ip = A[:200] @ A[0]
    print("inner products from a fixed minimal vector are in %s (expect 0,+-8,+-16,+-32)"
          % sorted(set(int(x) for x in (A @ A[0]))))
    TAU = {1: 2, 2: 6, 3: 12}
    for k in (1, 2, 3):
        lam = TAU[k] // 2
        cls = strong_classes(A, 8, lam, 10 ** 9, rng)            # |<u,u'>| <= 8  <=>  |cos| <= 1/4
        tot, worst, ne, nc, npo = verify(k, cls, A)
        sizes = [len(c) for c in cls]
        pred = 196560 + 2 * sum(sizes) + TAU[k]
        print("k=%d  dim %2d : %d classes of sizes %s" % (k, 24 + k, lam, sizes))
        print("        equator %d + caps %d + poles %d = %d ; predicted %d %s"
              % (ne, nc, npo, tot, pred, "OK" if tot == pred else "MISMATCH"))
        print("        largest inner product among all non-equator-equator pairs: %.15f  %s"
              % (worst, "OK (<= 1/2)" if worst <= 0.5 + 1e-12 else "VIOLATION"))
