#!/usr/bin/env python3
"""The cap directions of dimensions 84 and 85: the Kappa sections K_12 and K_13.

    python kappa_caps.py             re-derive and report against the shipped files
    python kappa_caps.py --write     re-derive and write data/kap1[23]_W.npy and _parts.npz

WHY THIS FILE EXISTS.  The cap layer of dimension 72 + k needs a 60-degree code W of R^k
partitioned into parts that are internally at cosine <= -1/2, and a part of three is a
ZERO-SUM TRIPLE worth 4 points per class line against a pair's 2 (KNOWLEDGE.md section 50).
Bigger W with a perfect triple partition is therefore strictly better, and at k = 12 and 13
the shipped W was NOT the biggest available:

    k     shipped W          |W|      best lattice     |W|     triples
    12    Lambda_12          648      K_12             756     216 -> 252
    13    Lambda_13          906      K_13             918     302 -> 306

Both numbers were already in this repository -- final81-95.py's own ANTIP_FULL table says
tau_lat(12) = 756 and tau_lat(13) = 918, and research/dim81-95/derive.py names them K_12 and
K_13 -- and were used only to CHECK that the triple scheme beat antipodal pairs, never as
configurations.  The reason they were not used is on record in research/dim81-95/richspan2.py:
a greedy hunt for a twelve-space of the Leech holding 756 minimal vectors, which did not find
one.  It was looking for a subspace; the right object is a SYMMETRY.

THE CONSTRUCTION, WITH NO SEARCH IN IT.  Let pi be an automorphism of the Golay code of cycle
shape 1^6 3^6 -- an element of order 3 in M_24 -- acting on the Leech lattice by permuting
coordinates.  It splits R^24 into its fixed space and the orthogonal complement, both of
dimension 12 (six fixed points and six 3-cycles give twelve orbits).  On the COMPLEMENT pi
has no eigenvalue 1, so

    I + pi + pi^2 = 0   there,

and therefore every orbit {v, pi v, pi^2 v} of a vector of the complement SUMS TO ZERO.  The
Leech minimal vectors lying in the complement are exactly those with v + pi v + pi^2 v = 0,
there are 756 of them, they span 12 dimensions, and their zero-sum triple partition is the
orbit set of pi.  Not a packing search: a group acting.

Both halves hold 756.  The complement is the one used, because there the partition is free;
the fixed space would need the same search the laminated sections need.

K_13 IS ONE DIRECTION MORE.  A Leech minimal vector lies in span(K_12) + Ru exactly when the
part of it orthogonal to span(K_12) is parallel to u, so the best thirteenth direction is not
searched for either: project all 196560 minimal vectors onto the complement of the span and
take the largest parallel class.  378 of them tie at 162, and 756 + 162 = 918.  Thirteen is
odd, so no order-3 isometry of R^13 can act freely -- one must fix a vector -- and those last
54 triples ARE searched for, by the same most-constrained-first greedy as lamtriples.py.

IN EXACTLY k COORDINATES.  A rotated copy of the cap set is built by a Cayley transform on
the coordinates W is written in, and that is an isometry of the AXIS SPACE only when W has
full rank in them -- 756 vectors written in 24 coordinates span 12 and the transform would
leave the space.  So both sets are reduced to exactly k coordinates with an integer diagonal
metric, which they admit for a reason: the complement of Fix(pi) is the direct sum, over the
six 3-cycles, of the plane x + y + z = 0, and that plane has the ORTHOGONAL integer basis
(1,-1,0) of norm 2 and (1,1,-2) of norm 6.  Integrality needs x = y mod 2, which holds
because the (-3, 1^23) Leech vectors have every coordinate odd and three odd numbers cannot
sum to zero -- that shape never meets the complement at all.  K_13 adds the one direction u,
orthogonal to the twelve by construction, so its basis is orthogonal outright.

    k = 12    756 vectors, metric diag(1,3)^6, norm 16
    k = 13    918 vectors, metric diag(3,1,3,1,3,1,3,1,3,1,3,1,1), norm 48

WHAT IS CHECKED HERE, AND WHAT IS NAMED.  A count does not identify a lattice.  The
determinant does: at minimal norm 4 these two come out at 729 = 3^6 and 972, which are the
determinants of the Coxeter-Todd lattice K_12 and of K_13 (Conway and Sloane, SPLAG Table
6.1).  Both are asserted below.  scripts/verify_parts.py then re-checks the shipped files
from scratch, without this script, against the conditions the construction actually needs.
"""
import os
import sys
from fractions import Fraction as F
from math import gcd

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'common'))
from golay import golay24                                              # noqa: E402
from m24 import find_autos                                             # noqa: E402
import leech                                                           # noqa: E402

# minimal norm of the Leech vectors in this scaling, and the determinants at minimal norm 4
NORM = 32
DET4 = {12: 729, 13: 972}
WANT = {12: 756, 13: 918}
TRIPLES = {12: 252, 13: 306}


def cycles_of(p):
    """the orbits of a permutation of 24 points, each sorted"""
    seen, out = set(), []
    for i in range(24):
        if i in seen:
            continue
        cyc, j = [], i
        while j not in seen:
            seen.add(j)
            cyc.append(j)
            j = int(p[j])
        out.append(sorted(cyc))
    return out


def order3(seed=1, tries=200, walk=6000):
    """an element of M_24 of cycle shape 1^6 3^6

    find_autos returns whatever the backtracking reaches, which is never of this shape
    directly; products of a handful of them are, and the walk is seeded so it is the same
    element on every run."""
    G = golay24()
    gens = find_autos(G, n_want=6, seed=seed, tries=tries)
    rng = np.random.default_rng(0)
    pool = [np.arange(24)] + list(gens)
    for _ in range(walk):
        c = pool[rng.integers(len(pool))][pool[rng.integers(len(pool))]]
        if len(pool) < 300:
            pool.append(c)
        if sorted(len(o) for o in cycles_of(c)) == [1] * 6 + [3] * 6:
            return np.asarray(c, dtype=np.int64), G
    raise SystemExit('no element of cycle shape 1^6 3^6 was reached')


def hnf_basis(M):
    """a basis of the lattice generated by the integer rows of M"""
    A = [[int(x) for x in r] for r in M]
    rows, cols = len(A), len(A[0])
    r = 0
    for c in range(cols):
        while True:
            nz = [i for i in range(r, rows) if A[i][c]]
            if len(nz) <= 1:
                break
            nz.sort(key=lambda i: abs(A[i][c]))
            p = nz[0]
            A[r], A[p] = A[p], A[r]
            for i in nz[1:]:
                if i == r:
                    continue
                q = A[i][c] // A[r][c]
                if q:
                    A[i] = [x - q * y for x, y in zip(A[i], A[r])]
        nz = [i for i in range(r, rows) if A[i][c]]
        if not nz:
            continue
        A[r], A[nz[0]] = A[nz[0]], A[r]
        r += 1
        if r == rows:
            break
    return [row for row in A[:r] if any(row)]


def det_gram(B):
    """the exact Gram determinant of a list of integer rows"""
    n = len(B)
    M = [[F(sum(B[i][t] * B[j][t] for t in range(len(B[0])))) for j in range(n)]
         for i in range(n)]
    d = F(1)
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c]), None)
        if p is None:
            return F(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for i in range(c + 1, n):
            if M[i][c]:
                f = M[i][c] / M[c][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return d


def check_code(V):
    """60-degree code, exactly: 4<x,y>^2 <= |x|^2|y|^2 whenever <x,y> > 0"""
    G = V @ V.T
    assert (np.diag(G) == NORM).all(), 'not all of one norm'
    off = G[~np.eye(len(G), dtype=bool)]
    assert 4 * int(off.max()) ** 2 <= NORM * NORM, 'not a 60-degree code'
    return G


def blocks(P):
    """the six 3-cycles of pi and its six fixed points"""
    cyc = [o for o in cycles_of(P) if len(o) == 3]
    fixed = [o[0] for o in cycles_of(P) if len(o) == 1]
    assert len(cyc) == 6 and len(fixed) == 6, 'pi is not of shape 1^6 3^6'
    return cyc, fixed


def ortho_basis(P, extra=None):
    """an orthogonal integer basis of span(K_12), and of span(K_13) when u is supplied"""
    cyc, _fixed = blocks(P)
    B = []
    for (a, b, d) in cyc:
        f = np.zeros(24, np.int64)
        f[a], f[b] = 1, -1
        B.append(f)
        g = np.zeros(24, np.int64)
        g[a], g[b], g[d] = 1, 1, -2
        B.append(g)
    if extra is not None:
        B.append(extra)
    B = np.array(B, dtype=np.int64)
    nrm = np.array([int(f @ f) for f in B], dtype=np.int64)
    assert (B @ B.T == np.diag(nrm)).all(), 'the basis is not orthogonal'
    return B, nrm


def reduce_coords(V, B, nrm):
    """(X, c, m): V in the coordinates of the orthogonal basis B, with an integer metric

    Each coordinate is cleared of denominators on its own and then of any square factor in
    its weight, which is what turns a weight of 9 into 1 rather than leaving it there."""
    k = len(B)
    A = [[F(int(w @ B[j]), int(nrm[j])) for j in range(k)] for w in V]
    den = [1] * k
    for row in A:
        for j, v in enumerate(row):
            den[j] = den[j] * v.denominator // gcd(den[j], v.denominator)
    X = np.array([[int(v * den[j]) for j, v in enumerate(row)] for row in A], dtype=np.int64)
    q = [F(int(nrm[j]), den[j] ** 2) for j in range(k)]
    lcm = 1
    for v in q:
        lcm = lcm * v.denominator // gcd(lcm, v.denominator)
    c = [int(v * lcm) for v in q]
    g = 0
    for v in c:
        g = gcd(g, v)
    while g > 1:                                   # a common factor is a scaling, not a metric
        d = next(p for p in (2, 3, 5, 7, 11, 13) if g % p == 0)
        c = [v // d for v in c]
        g //= d
    for j in range(k):                             # a square weight belongs in the coordinate
        for p in (2, 3, 5, 7, 11, 13):
            while c[j] % (p * p) == 0:
                c[j] //= p * p
                X[:, j] *= p
    c = np.array(c, dtype=np.int64)
    G = (X * c) @ X.T
    m = int(G[0, 0])
    assert (np.diag(G) == m).all(), 'the reduced vectors are not of one norm'
    G24 = V @ V.T
    assert (G * int(G24[0, 0]) == G24 * m).all(), 'the metric is not the Euclidean one'
    assert np.linalg.matrix_rank(X.astype(np.float64)) == k, 'the reduction lost rank'
    assert len({tuple(r) for r in X.tolist()}) == len(X), 'the reduction merged vectors'
    off = G[~np.eye(len(G), dtype=bool)]
    assert 4 * int(off.max()) ** 2 <= m * m, 'the reduction broke the 60-degree property'
    return X, c, m


def kappa12(P):
    """the 756 Leech minimal vectors of the complement of Fix(pi), and pi's orbits"""
    A = leech.load().astype(np.int64)
    perp = (A + A[:, P] + A[:, P][:, P] == 0).all(1)
    V = A[perp]
    assert len(V) == WANT[12], 'the complement holds %d minimal vectors, not %d' \
        % (len(V), WANT[12])
    assert np.linalg.matrix_rank(V.astype(np.float64)) == 12, 'the complement is not 12-dim'
    G = check_code(V)
    idx = {tuple(r): i for i, r in enumerate(V.tolist())}
    mp = np.array([idx[tuple(r)] for r in V[:, P].tolist()])
    assert (mp != np.arange(len(V))).all(), 'pi fixes a vector of the complement'
    seen = np.zeros(len(V), bool)
    T = []
    for i in range(len(V)):
        if seen[i]:
            continue
        a, b = int(mp[i]), int(mp[mp[i]])
        assert not (seen[a] or seen[b]), 'the orbits of pi are not disjoint'
        seen[[i, a, b]] = True
        T.append((i, a, b))
    T = np.array(T, dtype=np.int64)
    assert len(T) == TRIPLES[12] and seen.all(), 'the orbits do not cover'
    assert (V[T[:, 0]] + V[T[:, 1]] + V[T[:, 2]] == 0).all(), 'an orbit does not sum to zero'
    for a in range(3):
        for b in range(a + 1, 3):
            assert (G[T[:, a], T[:, b]] == -NORM // 2).all(), 'an orbit is not at cosine -1/2'
    return V, T, A


def kappa13(V12, A):
    """the 918 of the richest 13-space through span(K_12), and a triple partition"""
    B = []
    for r in V12:
        t = B + [r]
        if np.linalg.matrix_rank(np.array(t, dtype=np.float64)) == len(t):
            B.append(r)
        if len(B) == 12:
            break
    B = np.array(B, dtype=np.int64)
    Gb = B @ B.T
    M = [[F(int(Gb[i][j])) for j in range(12)] + [F(int(i == j)) for j in range(12)]
         for i in range(12)]
    for c in range(12):
        p = next(i for i in range(c, 12) if M[i][c])
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for i in range(12):
            if i != c and M[i][c]:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    den = 1
    for i in range(12):
        for j in range(12):
            den = den * M[i][12 + j].denominator // gcd(den, M[i][12 + j].denominator)
    adj = np.array([[int(M[i][12 + j] * den) for j in range(12)] for i in range(12)],
                   dtype=object)
    assert (adj @ Gb.astype(object) == den * np.eye(12, dtype=object)).all(), 'inverse wrong'
    Ao = A.astype(object)
    R = den * Ao - (Ao @ B.T.astype(object)) @ adj @ B.astype(object)
    assert (R @ B.T.astype(object) == 0).all(), 'residuals are not orthogonal to the span'
    inV = np.array([not any(r) for r in R])
    assert int(inV.sum()) == WANT[12], 'the span holds %d, not %d' % (inV.sum(), WANT[12])

    # bucket the rest by direction: divide out the gcd, then fix the sign of the first entry
    own = np.where(~inV)[0]
    keys = {}
    for i in own:
        r = [int(v) for v in R[i]]
        g = 0
        for v in r:
            g = gcd(g, abs(v))
        r = [v // g for v in r]
        j = next(t for t in range(24) if r[t])
        if r[j] < 0:
            r = [-v for v in r]
        keys.setdefault(tuple(r), []).append(int(i))
    top = max(len(v) for v in keys.values())
    best = min(k for k, v in keys.items() if len(v) == top)      # ties broken deterministically
    print('   largest parallel classes tie at %d (%d of them); %d + %d = %d'
          % (top, sum(1 for v in keys.values() if len(v) == top), WANT[12], top,
             WANT[12] + top))
    V = np.vstack([A[inV], A[np.array(keys[best], dtype=np.int64)]])
    assert len(V) == WANT[13], 'the section holds %d, not %d' % (len(V), WANT[13])
    assert np.linalg.matrix_rank(V.astype(np.float64)) == 13, 'the section is not 13-dim'
    G = check_code(V)
    return V, G, np.array(best, dtype=np.int64)


def triple_partition(V, G, target, seeds=600):
    """a perfect zero-sum triple partition, most-constrained-first with restarts"""
    n = len(V)
    idx = {tuple(r): i for i, r in enumerate(V.tolist())}
    tri = []
    for i in range(n):
        for j in np.where(G[i] == -NORM // 2)[0]:
            j = int(j)
            if j <= i:
                continue
            h = idx.get(tuple((-(V[i] + V[j])).tolist()))
            if h is not None and h > j:
                tri.append((i, j, h))
    tri = np.array(tri, dtype=np.int64)
    by = [[] for _ in range(n)]
    for ti, t in enumerate(tri):
        for v in t:
            by[v].append(ti)
    best = []
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        used = np.zeros(n, bool)
        dead = np.zeros(len(tri), bool)
        left = np.array([len(b) for b in by], dtype=np.int64)
        out = []
        while True:
            free = np.where(~used & (left > 0))[0]
            if not len(free):
                break
            pool = free[left[free] == left[free].min()]
            v = int(pool[rng.integers(len(pool))])
            opts = [t for t in by[v] if not dead[t]]
            if not opts:
                left[v] = 0
                continue
            t = int(opts[rng.integers(len(opts))])
            a, b, c = tri[t]
            used[[a, b, c]] = True
            out.append(t)
            for u in (a, b, c):
                for s in by[u]:
                    if not dead[s]:
                        dead[s] = True
                        for x in tri[s]:
                            left[x] -= 1
        if len(out) > len(best):
            best = out
        if 3 * len(best) == n:
            print('   perfect partition at restart %d' % seed)
            break
    assert 3 * len(best) == n, 'only %d triples, %d vectors left over' \
        % (len(best), n - 3 * len(best))
    T = tri[np.array(best, dtype=np.int64)]
    assert len(set(T.ravel().tolist())) == n, 'the triples are not disjoint'
    assert len(T) == target, 'got %d triples, expected %d' % (len(T), target)
    assert (V[T[:, 0]] + V[T[:, 1]] + V[T[:, 2]] == 0).all(), 'a triple does not sum to zero'
    for a in range(3):
        for b in range(a + 1, 3):
            assert (G[T[:, a], T[:, b]] == -NORM // 2).all(), 'a triple is not at cosine -1/2'
    return T


def main(write):
    print(__doc__)
    P, _G = order3()
    print('pi: cycle shape %s' % (sorted(len(o) for o in cycles_of(P)),))
    print('K_12:')
    V12, T12, A = kappa12(P)
    print('   %d minimal vectors, %d zero-sum triples from the orbits of pi'
          % (len(V12), len(T12)))
    print('K_13:')
    V13, G13, u = kappa13(V12, A)
    T13 = triple_partition(V13, G13, TRIPLES[13])
    print('   %d minimal vectors, %d zero-sum triples' % (len(V13), len(T13)))

    B12, n12 = ortho_basis(P)
    X12, c12, m12 = reduce_coords(V12, B12, n12)
    B13, n13 = ortho_basis(P, u)
    X13, c13, m13 = reduce_coords(V13, B13, n13)

    print()
    print('   k   |W|   triples   coords   norm   metric diag(c)')
    out = {12: (V12, T12, X12, c12, m12), 13: (V13, T13, X13, c13, m13)}
    for k in (12, 13):
        V, T, X, c, m = out[k]
        print('  %2d  %5d   %7d   %6d %6d   %s' % (k, len(V), len(T), X.shape[1], m,
                                                   ','.join(map(str, c.tolist()))))

    print()
    print('   k   det at minimal norm 4   expected   (a count does not identify a lattice)')
    for k in (12, 13):
        V = out[k][0]
        d = det_gram(hnf_basis(V)) / F(8) ** k
        assert d == DET4[k], 'k = %d: determinant %s, expected %d' % (k, d, DET4[k])
        print('  %2d  %21s   %8d' % (k, d, DET4[k]))

    ok = True
    for k in (12, 13):
        V, T, X, c, m = out[k]
        wf = os.path.join(DATA, 'kap%d_W.npy' % k)
        cf = os.path.join(DATA, 'kap%d_c.npy' % k)
        pf = os.path.join(DATA, 'kap%d_parts.npz' % k)
        if write:
            np.save(wf, X.astype(np.int64))
            np.save(cf, c.astype(np.int64))
            np.savez(pf, triples=T.astype(np.int64),
                     pairs=np.zeros((0, 2), dtype=np.int64))
            print('   wrote data/kap%d_W.npy (%d x %d), _c.npy and _parts.npz (%d triples)'
                  % (k, len(X), X.shape[1], len(T)))
        elif all(os.path.exists(q) for q in (wf, cf, pf)):
            S, sc, Zp = np.load(wf), np.load(cf), np.load(pf)
            same = (len(S) == len(X) and S.shape[1] == X.shape[1] and (sc == c).all()
                    and len(Zp['triples']) == len(T))
            print('   shipped kap%d: %d x %d, %d triples  %s'
                  % (k, len(S), S.shape[1], len(Zp['triples']), 'OK' if same else 'DISAGREES'))
            ok = ok and same
        else:
            print('   kap%d is not shipped yet; run with --write' % k)
            ok = False
    return ok


if __name__ == '__main__':
    sys.exit(0 if main('--write' in sys.argv) else 1)
