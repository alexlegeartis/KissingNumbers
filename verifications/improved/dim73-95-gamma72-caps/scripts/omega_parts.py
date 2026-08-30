#!/usr/bin/env python3
"""Perfect zero-sum triple partitions, from an order-3 isometry instead of a packing search.

    python omega_parts.py             re-derive and report against the shipped partitions
    python omega_parts.py --write     write the partition wherever it beats what is shipped
    python omega_parts.py --write 20  just that k

WHAT THE CAP LAYER IS PAID FOR.  Dimension 72 + k partitions its cap direction set W into
parts internally at cosine <= -1/2; a part of three is a ZERO-SUM TRIPLE and is worth 4 per
class line where a pair is worth 2 (KNOWLEDGE.md section 50).  The gain is 4 S(T) + 2 (S(T+P)
- S(T)) over the largest classes, so every triple that a partition fails to find is paid for
at roughly four times the size of one class -- eight thousand spheres apiece.

WHERE THE SHIPPED PARTITIONS CAME FROM, AND WHY THEY WERE SHORT.  scripts/lamtriples.py packs
triples greedily, most-constrained-first with restarts.  That is a good heuristic and it was
never perfect: at k = 16 it found 1434 triples of a possible 1440, at k = 18 2457 of 2466.
Six and nine triples do not sound like much.  They are 20 820 and 29 092 spheres.

A GROUP DOES IT EXACTLY.  Suppose S is an isometry of the axis space with

    I + S + S^2 = 0.

Then S has eigenvalues omega and omega-bar only, so it has no fixed vector and acts FREELY;
and for every v, the orbit {v, Sv, S^2 v} sums to (I + S + S^2) v = 0, so it is a zero-sum
triple outright.  If W is S-invariant, its partition into orbits is a PERFECT triple
partition, found by a group rather than packed.  This is what makes K_12 free in
scripts/kappa_caps.py, where S is a Leech automorphism of order 3; here S is searched for
directly, inside the axis space, for cap sets that were never suspected of having one.

ONLY IN EVEN DIMENSION.  Those eigenvalues come in conjugate pairs, so a space carrying such
an S has even dimension.  k = 9, 11, 13, 15, 17, 19, 21 and 23 are ruled out by that alone and
are not searched here.  They are not open either: scripts/oddparts.py closes k = 17, 19, 21
and 23 by a local search whose move set is derived rather than chosen, and the rest were
already perfect.  There is a second, independent reason the group cannot serve those k, which
that file states: an S-invariant partition is antipodally SYMMETRIC, and a mod-2 test rules a
symmetric perfect partition out for Lambda_k at k = 11, 13, 15, 17, 19 and 23.

THE SEARCH.  S is determined by its action on a basis, and choosing S(v) = w forces
S(w) = -v-w and S(-v-w) = v, so one decision fixes a whole triple.  Every later image must be
an isometry against all of them -- <S v_i, S v_j> = <v_i, v_j> -- which is what collapses a
space of 82^k assignments to a few dozen nodes.

THE TRIPLE-DEGREE FILTER, which is what makes rank 20 and 22 reachable at all.  The isometry
constraint only starts biting at the SECOND level, so at k = 20 the first level branches 664
ways and the search reached zero complete assignments in 900 seconds.  A second invariant was
sitting there unused: S maps zero-sum triples to zero-sum triples, so it preserves the number
of them through a vector,

    deg(S v) = deg(v),

and those degrees are not constant -- 664, some intermediate value and 712 at k = 20, three
values in all.  Filtering candidate images by degree takes the branching from 664 to EIGHT,
and k = 20 falls in 607 complete assignments.  The basis is ordered rarest-degree-first, so
the narrowest levels are decided first.  Computing the degrees is one blocked pass over the
Gram, which is the same pass that finds the triples.

CONSISTENCY ON A BASIS IS NOT ENOUGH, and this file was wrong about that once.  A basis
assignment fixes S on the span but says nothing about whether S maps W into W: the first
assignment reached is set-preserving at k = 18 and is NOT at k = 12 or 16, where the search
has to continue to the 18th and 77th.  So every complete assignment is tested by applying S
to the whole of W, in integer arithmetic -- coordinates against the basis are cleared of
denominators once and the image is three matrix products, not |W| eliminations.

    k     |W|      packed       orbits    worth
    12    756      252          252       the same partition, reached two ways
    16    4320     1434         1440      +20 820
    18    7398     2457         2466      +29 092
    20    17400    5782 + 23    5800      +47 592

scripts/verify_parts.py re-checks whatever is written here from scratch, without this file:
every triple summing to zero, every inner product exactly -|x|^2/2, the parts disjoint.
"""
import os
import sys
import time
from fractions import Fraction as F
from math import gcd

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

# (cap file, metric file or None, parts file) per source kind, in the order tried
SOURCES = (('kap%d_W.npy', 'kap%d_c.npy', 'kap%d_parts.npz'),
           ('lam%d_W.npy', None, 'parts_%d.npz'),
           ('rec_%d_W.npy', 'rec_%d_c.npy', 'rec_%d_parts.npz'))


def capset(k):
    """(X, c, parts file, cap file) for the largest shipped cap set at this k"""
    best = None
    for wf, cf, pf in SOURCES:
        p = os.path.join(DATA, wf % k)
        if not os.path.exists(p) or not os.path.exists(os.path.join(DATA, pf % k)):
            continue
        X = np.load(p).astype(np.int64)
        c = (np.load(os.path.join(DATA, cf % k)).astype(np.int64) if cf
             else np.ones(X.shape[1], np.int64))
        keep = (np.abs(X).sum(0) > 0) | (c != 1)
        cand = (len(X), X[:, keep], c[keep], pf % k, wf % k)
        if best is None or cand[0] > best[0]:
            best = cand
    return best[1:] if best else (None, None, None, None)


def inv_int(G):
    """(adj, den) with adj @ G = den * I, exactly"""
    n = len(G)
    M = [[F(int(G[i][j])) for j in range(n)] + [F(int(i == j)) for j in range(n)]
         for i in range(n)]
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c]), None)
        if p is None:
            return None, None
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    den = 1
    for i in range(n):
        for j in range(n):
            d = M[i][n + j].denominator
            den = den * d // gcd(den, d)
    return np.array([[int(M[i][n + j] * den) for j in range(n)] for i in range(n)],
                    dtype=object), den


def omega(k, X, c, budget):
    """(T, note): the orbit partition of an S with I + S + S^2 = 0, or None"""
    n, d = X.shape
    rank = int(np.linalg.matrix_rank(X.astype(np.float64)))
    if rank % 2:
        return None, 'odd dimension: no such isometry exists'
    m = int((X[0] * c) @ X[0])
    Y = (X * c).astype(np.int64)
    idx = {tuple(r): i for i, r in enumerate(X.tolist())}
    t0 = time.time()

    # every vector's triple-degree, in one blocked pass; float32 holds these Gram entries
    # exactly, which is asserted rather than assumed
    part_all = [[] for _ in range(n)]
    Xf = X.astype(np.float32)
    Yf = Y.astype(np.float32)
    assert int(np.abs(X).max()) * int(np.abs(Y).max()) * d < 2 ** 24, 'float32 exactness'
    BS = max(1, int(3e7 // n))
    for s0 in range(0, n, BS):
        Bm = (Xf[s0:s0 + BS] @ Yf.T).astype(np.int64)
        rr, cc = np.nonzero(Bm == -m // 2)
        for a, b in zip(rr, cc):
            a = int(a) + s0
            b = int(b)
            if idx.get(tuple((-(X[a] + X[b])).tolist())) is not None:
                part_all[a].append(b)
    deg = np.array([len(p) for p in part_all], dtype=np.int64)
    if deg.min() == 0:
        return None, ('%d vectors lie in no zero-sum triple, so no such isometry exists'
                      % int((deg == 0).sum()))

    # rarest degree first: the narrowest levels are decided first
    rare = {int(v): int(c) for v, c in zip(*np.unique(deg, return_counts=True))}
    basis = []
    for i in sorted(range(n), key=lambda t: (rare[int(deg[t])], int(deg[t]), t)):
        if np.linalg.matrix_rank(X[basis + [i]].astype(np.float64)) == len(basis) + 1:
            basis.append(i)
        if len(basis) == rank:
            break
    B = X[basis].astype(np.int64)
    adj, den = inv_int(B @ B.T)
    if adj is None:
        return None, 'the basis is degenerate'
    CO = (X.astype(object) @ B.T.astype(object)) @ adj                  # n x rank, x den
    assert (CO @ B.astype(object) == den * X.astype(object)).all(), 'coordinate solve wrong'

    rows = {}

    def gram_row(i):
        if i not in rows:
            rows[i] = X @ Y[i]
        return rows[i]

    # candidate images: a triple partner OF THE SAME DEGREE
    part = {i: [j for j in part_all[i] if deg[j] == deg[i]] for i in basis}
    src, img = [], []
    tried = [0]

    def push(u, z):
        add = 0
        h = idx[tuple((-(X[u] + X[z])).tolist())]
        for p, q in ((u, z), (z, h), (h, u)):
            if p in src:
                if img[src.index(p)] != q:
                    for _ in range(add):
                        src.pop()
                        img.pop()
                    return False
                continue
            if deg[q] != deg[p]:
                for _ in range(add):
                    src.pop()
                    img.pop()
                return False
            gq, gp = gram_row(q), gram_row(p)
            if any(int(gq[img[t]]) != int(gp[src[t]]) for t in range(len(src))):
                for _ in range(add):
                    src.pop()
                    img.pop()
                return False
            src.append(p)
            img.append(q)
            add += 1
        return add

    def check():
        """S on a basis is not S on the set: apply it to all of W and look"""
        tried[0] += 1
        Im = X[[img[src.index(i)] for i in basis]].astype(object)
        Q = CO @ Im
        if not (Q % den == 0).all():
            return None
        out = [tuple(int(v) for v in r) for r in Q // den]
        if set(out) != set(map(tuple, X.tolist())):
            return None
        mp = np.array([idx[o] for o in out])
        if (mp == np.arange(n)).any():
            return None                                    # S fixes a vector: not free
        seen = np.zeros(n, bool)
        T = []
        for i in range(n):
            if seen[i]:
                continue
            a, b = int(mp[i]), int(mp[mp[i]])
            if seen[a] or seen[b]:
                return None
            seen[[i, a, b]] = True
            T.append((i, a, b))
        T = np.array(T, dtype=np.int64)
        if not seen.all() or not (X[T[:, 0]] + X[T[:, 1]] + X[T[:, 2]] == 0).all():
            return None
        return T

    res = [None]

    def bt(t):
        if time.time() - t0 > budget:
            raise TimeoutError
        if t == rank:
            T = check()
            if T is None:
                return False
            res[0] = T
            return True
        for z in part[basis[t]]:
            r = push(basis[t], z)
            if r is False:
                continue
            if bt(t + 1):
                return True
            for _ in range(r):
                src.pop()
                img.pop()
        return False

    try:
        bt(0)
    except TimeoutError:
        return None, 'no isometry within the budget (%d assignments tried)' % tried[0]
    if res[0] is None:
        return None, 'none exists for this cap set (%d assignments tried)' % tried[0]
    return res[0], '%d assignments tried' % tried[0]


def main(write, ks, budget):
    print(__doc__)
    print('   k     |W|   dim   shipped T + P    orbits    2T+P  was   note')
    print('   ' + '-' * 92)
    gained = 0
    for k in ks:
        X, c, pf, wf = capset(k)
        if X is None:
            continue
        n = len(X)
        rank = int(np.linalg.matrix_rank(X.astype(np.float64)))
        z = np.load(os.path.join(DATA, pf))
        T0, P0 = len(z['triples']), len(z['pairs'])
        was = 2 * T0 + P0
        T, note = omega(k, X, c, budget)
        if T is None:
            print('   %2d %7d %5d   %5d + %-5d    %-9s %5s %5d   %s'
                  % (k, n, rank, T0, P0, '-', '-', was, note))
            continue
        now = 2 * len(T)
        print('   %2d %7d %5d   %5d + %-5d    %-9d %5d %5d   %s'
              % (k, n, rank, T0, P0, len(T), now, was, note))
        if now > was:
            gained += 1
            if write:
                np.savez(os.path.join(DATA, pf), triples=T.astype(np.int64),
                         pairs=np.zeros((0, 2), dtype=np.int64))
                print('        wrote data/%s: %d triples, no pairs, no singletons'
                      % (pf, len(T)))
    print()
    if write:
        print('   %d partitions improved; run final81-95.py and scripts/verify_parts.py'
              % gained)
    else:
        print('   %d partitions would improve; re-run with --write' % gained)
    return 0


if __name__ == '__main__':
    # k = 22 needs the longer budget: its triple-degree is CONSTANT at 1492, so the filter
    # buys nothing there and the first level stays 1492 wide -- the isometry condition alone
    # still settles it, in twelve complete assignments, but the pass that finds the triples
    # is 1.2e9 pairs before the search starts.
    _ks = [int(a) for a in sys.argv[1:] if a.isdigit()] or list(range(9, 24))
    _b = next((float(a.split('=')[1]) for a in sys.argv[1:] if a.startswith('--budget=')),
              3600.0)
    sys.exit(main('--write' in sys.argv, _ks, _b))
