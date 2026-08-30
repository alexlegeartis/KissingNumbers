#!/usr/bin/env python3
"""Perfect zero-sum triple partitions of Lambda_k for the ODD k, by a complete plateau walk.

    python oddparts.py                re-derive k = 17, 19, 21, 23 and report
    python oddparts.py 17 19          just those
    python oddparts.py --write        write wherever the result beats what is shipped
    python oddparts.py 23 --seeds=12  the default; k = 23 closes on the fourth seed, the
                                      others on the first

WHAT A TRIPLE IS WORTH.  Dimension 72 + k partitions its cap direction set W into parts
internally at cosine <= -1/2; a part of three is a ZERO-SUM TRIPLE and is worth 4 per class
line where a pair is worth 2 (KNOWLEDGE.md section 50).  A class here is about 2100 lines, so
a triple the partition fails to find costs some 8400 spheres.

WHY THE SHIPPED PARTITIONS WERE SHORT.  scripts/lamtriples.py packs greedily with restarts.
It never revisits a placed triple, so it stops at a MAXIMAL packing rather than a maximum one:
1773 triples of a possible 1782 at k = 17, and 30990 of 31050 at k = 23.  Nine and sixty
triples do not sound like much.  They are 29 120 and 150 328 spheres.  The even k were closed
by an order-3 isometry (scripts/omega_parts.py), which cannot exist in odd dimension -- its
eigenvalues are omega and omega-bar and come in conjugate pairs -- so these four had been left
as the standing open item.

THE MOVE, AND WHY IT IS THE ONLY ONE.  The hypergraph of zero-sum triples is LINEAR: two
vectors a and b lie in at most one triple, since the third is forced to be -(a+b).  Take a
triangle {u,v,w} with u uncovered.  If v and w belonged to a common packed triple t, then t's
third vector would be -(v+w) = u, which is uncovered -- so they never do, and evicting the
triples that meet {u,v,w} frees either 3 vectors (one eviction) or 6 (two).  One eviction
leaves the uncovered count UNCHANGED, and it happens exactly when one of v, w is itself
uncovered.  So

    every neutral or improving move is a pair of uncovered vectors at inner product -|x|^2/2,

and with |U| in the tens the whole move set is a few hundred pairs, enumerated in one pass.
That is the whole algorithm: enumerate the pairs, take a move that covers three uncovered
vectors if one exists, else a neutral one at random, and when neither exists evict two triples
for one and carry on from there, keeping the best packing seen.  A sampler that draws a
random triangle through a random uncovered vector instead is looking for a needle: at k = 23
each vector lies in 1232 triples and all but a handful of them lead nowhere.

WHAT IT REACHES.  All four close: 1782, 3556, 9240 and 31050 triples, perfect partitions of
every one of the four shells.

THE CAP, WHICH IS THE ONE PARAMETER THAT MATTERS.  The pair scan is quadratic in |U|, so below
some threshold the walk enumerates and above it it samples.  At 150 -- the first value tried --
k = 23 stalled at 31047 or 31048 through cold restarts, iterated local search with kicks of up
to forty triples, and a tabu on recently placed triples -- four strategies, one answer, which
reads as a fact about the lattice.  It was the threshold: at 800 k = 23 closes in 55 seconds.
Sampling once the packing is nearly full is exactly where sampling is worst, and 150 was
switching to it far too early.

THE WARM START.  Lambda_(k-1)'s minimal vectors sit inside Lambda_k's as the level-0 set of
the functional orthogonal to their span, and for even k-1 that partition is already perfect --
scripts/omega_parts.py made it with an order-3 isometry.  Starting the walk there rather than
empty is a different basin: it must then break level-0 triples to serve the vectors at the
other levels, which is the trade the counting says a perfect partition makes.  Both starts are
run at every seed and the better is kept.

THE TRAP THIS AVOIDS.  W is closed under negation and a triple's sign pattern is unique up to
a global sign (a+b+c = 0 and a-b+c = 0 would give b = 0), so triples come in antipodal pairs
and it is tempting to halve the problem by packing LINES instead -- which is also exactly what
an order-3 isometry would deliver, since the orbit of -v is the negative of the orbit of v.
It is not available here.  Summing [a] + [b] + [c] = 0 in Lambda / 2 Lambda over the parts of
an antipodally symmetric perfect partition gives

    sum over the |W|/2 minimal LINES of [v]  =  0   in  Lambda / 2 Lambda,

and that fails at k = 11, 13, 15, 17, 19 and 23 -- the `sym` column below, computed by an
exact lattice membership test.  On the vectors the same sum is doubled and the condition is
vacuous.  A line-level walk duly stalls one triangle short at k = 17, at 890 of 891, and the
reason is arithmetic rather than search.  The k = 17, 19 and 23 partitions below are therefore
NOT antipodally symmetric and cannot be; k = 21 passes the test, and there it is the
eigenvalue argument alone that rules the group construction out.

    python oddparts.py
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
KS = (17, 19, 21, 23)
# fixed multipliers below 2^50: coordinates of Leech minimal vectors at squared norm 32 are
# at most 4 in absolute value over 24 coordinates, so the dot product stays under 2^58 and
# cannot wrap.  Uniqueness of the result is asserted, never assumed.
_R = np.random.default_rng(20260830).integers(1, 1 << 50, 32).astype(np.int64)


def _hash(A):
    return (A.astype(np.int64) * _R[:A.shape[1]]).sum(1)


class Lookup(object):
    """index of a row of V, or -1, for many rows at once"""

    def __init__(self, V):
        self.V = np.ascontiguousarray(V)
        h = _hash(self.V)
        o = np.argsort(h, kind='stable')
        self.h, self.o = h[o], o
        assert len(np.unique(self.h)) == len(h), 'hash collision'

    def find(self, Q):
        hq = _hash(Q)
        p = np.searchsorted(self.h, hq)
        p = np.where(p >= len(self.h), 0, p)
        idx = self.o[p].copy()
        idx[(self.h[p] != hq) | (self.V[idx] != Q).any(1)] = -1
        return idx


class Triples(object):
    """for every vector, the zero-sum triples through it, as a CSR

    The blocks come out ordered by source vector, so no global sort is needed: at k = 23 the
    table is 115 million rows and sorting it would double a 900 MB peak.
    """

    def __init__(self, V, block=64):
        self.V = V
        # float32 is EXACT here: |coordinate| <= 4 over 24 coordinates gives |<a,b>| <= 384,
        # far below 2^24, and numpy's int64 matmul is not BLAS-backed and twenty times slower
        self.Vf = np.ascontiguousarray(V.astype(np.float32))
        self.lk = Lookup(V)
        nm = int(V[0] @ V[0])
        assert (np.einsum('ij,ij->i', V, V) == nm).all(), 'unequal norms'
        n = len(V)
        cnt = np.zeros(n, np.int64)
        pieces = []
        for s0 in range(0, n, block):
            b = np.arange(s0, min(s0 + block, n))
            IP = self.Vf[b] @ self.Vf.T
            r, c = np.nonzero(IP == -0.5 * nm)
            m = self.lk.find(-(V[b[r]] + V[c]))
            ok = (m >= 0) & (c < m)
            r, c, m = r[ok], c[ok], m[ok]
            pieces.append(np.stack([c, m], 1).astype(np.int32))
            np.add.at(cnt, b[r], 1)
        self.pair = np.concatenate(pieces, 0)
        self.off = np.concatenate([[0], np.cumsum(cnt)]).astype(np.int64)
        assert len(self.pair) == self.off[-1]

    def __getitem__(self, i):
        return self.pair[self.off[i]:self.off[i + 1]]


def warm_start(k, tri):
    """Lambda_(k-1)'s perfect partition, as triples of Lambda_k, or None

    Lambda_(k-1)'s minimal vectors sit inside Lambda_k's as the level-0 set of the functional
    orthogonal to their span, and for even k-1 that partition is already perfect (an order-3
    isometry made it).  Starting the walk there is a different basin from starting empty: the
    walk must then break level-0 triples to serve the vectors at the other levels, which is
    the trade the counting says a perfect partition makes.  At k = 23 that is the difference
    between 31047 and 31048; below k = 23 the cold start closes on its own.
    """
    wf = os.path.join(DATA, 'lam%d_W.npy' % (k - 1))
    pf = os.path.join(DATA, 'parts_%d.npz' % (k - 1))
    if not (os.path.exists(wf) and os.path.exists(pf)):
        return None
    W = np.load(wf).astype(np.int64)
    if W.shape[1] != tri.V.shape[1]:         # written in different coordinates: k = 17 only
        return None
    idx = tri.lk.find(np.ascontiguousarray(W))
    if (idx < 0).any():                      # not a subset in these coordinates
        return None
    T = np.load(pf)['triples']
    if 3 * len(T) != len(W):                 # only a perfect one is worth seeding with
        return None
    out = [tuple(sorted(int(idx[x]) for x in t)) for t in T.tolist()]
    for a, b, c in out:
        assert (tri.V[a] + tri.V[b] + tri.V[c] == 0).all(), 'seed triple does not sum to zero'
    return out


def walk(V, tri, rng, iters, cap=800, start=None):
    """the plateau walk described above; returns (uncovered, triples)"""
    n = len(V)
    Vf, lk = tri.Vf, tri.lk
    nm = float(V[0] @ V[0])
    part = np.full(n, -1, np.int32)
    slots, free = [], []
    if start:
        slots = [tuple(t) for t in start]
        for _i, _t in enumerate(slots):
            part[list(_t)] = _i
    unc = [i for i in range(n) if part[i] < 0]
    best = (len(unc), [t for t in slots if t is not None])
    # the uncovered count is maintained incrementally: a move always places three vectors
    # that are uncovered at the moment it places them, and an eviction always frees three
    state = [len(unc)]

    def evict(t):
        for y in slots[t]:
            part[y] = -1
            unc.append(y)
        slots[t] = None
        free.append(t)
        state[0] += 3

    def place(a, b, c):
        sl = free.pop() if free else len(slots)
        if sl == len(slots):
            slots.append(None)
        slots[sl] = (a, b, c)
        for y in (a, b, c):
            assert part[y] < 0
            part[y] = sl
        state[0] -= 3

    for _ in range(iters):
        if len(unc) > 4 * cap + 4 * state[0]:
            unc = [i for i in unc if part[i] < 0]
        nu = state[0]
        if nu < best[0]:
            best = (nu, [t for t in slots if t is not None])
            if nu == 0:
                break
        moved = False
        if nu <= cap:
            unc = [i for i in unc if part[i] < 0]
            U = np.array(sorted(set(unc)), np.int64)
            IP = Vf[U] @ Vf[U].T
            i, j = np.nonzero(IP == -0.5 * nm)
            m = i < j
            i, j = i[m], j[m]
            if len(i):
                third = lk.find(-(V[U[i]] + V[U[j]]))
                ok = third >= 0
                i, j, third = i[ok], j[ok], third[ok]
                if len(i):
                    p3 = part[third]
                    good = np.nonzero(p3 < 0)[0]
                    pick = good if len(good) else np.arange(len(i))
                    r = int(pick[rng.integers(len(pick))])
                    if p3[r] >= 0:
                        evict(int(p3[r]))
                    place(int(U[i[r]]), int(U[j[r]]), int(third[r]))
                    moved = True
        if not moved:                       # no plateau move: two triples out, one in
            u = unc[int(rng.integers(len(unc)))]
            while part[u] >= 0:
                u = unc[int(rng.integers(len(unc)))]
            d = tri[u]
            pv, pw = part[d[:, 0]], part[d[:, 1]]
            c0 = np.nonzero((pv < 0) & (pw < 0))[0]
            if not len(c0):
                c0 = np.nonzero((pv < 0) ^ (pw < 0))[0]
            if not len(c0):
                c0 = np.arange(len(d))
            r = int(c0[rng.integers(len(c0))])
            v, w = int(d[r, 0]), int(d[r, 1])
            for x in (v, w):
                if part[x] >= 0:
                    evict(int(part[x]))
            place(u, v, w)
    return best


def hnf(rows):
    """row Hermite normal form over Z, exact in Python ints"""
    M = [list(map(int, r)) for r in rows]
    m, n = len(M), len(M[0])
    r = 0
    for c in range(n):
        piv = None
        for i in range(r, m):
            if M[i][c] and (piv is None or abs(M[i][c]) < abs(M[piv][c])):
                piv = i
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        again = True
        while again:
            again = False
            for i in range(r + 1, m):
                if M[i][c]:
                    q = M[i][c] // M[r][c]
                    M[i] = [a - q * b for a, b in zip(M[i], M[r])]
                    if M[i][c]:
                        M[r], M[i] = M[i], M[r]
                        again = True
        r += 1
    return [row for row in M[:r] if any(row)]


def in_lattice(v, basis):
    v = list(map(int, v))
    for row in basis:
        c = next(i for i, x in enumerate(row) if x)
        if v[c] % row[c]:
            return False
        q = v[c] // row[c]
        v = [a - q * b for a, b in zip(v, row)]
    return not any(v)


def symmetric_possible(W):
    """can a perfect partition be antipodally symmetric?  sum of the line classes mod 2"""
    first = np.argmax(W != 0, axis=1)
    C = W * np.sign(W[np.arange(len(W)), first])[:, None]
    keys, L = set(), []
    for row in C:
        t = tuple(int(x) for x in row)
        if t not in keys:
            keys.add(t)
            L.append(t)
    assert 2 * len(L) == len(W), 'the shell is not antipodally closed'
    S = [sum(col) for col in zip(*L)]
    return all(x % 2 == 0 for x in S) and in_lattice([x // 2 for x in S], hnf(W.tolist()))


def pair_up(V, free):
    """antipodal pairs on whatever the walk left over"""
    out, used = [], set()
    for a in free:
        if a in used:
            continue
        for b in free:
            if b != a and b not in used and 2 * int(V[a] @ V[b]) <= -int(V[a] @ V[a]):
                used.add(a)
                used.add(b)
                out.append((a, b))
                break
    return out


def main(ks, write, seeds, iters):
    print(__doc__)
    print('    k     |W|   perfect   shipped T + P    found T + P     2T+P    was   sym   note')
    print('    ' + '-' * 90)
    gained = 0
    for k in ks:
        V = np.load(os.path.join(DATA, 'lam%d_W.npy' % k)).astype(np.int64)
        pf = 'parts_%d.npz' % k
        z = np.load(os.path.join(DATA, pf))
        T0, P0 = len(z['triples']), len(z['pairs'])
        sym = symmetric_possible(V)
        t0 = time.time()
        tri = Triples(V)
        seedT = warm_start(k, tri)
        best = (len(V), [])
        for sd in range(1, seeds + 1):
            for st in ((None, seedT) if seedT else (None,)):
                b = walk(V, tri, np.random.default_rng(sd), iters, start=st)
                if b[0] < best[0]:
                    best = b
                if best[0] == 0:
                    break
            if best[0] == 0:
                break
        nu, sol = best
        cover = np.zeros(len(V), bool)
        for t in sol:
            cover[list(t)] = True
        pr = pair_up(V, np.nonzero(~cover)[0].tolist())
        for a, b, c in sol:                       # exact, every one of them
            assert (V[a] + V[b] + V[c] == 0).all(), 'not a zero-sum triple'
        assert len(set(x for t in sol for x in t)) == 3 * len(sol), 'parts overlap'
        now, was = 2 * len(sol) + len(pr), 2 * T0 + P0
        print('   %2d %7d   %7d   %5d + %-5d    %5d + %-5d  %7d %6d   %-5s %s'
              % (k, len(V), len(V) // 3, T0, P0, len(sol), len(pr), now, was,
                 'yes' if sym else 'no', 'PERFECT' if nu == 0 else '%d uncovered' % nu),
              flush=True)
        if now > was:
            gained += 1
            if write:
                np.savez(os.path.join(DATA, pf),
                         triples=np.array(sorted(tuple(sorted(t)) for t in sol), np.int64),
                         pairs=np.array(sorted(pr), np.int64).reshape(-1, 2))
                print('        wrote data/%s: %d triples, %d pairs  (%.0f s)'
                      % (pf, len(sol), len(pr), time.time() - t0))
    print()
    print('   the sym column is the mod-2 test: "no" means NO perfect partition of that shell')
    print('   is antipodally symmetric, so no order-3 isometry can produce one either.')
    if write:
        print('   %d partitions improved; run final81-95.py and scripts/verify_parts.py' % gained)
    else:
        print('   %d partitions would improve; re-run with --write' % gained)
    return 0


if __name__ == '__main__':
    _ks = [int(a) for a in sys.argv[1:] if a.isdigit()] or list(KS)
    _s = next((int(a.split('=')[1]) for a in sys.argv[1:] if a.startswith('--seeds=')), 12)
    _i = next((int(a.split('=')[1]) for a in sys.argv[1:] if a.startswith('--iters=')), 1200000)
    sys.exit(main(_ks, '--write' in sys.argv, _s, _i))
