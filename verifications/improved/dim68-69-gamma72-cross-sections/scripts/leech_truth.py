#!/usr/bin/env python3
"""Exact ground truth for the k-point LP's Leech validation cases.

The LP in ../../../common/kpoint_lp.py certifies a LOWER bound on the number of minimal
vectors orthogonal to a
k-tuple of prescribed Gram.  For the Leech lattice that number can be counted outright, so the
validation targets need not be quoted from the literature: this script builds all 196 560
minimal vectors, finds an actual tuple realising each Gram, and counts.

Coordinates are the standard integral model, Lambda_24 scaled so that a minimal vector has
sum of squares 32; the lattice inner product is then <x,y> = x.y / 8, giving minimum norm 4.
The three shapes are (4^2, 0^22), (2^8, 0^16) on a Golay octad, and (-3, 1^23) up to sign.
"""
import itertools
import numpy as np

# ---- the extended binary Golay code -----------------------------------------------------
# The cyclic [23,12,7] Golay code has generator polynomial
#     g(x) = x^11 + x^10 + x^6 + x^5 + x^4 + x^2 + 1,
# a divisor of x^23 - 1 over GF(2).  Its 12 cyclic shifts generate the code; appending an
# overall parity bit gives the extended [24,12,8] code, whose 759 weight-8 words are the
# octads.  The weight distribution 1/759/2576/759/1 is asserted below, so a wrong polynomial
# cannot pass silently.
GPOLY = (0, 2, 4, 5, 6, 10, 11)


def golay():
    G = np.zeros((12, 24), dtype=np.int8)
    for i in range(12):
        for e in GPOLY:
            G[i, i + e] = 1
        G[i, 23] = G[i, :23].sum() % 2      # overall parity bit
    words = np.zeros((4096, 24), dtype=np.int8)
    for n in range(4096):
        w = np.zeros(24, dtype=np.int8)
        for i in range(12):
            if n >> i & 1:
                w ^= G[i]
        words[n] = w
    wt = words.sum(axis=1)
    dist = {int(k): int((wt == k).sum()) for k in np.unique(wt)}
    assert dist == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}, dist
    return words


def octads(words):
    return words[words.sum(axis=1) == 8]


def in_leech(V, codes):
    """the standard membership test, applied row-wise to a candidate array

    x lies in Lambda_24 (this scaling) iff all coordinates are congruent mod 2, say to m, the
    set {i : x_i = m+2 mod 4} is a Golay codeword, and sum x_i = 4m mod 8."""
    V = V.astype(np.int64)
    m = V[:, 0] & 1
    par = ((V & 1) == m[:, None]).all(axis=1)
    S = ((V - (m[:, None] + 2)) % 4) == 0
    bits = (S.astype(np.int64) * (1 << np.arange(24, dtype=np.int64))).sum(axis=1)
    isc = np.isin(bits, codes)
    su = (V.sum(axis=1) - 4 * m) % 8 == 0
    return par & isc & su


def minimal_vectors():
    """all 196 560 minimal vectors, generated as candidates and filtered by in_leech

    Nothing here relies on getting the sign conventions of the three shapes right: candidates
    are produced generously and the membership test decides.  The shape counts are asserted."""
    words = golay()
    codes = np.sort((words.astype(np.int64) * (1 << np.arange(24, dtype=np.int64))).sum(axis=1))
    cand = []
    # shape (4^2, 0^22)
    for i, j in itertools.combinations(range(24), 2):
        for si in (4, -4):
            for sj in (4, -4):
                v = np.zeros(24, dtype=np.int16)
                v[i], v[j] = si, sj
                cand.append(v)
    n1 = len(cand)
    # shape (2^8, 0^16) on the octads, every sign pattern
    for w in words[words.sum(axis=1) == 8]:
        sup = np.nonzero(w)[0]
        for mask in range(256):
            v = np.zeros(24, dtype=np.int16)
            for t, q in enumerate(sup):
                v[q] = -2 if mask >> t & 1 else 2
            cand.append(v)
    n2 = len(cand) - n1
    # shape (+-3, +-1^23): the +-1 pattern is a Golay word, both signs of the 3 offered
    for w in words:
        base = np.where(w == 1, -1, 1).astype(np.int16)
        for q in range(24):
            for sg in (3, -3):
                v = base.copy()
                v[q] = sg
                cand.append(v)
    n3 = len(cand) - n1 - n2
    V = np.array(cand, dtype=np.int16)
    ok = in_leech(V, codes) & ((V.astype(np.int32) ** 2).sum(axis=1) == 32)
    V = V[ok]
    got = (int(ok[:n1].sum()), int(ok[n1:n1 + n2].sum()), int(ok[n1 + n2:].sum()))
    assert got == (1104, 97152, 98304), got
    assert len(V) == 196560, len(V)
    return V


def find_tuple(V, G, mu=4, tries=200, seed=1):
    """a k-tuple of minimal vectors with the prescribed Gram, or None

    G is in lattice units, in which a minimal vector has norm mu.  These coordinates have raw
    square length 32 for a minimal vector, so the raw dot product of a lattice inner product g
    is g * 32/mu; the Gram of the tuple found is asserted against G before returning."""
    K = len(G)
    scale = 32 // mu
    rng = np.random.default_rng(seed)
    Vi = V.astype(np.int32)
    for _ in range(tries):
        pick = [int(rng.integers(len(V)))]
        ok = True
        for l in range(1, K):
            want = np.array([G[l][t] * scale for t in range(l)], dtype=np.int32)
            ips = Vi[pick] @ Vi.T
            cand = np.nonzero((ips == want[:, None]).all(axis=0))[0]
            if not len(cand):
                ok = False
                break
            pick.append(int(cand[rng.integers(len(cand))]))
        if ok:
            U = Vi[pick]
            gram = (U @ U.T) // scale
            assert (gram == np.array(G)).all(), (gram, G)
            return U
    return None


if __name__ == '__main__':
    V = minimal_vectors()
    print("minimal vectors of Lambda_24: %d" % len(V))
    assert len(V) == 196560
    Vi = V.astype(np.int32)
    n0 = (Vi @ Vi[0]) == 0
    print("marginals about one vector: <v,u> = 0: %d, 1: %d, 2: %d, 4: %d" % (
        n0.sum(), ((Vi @ Vi[0]) == 8).sum(), ((Vi @ Vi[0]) == 16).sum(),
        ((Vi @ Vi[0]) == 32).sum()))
    CASES = [([[4, 0], [0, 4]], 'k=2 orth'), ([[4, 1], [1, 4]], 'k=2 cos 1/4'),
             ([[4, 2], [2, 4]], 'k=2 60deg'),
             ([[4, 0, 0], [0, 4, 0], [0, 0, 4]], 'k=3 orth'),
             ([[4, 2, 2], [2, 4, 2], [2, 2, 4]], 'k=3 A_3'),
             ([[4, 0, 0, 0], [0, 4, 0, 0], [0, 0, 4, 0], [0, 0, 0, 4]], 'k=4 orth'),
             ([[4, 2, 2, 2], [2, 4, 2, 2], [2, 2, 4, 2], [2, 2, 2, 4]], 'k=4 A_4'),
             ([[4, -2, -2, -2], [-2, 4, 0, 0], [-2, 0, 4, 0], [-2, 0, 0, 4]], 'k=4 D_4')]
    print()
    print("exact count of minimal vectors orthogonal to a tuple of the given Gram:")
    for G, name in CASES:
        counts = set()
        for seed in range(1, 13):
            U = find_tuple(V, G, seed=seed)
            if U is None:
                counts.add(None)
                continue
            counts.add(int(((Vi @ U.T) == 0).all(axis=1).sum()))
        c = sorted(counts, key=lambda x: (x is None, x))
        print("   %-10s %s%s" % (name, c[0] if len(c) == 1 else c,
                                 '' if len(c) == 1 else '   <- NOT a single orbit'))

    # The sampler above needs luck: at k = 3 the second orbit turns up in about 2.8% of
    # random orthogonal triples, so twelve seeds miss it roughly seven times in ten and
    # the loop then prints a single number and looks conclusive.  Exhibit the two triples
    # instead, which settles it without sampling.  Both are pairwise orthogonal with every
    # norm mu, so both have Gram mu*I_3 -- and their cross-sections differ.  The count is
    # therefore not a function of the Gram matrix, and any single 'truth' for this row is
    # one embedding among several.
    def vec(pairs):
        x = np.zeros(24, dtype=np.int64)
        for idx, s in pairs:
            x[idx] = 4 * s
        return x

    TRIPLES = [
        ('disjoint pairs', [[(0, 1), (1, 1)], [(2, 1), (3, 1)], [(4, 1), (5, 1)]], 19530),
        ('one shared pair', [[(0, 1), (1, 1)], [(0, 1), (1, -1)], [(2, 1), (3, 1)]], 19962),
    ]
    print()
    print('two orthogonal triples of the SAME Gram matrix, exhibited:')
    seen = []
    for name, spec, want in TRIPLES:
        U = np.array([vec(p) for p in spec], dtype=np.int64)
        assert ((U * U).sum(1) == 32).all(), 'not minimal vectors'
        for r in U:
            assert (Vi == r).all(axis=1).any(), 'not in the shell: %s' % r
        G = (U @ U.T) // 8
        assert (G == 4 * np.eye(3, dtype=np.int64)).all(), G
        got = int(((Vi @ U.T) == 0).all(axis=1).sum())
        assert got == want, (name, got, want)
        seen.append(got)
        print('   %-16s Gram = mu*I_3, cross-section %d' % (name, got))
    assert seen == [19530, 19962], seen
    print('   the two ends of the certified interval [19530, 19962], both attained,')
    print('   so the cross-section is NOT determined by the Gram matrix here.')
