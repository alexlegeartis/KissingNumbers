#!/usr/bin/env python3
"""Dimensions 62 and 63 from level 0 of Edel-Rains-Sloane, with nothing added.

    tau(62) >=  67 108 864 = 2^26      (previously 52 417 932)
    tau(63) >= 134 217 728 = 2^27      (previously 52 418 564)

THE CONSTRUCTION.  Let C be a binary linear code of length n_0, dimension k, minimum distance
d.  Map each codeword c to the sign vector

    x(c) = ( (-1)^c_1, ..., (-1)^c_n0 ) / sqrt(n_0)   in R^n0,

a unit vector.  Two codewords at Hamming distance t give

    <x(c), x(c2)> = (n_0 - 2t) / n_0 ,

so every pair is at 60 degrees or more exactly when t >= n_0/4 for every attainable t, that is
when

    d >= ceil(n_0 / 4).

The 2^k sign vectors are then a kissing configuration of R^n0, and of R^n for every n >= n_0
by inclusion, so

    tau(n) >= max over n_0 <= n of A(n_0, ceil(n_0/4)).

This is level 0 of the Edel-Rains-Sloane construction with a single support: no chain, no
support codes, no glue.  It is one line, and it is why dimensions 62 and 63 are here rather
than in ../dim49-63-p48-caps/, whose construction cannot reach these values -- not with a
better class family, not with any class family; see that package's scripts/ceiling.py.

WHAT IS VERIFIED HERE.  The geometry is checked in exact rational arithmetic for every
attainable Hamming distance, with a negative control one below the threshold.  The whole
mechanism is then run END TO END on three codes that ARE exhibited in coordinates and whose
minimum distance is verified by enumerating every codeword -- including Reed-Muller
RM(2,5) = [32,16,8], where d = ceil(32/4) = 8 holds with equality, so the boundary of the
condition is exercised.  Every PAIR is covered too: for a linear code the differences are
again codewords, so the largest inner product is exactly (n_0 - 2d)/n_0 and the minimum
distance settles all pairs at once.  That identity is itself confirmed against brute force
over all 4096^2 pairs of the Golay code, where the quadratic form still fits in memory.

WHAT IS TAKEN FROM THE LITERATURE.  Only the two code parameters

    [62,26,16]  and  [63,27,16]      (Grassl, codetables.de, read 2026-08-23)

Both are shortenings of one [64,28,16], itself built as: cyclic [73,36,16] with generator
polynomial x^37 + x^36 + x^34 + x^33 + x^32 + x^27 + x^25 + x^24 + x^22 + x^21 + x^19 + x^18
+ x^15 + x^11 + x^10 + x^8 + x^7 + x^5 + x^3 + 1, punctured to [72,36,15], Construction B2 to
[63,28,15], extended to [64,28,16]; shortening once gives [63,27,16] and twice [62,26,16].
Shortening an [n,k,d] code gives [n-1,k-1,>=d], so d = 16 survives both times, and
ceil(62/4) = ceil(63/4) = 16.  Exhibiting that [64,28,16] in coordinates is the one thing this
package does not do, and the chain is recorded so that it can be.  That is the same standing
as A(38,10) and A(39,10) in ../dim39-ers-constant-weight/, which are cited and not rebuilt.
"""
import itertools
from fractions import Fraction as F

import numpy as np

FAIL = []


def check(label, ok, detail=''):
    print("  %-4s %s%s" % ('PASS' if ok else 'FAIL', label, ('   ' + detail) if detail else ''))
    if not ok:
        FAIL.append(label)


def geometry(n0, d):
    """the worst inner product over every attainable Hamming distance, exactly"""
    worst = None
    firstbad = None
    for t in range(d, n0 + 1):
        ip = F(n0 - 2 * t, n0)
        if worst is None or ip > worst:
            worst = ip
        if ip > F(1, 2) and firstbad is None:
            firstbad = t
    return firstbad is None, firstbad, worst


def rref(G):
    G = G.copy() % 2
    r = 0
    for c in range(G.shape[1]):
        piv = next((i for i in range(r, G.shape[0]) if G[i, c]), None)
        if piv is None:
            continue
        G[[r, piv]] = G[[piv, r]]
        for i in range(G.shape[0]):
            if i != r and G[i, c]:
                G[i] ^= G[r]
        r += 1
        if r == G.shape[0]:
            break
    return G[:r]


def rm(r, m):
    """Reed-Muller RM(r,m): all monomials of degree <= r evaluated on GF(2)^m"""
    pts = np.array(list(itertools.product((0, 1), repeat=m)), dtype=np.int8)
    rows = []
    for deg in range(r + 1):
        for S in itertools.combinations(range(m), deg):
            v = np.ones(len(pts), dtype=np.int8)
            for i in S:
                v = v & pts[:, i]
            rows.append(v)
    return rref(np.array(rows, dtype=np.int8))


def golay24():
    """the extended binary Golay code [24,12,8]: cyclic [23,12,7] plus an overall parity bit"""
    GPOLY = (0, 2, 4, 5, 6, 10, 11)
    G = np.zeros((12, 24), dtype=np.int8)
    for i in range(12):
        for e in GPOLY:
            G[i, i + e] = 1
        G[i, 23] = G[i, :23].sum() % 2
    return G


def codewords(G):
    k = G.shape[0]
    assert k <= 20, 'refusing to enumerate more than a million codewords'
    idx = np.arange(1 << k, dtype=np.uint32)
    bits = ((idx[:, None] >> np.arange(k, dtype=np.uint32)) & 1).astype(np.int8)
    return (bits @ G) % 2


def min_distance(C):
    w = C.sum(axis=1)
    return int(w[w > 0].min())


def worst_pair_linear(n0, d):
    """max over distinct pairs of n_0 * <x(c), x(c2)>, for a LINEAR code

    The differences of codewords are again codewords, so the largest inner product comes from
    the smallest nonzero weight: n_0 - 2*d.  This is exact and costs nothing, which matters --
    the quadratic form for RM(2,5) would be 65536^2 entries."""
    return n0 - 2 * d


def worst_pair_exhaustive(C, n0):
    """the same quantity by brute force over all pairs, for small codes, to check the identity"""
    S = 1 - 2 * C.astype(np.int32)
    best = -n0
    for i in range(0, len(S), 512):
        blk = S[i:i + 512] @ S.T
        for r, row in enumerate(blk):
            row[i + r] = -n0
        best = max(best, int(blk.max()))
    return best



# The neighbourhood of Grassl's table, read 2026-08-23: best known minimum distance (lower
# bound) for a binary linear code of length n and dimension k.  Recorded so the two cited
# parameters can be checked in context, and so it is visible that one more dimension is not
# available at either length.
#          k =  22  23  24  25  26  27  28  29  30
BKLC_D = {
    56: (14, 14, 13, 12, 12, 12, 12, 11, 10),
    57: (15, 14, 14, 12, 12, 12, 12, 12, 11),
    58: (16, 15, 14, 13, 12, 12, 12, 12, 12),
    59: (16, 16, 15, 14, 13, 12, 12, 12, 12),
    60: (16, 16, 16, 15, 14, 13, 12, 12, 12),
    61: (16, 16, 16, 16, 15, 14, 13, 12, 12),
    62: (16, 16, 16, 16, 16, 15, 14, 13, 12),
    63: (17, 16, 16, 16, 16, 16, 15, 14, 13),
    64: (18, 17, 16, 16, 16, 16, 16, 14, 14),
}
KMIN = 22


def best_k(n0):
    """the largest dimension whose best known minimum distance still reaches ceil(n0/4)"""
    need = -(-n0 // 4)
    ks = [KMIN + i for i, d in enumerate(BKLC_D[n0]) if d >= need]
    return (max(ks) if ks else None), need


if __name__ == '__main__':
    print(__doc__)
    print("=" * 88)
    print("1. the geometry, in exact rational arithmetic")
    print("=" * 88)
    for n0, d in ((24, 8), (32, 8), (62, 16), (63, 16)):
        ok, bad, worst = geometry(n0, d)
        check("n_0 = %-3d d = %-3d  ceil(n_0/4) = %-3d" % (n0, d, -(-n0 // 4)), ok,
              "worst inner product %s, at most 1/2" % worst)
    for n0, d in ((62, 15), (63, 15)):
        ok, bad, worst = geometry(n0, d)
        check("n_0 = %-3d d = %-3d  NEGATIVE CONTROL, must fail" % (n0, d), not ok,
              "fails at Hamming distance %s, as it must" % bad)

    print()
    print("=" * 88)
    print("2. the whole mechanism, end to end, on codes exhibited in coordinates")
    print("=" * 88)
    for G, name in ((golay24(), 'extended Golay'),
                    (rm(1, 5), 'Reed-Muller RM(1,5)'),
                    (rm(2, 5), 'Reed-Muller RM(2,5)')):
        k, n0 = G.shape[0], G.shape[1]
        C = codewords(G)
        d = min_distance(C)
        need = -(-n0 // 4)
        mx = worst_pair_linear(n0, d)
        ok = d >= need and mx * 2 <= n0
        note = ''
        if len(C) <= 4096:                     # small enough to also do it the slow way
            ex = worst_pair_exhaustive(C, n0)
            ok = ok and ex == mx
            note = ', brute force over all %d^2 pairs agrees' % len(C)
        check("%-20s [%d,%d,%d], needs d >= %-2d" % (name, n0, k, d, need), ok,
              "worst pair %d, at most n_0/2 = %d, so tau(%d) >= %d%s"
              % (mx, n0 // 2, n0, len(C), note))

    print()
    print("=" * 88)
    print("3. dimensions 62 and 63")
    print("=" * 88)
    CITED = {62: (26, 16), 63: (27, 16)}
    PREV = {62: 52417932, 63: 52418564}
    CAP = {62: 66075240, 63: 70543480}          # ../dim49-63-p48-caps/scripts/ceiling.py
    vals = {}
    for n0 in (62, 63):
        k, d = CITED[n0]
        need = -(-n0 // 4)
        ok, _, worst = geometry(n0, d)
        v = 1 << k
        vals[n0] = v
        check("dim %d: [%d,%d,%d] gives 2^%d = %d" % (n0, n0, k, d, k, v), ok and d >= need,
              "d = %d >= ceil(%d/4) = %d, worst inner product %s" % (d, n0, need, worst))
        check("dim %d: beats the previous value %d" % (n0, PREV[n0]), v > PREV[n0],
              "a factor %.2f" % (v / float(PREV[n0])))
        check("dim %d: beats the P_48 cap CEILING %d" % (n0, CAP[n0]), v > CAP[n0],
              "by %d, so that construction cannot reach this" % (v - CAP[n0]))
    print()
    print("=" * 88)
    print("4. the cited parameters, against the table around them")
    print("=" * 88)
    for n0 in sorted(BKLC_D):
        k, need = best_k(n0)
        star = ''
        if n0 in CITED:
            star = '   <- cited: [%d,%d,%d]' % (n0, CITED[n0][0], CITED[n0][1])
            check("n_0 = %d: largest k with d >= %d is %d, giving 2^%d = %d%s"
                  % (n0, need, k, k, 1 << k, star), k == CITED[n0][0],
                  "one more, k = %d, has d = %d < %d"
                  % (k + 1, BKLC_D[n0][k + 1 - KMIN], need))
        else:
            print("  ---- n_0 = %2d: needs d >= %2d, largest such k = %2d, gives 2^%-2d = %-11d"
                  % (n0, need, k, k, 1 << k))
    # nothing shorter does better for these two dimensions
    for n0 in (62, 63):
        alt = max((1 << best_k(m)[0]) for m in BKLC_D if m <= n0 and best_k(m)[0])
        check("dim %d: no length <= %d does better than 2^%d" % (n0, n0, CITED[n0][0]),
              alt == vals[n0], "best over all lengths in the table is %d" % alt)

    print()
    check("monotone: tau(62) <= tau(63)", vals[62] <= vals[63])
    check("monotone: tau(63) <= tau(64) = 331737984", vals[63] <= 331737984)

    print()
    print("=" * 88)
    if FAIL:
        print("*** %d FAILURES: %s ***" % (len(FAIL), ', '.join(FAIL)))
    else:
        print("ALL CHECKS PASSED -- tau(62) >= %d, tau(63) >= %d" % (vals[62], vals[63]))
    print("=" * 88)
    raise SystemExit(1 if FAIL else 0)
