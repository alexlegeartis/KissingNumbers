#!/usr/bin/env python3
"""Dimensions 62 and 63: the FULL Edel-Rains-Sloane chain, not only its level 0.

    python scripts/verify_chain.py        # about 20 s

verify.py establishes level 0 -- the 2^26 and 2^27 sign vectors of [62,26,16] and [63,27,16].
That is a complete kissing configuration on its own, and it is what this package used to claim.
It is not the whole of the Edel-Rains-Sloane construction.  Their count is

    N(n) = sum_nu  A(n, n_nu, n_nu) * A(n_nu, ceil(n_nu/4)),   n >= n_0 >= 4 n_1 >= 16 n_2 ...

and the levels above 0 are what carry their published dimension-64 total from 2^28 = 268 435 456
to 331 737 984.  Adding them here gives

    tau(62) >=  67 108 864 + 4 194 304 + 7 564 =  71 310 732     (+4 201 868)
    tau(63) >= 134 217 728 + 4 194 304 + 7 812 = 138 419 844     (+4 202 116)

EVERYTHING ADDED IS BUILT HERE, NOTHING NEW IS CITED.  The level-1 layer needs two codes and
both are constructed and their minimum distances enumerated in this file:

    a cyclic [15,6,8]_4      -> 4096 supports of weight 15 inside 60 of the coordinates,
                                pairwise meeting in at most 7, by the grid map
    a binary [15,10,4]       -> 1024 sign patterns on each such support

and the level-2 layer is all C(n,2) pairs of coordinates with all four sign patterns.  The two
level-0 code parameters [62,26,16] and [63,27,16] are cited exactly as before, and nothing else
is.  (codetables.de also lists [15,6,8]_4 as the best known quaternary code at those
parameters, so the grid map is not being handed anything unpublished either -- but it does not
have to be taken on trust, because the code is here.)

THE GRID MAP.  Index the cells of a w x q grid and let a support take one cell per row.  Two
supports then meet in exactly w - d_H of their choice words, so a q-ary code of length w and
minimum distance d gives supports meeting in at most w - d.  Level nu needs |S ^ S'| <= n_nu/2,
so d >= ceil(w/2) suffices.  Here w = 15, q = 4, d = 8 > 15/2, and 15 * 4 = 60 <= 62.

THE MECHANISM IS CHECKED END TO END IN COORDINATES, at the two dimensions where the chain
happens to reproduce the exact kissing number:

    n =  8, chain (8, 2)      240 points  = tau(8), the E_8 root system
    n = 16, chain (16, 4, 1)  4320 points = tau(16), the Barnes-Wall value and the record

Every pair of both configurations is checked directly.  Getting 240 and 4320 on the nose, from
a construction that knows nothing of E_8 or of BW_16, is the validation.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import numpy as np

FAIL = []


def check(cond, msg):
    print("   %-74s %s" % (msg, "ok" if cond else "*** FAILED ***"))
    if not cond:
        FAIL.append(msg)


print(__doc__)

# =========================================================================================
print("=" * 94)
print("1.  THE GEOMETRY, in exact rational arithmetic")
print("=" * 94)
# Level nu holds the vectors  a_nu * (+-1 on a support of size n_nu, 0 elsewhere), with
# a_nu^2 = n_0 / n_nu, so every point has squared norm n_0 and the 60-degree condition is
# <x,y> <= n_0/2.


def conditions(n0, chain):
    """(description, holds) for every pairwise condition the chain has to satisfy"""
    out = []
    for i, w in enumerate(chain):
        a2 = F(n0, w)                                   # a_nu^2
        d = -((-w) // 4)                                # ceil(w/4), the sign-code distance
        out.append(("level %d, same support, sign distance %d: %s <= %s"
                    % (i, d, a2 * (w - 2 * d), F(n0, 2)), a2 * (w - 2 * d) <= F(n0, 2)))
        m = w // 2
        out.append(("level %d, supports meeting in %d: %s <= %s"
                    % (i, m, a2 * m, F(n0, 2)), a2 * m <= F(n0, 2)))
        for j, v in enumerate(chain[i + 1:], start=i + 1):
            # levels i < j: <x,y> <= a_i a_j v, and (a_i a_j v)^2 = n0^2 * v / w
            lhs2 = F(n0) ** 2 * F(v, w)
            out.append(("levels %d and %d (%d >= 4*%d?): %s <= %s"
                        % (i, j, w, v, lhs2, F(n0, 2) ** 2), lhs2 <= F(n0, 2) ** 2))
    return out


for n0, chain in ((62, (62, 15, 2)), (63, (63, 15, 2))):
    print("   n_0 = %d, chain %s" % (n0, chain))
    for msg, ok in conditions(n0, chain):
        check(ok, "      " + msg)
check(not all(ok for _, ok in conditions(62, (62, 16, 2))),
      "NEGATIVE CONTROL: chain (62, 16, 2) violates a condition, since 62 < 4*16")

# =========================================================================================
print()
print("=" * 94)
print("2.  THE LEVEL-1 SUPPORT CODE: a cyclic [15,6,8]_4, built, distance enumerated")
print("=" * 94)
# GF(4) = {0,1,2,3}: addition is XOR of two bits, 2 = w and 3 = w^2 = w + 1.
MUL = np.zeros((4, 4), dtype=np.int8)
LOG = {1: 0, 2: 1, 3: 2}
ANTI = {0: 1, 1: 2, 2: 3}
for a in range(1, 4):
    for b in range(1, 4):
        MUL[a, b] = ANTI[(LOG[a] + LOG[b]) % 3]
# GF(16) = GF(2)[x]/(x^4 + x + 1), used only to build the generator polynomial from its roots
E, L, x = [0] * 32, [0] * 16, 1
for i in range(15):
    E[i], L[x] = x, i
    x <<= 1
    if x & 16:
        x ^= 0b10011
for i in range(15, 30):
    E[i] = E[i - 15]


def m16(a, b):
    return 0 if a == 0 or b == 0 else E[L[a] + L[b]]


EMB = {0: 0, 1: 1, E[5]: 2, E[10]: 3}          # GF(4) sitting inside GF(16)

DEFSET = [0, 1, 2, 3, 4, 5, 8, 10, 12]
cos, seen = [], set()
for s in range(15):
    if s in seen:
        continue
    c, y = [], s
    while y not in c:
        c.append(y)
        seen.add(y)
        y = (y * 4) % 15
    cos.append(sorted(c))
check(all(set(c) <= set(DEFSET) or not (set(c) & set(DEFSET)) for c in cos),
      "the defining set %s is a union of cyclotomic cosets mod 15 under x -> 4x" % DEFSET)

g16 = [1]
for i in DEFSET:
    r, ng = E[i % 15], [0] * (len(g16) + 1)
    for j, c in enumerate(g16):
        ng[j] ^= m16(c, r)
        ng[j + 1] ^= c
    g16 = ng
check(all(c in EMB for c in g16), "its generator polynomial has all coefficients in GF(4)")
g = [EMB[c] for c in g16]
K = 15 - len(DEFSET)
check(K == 6, "dimension k = 15 - %d = %d" % (len(DEFSET), K))

M = np.array(list(product(range(4), repeat=K)), dtype=np.int8)
C4 = np.zeros((len(M), 15), dtype=np.int8)
for j in range(K):
    for i, gc in enumerate(g):
        if gc:
            C4[:, j + i] ^= MUL[gc][M[:, j]]
wts = (C4 != 0).sum(axis=1)
DMIN = int(wts[wts > 0].min())
check(len(C4) == 4096, "4^6 = 4096 codewords enumerated")
check(DMIN == 8, "minimum distance over all 4096 codewords, by enumeration: %d" % DMIN)
check(DMIN >= -((-15) // 2), "8 >= ceil(15/2) = 8, so the grid map gives what level 1 needs")

# ---- the grid map: 4096 supports of weight 15 inside a 15 x 4 = 60 cell grid -------------
CELL = C4.astype(np.int64) + 4 * np.arange(15, dtype=np.int64)[None, :]
BITS = np.zeros(4096, dtype=object)
for r in range(15):
    for i in range(4096):
        BITS[i] |= 1 << int(CELL[i, r])
check(len(set(BITS.tolist())) == 4096, "the 4096 supports are distinct")
check(all(bin(b).count('1') == 15 for b in BITS), "every support has weight exactly 15")
worst = 0
for i in range(4096):
    bi = BITS[i]
    for j in range(i + 1, 4096):
        o = bin(bi & BITS[j]).count('1')
        if o > worst:
            worst = o
check(worst == 15 - DMIN,
      "largest overlap over all 8 386 560 pairs: %d = 15 - %d" % (worst, DMIN))
check(worst <= 15 // 2, "%d <= floor(15/2) = 7, which is what level 1 needs" % worst)

# =========================================================================================
print()
print("=" * 94)
print("3.  THE LEVEL-1 SIGN CODE: a binary [15,10,4], built, distance enumerated")
print("=" * 94)
Hm = np.array([[(j >> i) & 1 for j in range(1, 16)] for i in range(4)], dtype=np.int8)
Aw = Hm.copy() % 2
piv, r = [], 0
for c in range(15):
    pr = next((i for i in range(r, 4) if Aw[i, c]), None)
    if pr is None:
        continue
    Aw[[r, pr]] = Aw[[pr, r]]
    for i in range(4):
        if i != r and Aw[i, c]:
            Aw[i] ^= Aw[r]
    piv.append(c)
    r += 1
basis = []
for fc in [c for c in range(15) if c not in piv]:
    v = np.zeros(15, dtype=np.int8)
    v[fc] = 1
    for i, c in enumerate(piv):
        v[c] = Aw[i, fc]
    basis.append(v)
Bas = np.array(basis, dtype=np.int8)
check(len(Bas) == 11 and not (Hm.dot(Bas.T) % 2).any(),
      "a [15,11,3] Hamming code basis, in the kernel of its parity-check matrix")
Ext = np.hstack([Bas, (Bas.sum(axis=1) % 2)[:, None]])          # -> [16,11,4]
one = next((i for i in range(11) if Ext[i, 15] == 1), None)
Sh = []
for i in range(11):
    if i == one:
        continue
    v = Ext[i].copy()
    if v[15] == 1:
        v = (v + Ext[one]) % 2
    Sh.append(v[:15])
Sh = np.array(Sh, dtype=np.int8)
check(len(Sh) == 10, "extended to [16,11,4] and shortened once to [15,10,?]")
CW = np.zeros((1 << 10, 15), dtype=np.int8)
for j in range(10):
    CW ^= (((np.arange(1 << 10) >> j) & 1)[:, None] * Sh[j][None, :]).astype(np.int8)
w2 = CW.sum(axis=1)
D15 = int(w2[w2 > 0].min())
check(len(np.unique(CW, axis=0)) == 1024, "1024 codewords, all distinct")
check(D15 == 4, "minimum distance over all 1024 codewords, by enumeration: %d" % D15)
check(D15 >= -((-15) // 4), "4 >= ceil(15/4) = 4, which is what level 1 needs")

# =========================================================================================
print()
print("=" * 94)
print("4.  THE WHOLE MECHANISM, END TO END IN COORDINATES, WHERE THE ANSWER IS KNOWN")
print("=" * 94)


def build(n, n0, levels):
    """levels: (w, supports, sign code).  Integer points, all of squared norm n0."""
    P = []
    for w, sups, signs in levels:
        q = n0 // w
        assert q * w == n0 and int(round(q ** .5)) ** 2 == q, (n0, w)
        a = int(round(q ** .5))                          # a_nu = sqrt(n0/w), an integer here
        for S in sups:
            for s in signs:
                v = [0] * n
                for j, p in enumerate(S):
                    v[p] = a * (1 - 2 * int(s[j]))
                P.append(tuple(v))
    return P


def allsigns(w):
    """every binary vector of length w -- a valid sign code exactly when ceil(w/4) = 1"""
    assert -((-w) // 4) == 1
    return [tuple((m >> i) & 1 for i in range(w)) for m in range(1 << w)]


# ---- n = 8, chain (8, 2) ----------------------------------------------------------------
even8 = [tuple((m >> i) & 1 for i in range(8)) for m in range(256) if bin(m).count('1') % 2 == 0]
d8 = min(sum(a[i] ^ b[i] for i in range(8)) for a, b in combinations(even8, 2))
check(len(even8) == 128 and d8 == 2,
      "level 0 for n = 8: the even-weight [8,7,2] code, 128 words, d = 2 = ceil(8/4)")
P8 = build(8, 8, [(8, [tuple(range(8))], even8),
                  (2, list(combinations(range(8), 2)), allsigns(2))])
A8 = np.array(P8, dtype=np.int64)
check(len(set(P8)) == len(P8), "all %d points distinct" % len(P8))
check(set((A8 * A8).sum(axis=1).tolist()) == {8}, "every point has squared norm 8")
G8 = A8.dot(A8.T)
np.fill_diagonal(G8, -10 ** 9)
check(int(G8.max()) <= 4,
      "largest off-diagonal inner product %d <= 8/2 = 4, over all %d pairs"
      % (int(G8.max()), len(P8) * (len(P8) - 1) // 2))
check(len(P8) == 240, "TOTAL %d  ==  tau(8) = 240, the E_8 root system, exactly" % len(P8))

# ---- n = 16, chain (16, 4, 1) -----------------------------------------------------------
c16 = []
for m in range(1 << 15):
    v = np.array([(m >> i) & 1 for i in range(15)], dtype=np.int8)
    if not (Hm.dot(v) % 2).any():
        c16.append(tuple(np.concatenate([v, [v.sum() % 2]]).tolist()))
wt = np.array([sum(c) for c in c16])
check(len(c16) == 2048 and int(wt[wt > 0].min()) == 4,
      "level 0 for n = 16: the extended Hamming [16,11,4], 2048 words, d = 4 = ceil(16/4)")
flats = sorted({tuple(sorted((a, b, c, a ^ b ^ c)))
                for a in range(16) for b in range(16) for c in range(16)
                if len({a, b, c, a ^ b ^ c}) == 4})
check(len(flats) == 140,
      "level 1 supports: the %d affine 2-flats of AG(4,2), an S(3,4,16)" % len(flats))
check(max(len(set(u) & set(v)) for u, v in combinations(flats, 2)) <= 2,
      "they meet pairwise in at most 2 = floor(4/2)")
P16 = build(16, 16, [(16, [tuple(range(16))], c16),
                     (4, flats, allsigns(4)),
                     (1, [(p,) for p in range(16)], allsigns(1))])
A16 = np.array(P16, dtype=np.int64)
check(len(set(P16)) == len(P16), "all %d points distinct" % len(P16))
check(set((A16 * A16).sum(axis=1).tolist()) == {16}, "every point has squared norm 16")
G16 = A16.dot(A16.T)
np.fill_diagonal(G16, -10 ** 9)
check(int(G16.max()) <= 8,
      "largest off-diagonal inner product %d <= 16/2 = 8, over all %d pairs"
      % (int(G16.max()), len(P16) * (len(P16) - 1) // 2))
check(len(P16) == 4320,
      "TOTAL %d  ==  tau(16) = 4320, the Barnes-Wall value and the record, exactly" % len(P16))

# =========================================================================================
print()
print("=" * 94)
print("5.  DIMENSIONS 62 AND 63")
print("=" * 94)
LEVEL0 = {62: 1 << 26, 63: 1 << 27}          # [62,26,16] and [63,27,16], cited as before
TOTAL = {}
for n in (62, 63):
    l0, l1, l2 = LEVEL0[n], len(C4) * len(CW), comb(n, 2) * 4
    TOTAL[n] = l0 + l1 + l2
    print("   n = %d, chain (%d, 15, 2)" % (n, n))
    print("      level 0   A(%d,%d,%d) = 1      x A(%d,16) >= 2^%d           ->  %12d"
          % (n, n, n, n, 26 if n == 62 else 27, l0))
    print("      level 1   A(%d,15,15) >= 4096  x A(15,4)  >= 1024          ->  %12d" % (n, l1))
    print("      level 2   A(%d, 2, 2)  = %-5d x A( 2,1)   =    4          ->  %12d"
          % (n, comb(n, 2), l2))
    print("      %sTOTAL   tau(%d) >= %d" % (" " * 44, n, TOTAL[n]))
    print()
check(TOTAL[62] == 71310732,
      "tau(62) >= 71 310 732   (level 0 alone gives 67 108 864;  +4 201 868)")
check(TOTAL[63] == 138419844,
      "tau(63) >= 138 419 844  (level 0 alone gives 134 217 728; +4 202 116)")
for n in (62, 63):
    w3 = (n // 3) * ((n - 1) // 2) * 8       # a partial Steiner triple system, times A(3,1)
    check(comb(n, 2) * 4 > w3,
          "n = %d: level 2 at w = 2 gives %d, at w = 3 at most %d, so w = 2 is taken"
          % (n, comb(n, 2) * 4, w3))
check(4 * 16 > 63, "w = 15 is the largest level-1 weight the chain allows: 4*16 = 64 > 63")
check(TOTAL[62] > 70543480 and TOTAL[63] > 70543480,
      "both exceed the P_48 cap construction's absolute ceiling, 66 075 240 and 70 543 480")

print()
print("=" * 94)
if FAIL:
    print("*** %d CHECKS FAILED ***" % len(FAIL))
    raise SystemExit(1)
print("ALL CHECKS PASSED")
print("=" * 94)
