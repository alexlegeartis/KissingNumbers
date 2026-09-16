"""THE EQUATOR: its deep hole fixes the layer radii, and deleting from it buys nothing.

    python equator.py

Two statements about the 4320-point equator, both in exact integer arithmetic.

1.  THE DEEP-HOLE RULE.  In the scaling where everything has norm^2 = 8, a layer at radius r
    over the equator carries vectors x = (v, w) with |w|^2 = r^2 and |v|^2 = 8 - r^2.  Against
    an equator vector e (norm^2 = 8, w-part zero) the 60-degree condition is <v, e> <= 4, i.e.

        < v/|v| , e/|e| >  <=  4 / sqrt(8(8 - r^2)).

    So a layer at radius r exists only if the equator has a direction that far from all of its
    points, i.e. only if  m(E) := min_z max_e <z, e/|e|>  satisfies  m(E) <= 4/sqrt(8(8-r^2)),
    which rearranges to

        r^2  >=  8 - 2/m(E)^2 .

    An explicit direction is a certificate for the EASY side of this (a layer does exist), and
    both certificates are exact:

      * even Barnes-Wall shell, z = v/(2 sqrt 2) with v a +-1 vector on a six-set:
        |z|^2 = 3/4, so m <= 1/sqrt3 and r^2 >= 8 - 6 = 2 -- the TIER radius, and |v|^2 = 6 is
        exactly the v-part of a tier vector.  The deep holes ARE the tier directions.
      * odd Barnes-Wall shell, z = s sqrt2/6 with s in {+-1}^16:
        |z|^2 = 8/9, so m <= 3 sqrt2/8 and r^2 >= 8 - 64/9 = 8/9 -- the FLAT radius, and the
        layer vector's v-part is (2/3)(+-1)^16, which is the flat.

    That is the whole of Cohn-Li's sign flip: it deepens the equator's holes from 1/sqrt3 to
    3 sqrt2/8, 8.2%, and buys exactly the one extra layer that opens.

    NOT proved here: that those m are MINIMAL, i.e. that no layer sits closer in.  That is a
    minimax, measured by LP ascent (research/dim20-21/maximal2.py, research/noequator/) and
    reported as measured; the certificates below are one-sided.

2.  THE EQUATOR IS FREE FOR THE TIER.  Without an equator a tier of full parity classes could
    use any six-subset of [16], and would be 32 * A(16,8,6).

    Theorem.  A(16,8,6) = 16.  Let there be b blocks of size 6, pairwise meeting in <= 2, and
    let r_p count the blocks through p.  Then sum_p r_p = 6b and, for each block X,
    sum_{p in X} r_p - 6 = sum_{Y != X} |X ^ Y| <= 2(b-1), so sum_{p in X} r_p <= 2b + 4.
    Summing over blocks, sum_p r_p^2 <= b(2b+4), while Cauchy-Schwarz gives
    sum_p r_p^2 >= (6b)^2/16 = 9b^2/4.  Hence 9b^2/4 <= 2b^2 + 4b, i.e. b <= 16.  QED

    16 is attained, by the weight-6 words of a bent coset of RM(1,4) -- the 2-(16,6,2) biplane
    -- which is built and checked below.  So a tier is 512 with or without the equator: the
    4320 equator points cost a tier nothing, and no amount of equator deletion enlarges one.
"""
import os, sys, itertools

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from leech import load as leech_load
from golay import golay24

ok_all = []


def ok(cond, msg):
    ok_all.append(bool(cond))
    print(("  OK   " if cond else "  FAIL ") + msg)
    return bool(cond)


# ---------------------------------------------------------------- the two 4320-point shells
A = leech_load().astype(np.int64)
G = golay24()
O = G[G.sum(1) == 8][0]
oct_idx = np.nonzero(O)[0]
rest = [i for i in range(24) if O[i] == 0]
EVEN = A[(A[:, oct_idx] ** 2).sum(1) == 0][:, rest]        # norm^2 = 32 in these coordinates
ok(EVEN.shape == (4320, 16) and ((EVEN * EVEN).sum(1) == 32).all(),
   "even Barnes-Wall shell: 4320 vectors of norm^2 = 32 (Leech (^) octad complement)")

blockmask = (np.abs(EVEN) == 2).sum(1) == 8
ODD = EVEN.copy()
flip = np.nonzero(blockmask)[0]
for i in flip:                                              # make the minus count ODD
    ODD[i, np.nonzero(ODD[i])[0][0]] *= -1
ok((np.abs(ODD) == np.abs(EVEN)).all()
   and ((ODD[blockmask] < 0).sum(1) % 2 == 1).all()
   and ((EVEN[blockmask] < 0).sum(1) % 2 == 0).all(),
   "odd shell: the same 4320 supports, every (+-2^8) block vector given an ODD minus count")
for nm, S in (("even", EVEN), ("odd", ODD)):
    gram_max = -99
    for s in range(0, len(S), 500):
        C = S[s:s + 500] @ S.T
        for t in range(len(C)):
            C[t, s + t] = -99
        gram_max = max(gram_max, int(C.max()))
    ok(gram_max <= 16, "%s shell is a 60-degree code (max inner product %d <= 16 = 32/2)"
       % (nm, gram_max))

# ---------------------------------------------------------------- deep-hole certificates
# even: z = v/(2 sqrt2), v in {0,+-1}^16 of weight 6.  <z, e/|e|> = <v,e>/(2 sqrt2 * 4 sqrt2)
#       = <v,e>/16, so the condition <z, ehat> <= 1/2 is exactly the INTEGER <v,e> <= 8.
print("\neven shell, six-set directions:  <z, ehat> <= 1/2  <=>  <v, e> <= 8, v of weight 6")
six = []
for supp in itertools.combinations(range(16), 6):
    for m in range(64):
        v = np.zeros(16, dtype=np.int64)
        for t, p in enumerate(supp):
            v[p] = 1 - 2 * ((m >> t) & 1)
        six.append(v)
SIX = np.array(six)


def maxip(X, S):
    """max_e <x, e> for every row of X, in float32 -- every value is a small integer
    (|<x,e>| <= 6*4 = 24 here, 16*4 = 64 below), so the arithmetic is exact."""
    out = np.empty(len(X), dtype=np.int64)
    Sf = S.astype(np.float32)
    for s in range(0, len(X), 20000):
        M = X[s:s + 20000].astype(np.float32) @ Sf.T
        assert np.all(M == np.rint(M))
        out[s:s + 20000] = M.max(1).astype(np.int64)
    return out


best = maxip(SIX, EVEN)
good = SIX[best <= 8]
ok(len(good) > 0, "   %d of the %d weight-6 sign vectors have <v,e> <= 8 for all 4320 (the "
                  "deep holes)" % (len(good), len(SIX)))
supports = set(tuple(np.nonzero(v)[0].tolist()) for v in good)
ok(len(supports) == 448 and len(good) == 448 * 64,
   "   they are 448 six-sets carrying ALL 64 sign patterns each -- 448 is the number of "
   "weight-6 words of RM(2,4), and 28672 is the tier ground set of the dimension-18 work")
print("   one of them:", good[0].tolist())
print("   |z|^2 = 6/8 = 3/4, so m(even) <= 1/sqrt3 = 0.577350269 and a layer opens at")
print("   r^2 = 8 - 2/m^2 = 8 - 6 = 2: the TIER, whose v-part has norm 6 -- the six-set itself.")
ok(int((good[0] * good[0]).sum()) == 6, "   |v|^2 = 6 = 8 - 2, the tier's v-norm exactly")

# odd: z = s sqrt2/6, s in {+-1}^16.  <z, ehat> = <s,e> * (sqrt2/6) / (4 sqrt2) = <s,e>/24,
#      so the condition is the INTEGER <s,e> <= 12.
print("\nodd shell, all-+-1 directions:   <z, ehat> <= 1/2  <=>  <s, e> <= 12, s in {+-1}^16")
SGN = np.array([[1 - 2 * ((m >> t) & 1) for t in range(16)] for m in range(65536)],
               dtype=np.int64)
bo = maxip(SGN, ODD)
goodo = SGN[bo <= 12]
be = maxip(SGN, EVEN)
ok(len(goodo) == 2048, "   %d of the 65536 sign vectors have <s,e> <= 12 for all 4320 of the "
                       "ODD shell -- 2048 = |RM(2,4)|, the flat layer's word set" % len(goodo))
ok(int((be <= 12).sum()) == 0, "   and NONE of them does over the EVEN shell -- the flip is "
                               "what opens the layer")

# The two deep-hole sets are the same code, and naming it is checkable rather than asserted.
CODE = set(tuple(((1 - s) // 2).tolist()) for s in goodo)          # minus-sign indicators
# LINEARITY, completely: reduce to an F_2 basis, generate its whole span, compare the sets.
# (A sample of pairs would be a check that cannot fail on this data.)
rows = [list(c) for c in CODE]
piv, r = [], 0
for col in range(16):
    p = next((i for i in range(r, len(rows)) if rows[i][col]), None)
    if p is None:
        continue
    rows[r], rows[p] = rows[p], rows[r]
    for i in range(len(rows)):
        if i != r and rows[i][col]:
            rows[i] = [(a ^ b) for a, b in zip(rows[i], rows[r])]
    piv.append(col)
    r += 1
basis = [rows[i] for i in range(r)]
span = {(0,) * 16}
for b in basis:
    span |= set(tuple(x ^ y for x, y in zip(w, b)) for w in span)
closed = (span == CODE) and r == 11
wd = {}
for c in CODE:
    wd[sum(c)] = wd.get(sum(c), 0) + 1
ok(closed and len(CODE) == 2048 and wd == {0: 1, 4: 140, 6: 448, 8: 870, 10: 448, 12: 140, 16: 1},
   "   those 2048 words are LINEAR, of length 16 and dimension 11, with weight distribution "
   "1/140/448/870/448/140/1 -- so d = 4 and the code is the extended Hamming code of length "
   "16, which is unique: it is RM(2,4)")
ok(supports == set(tuple(np.nonzero(np.array(c))[0].tolist()) for c in CODE if sum(c) == 6),
   "   and the 448 six-sets that are deep holes of the EVEN shell are EXACTLY that code's "
   "weight-6 words -- one code, both shells")
print("   |z|^2 = 16*2/36 = 8/9, so m(odd) <= 3 sqrt2/8 = 0.530330086 and a layer opens at")
print("   r^2 = 8 - 2/m^2 = 8 - 64/9 = 8/9: the FLAT, whose v-part is (2/3)(+-1)^16.")
ok(abs(8 - 2 / (3 * 2 ** 0.5 / 8) ** 2 - 8.0 / 9) < 1e-12,
   "   8 - 2/(3 sqrt2/8)^2 = 8/9 exactly, and 16*(2/3)^2 = 64/9 = 8 - 8/9")
print("   measured (NOT proved here): the LP ascent of research/dim20-21/maximal2.py puts the")
print("   minimax at exactly these two values, 0.577350269 and 0.530330086.")

# ---------------------------------------------------------------- A(16,8,6) = 16
print("\nA(16,8,6) = 16: the bound is the two-line count above; here is the attainment.")


def evaluate(coeffs, deg2):
    """f(x) = sum c_S prod_{i in S} x_i over F_2^4, as a 16-bit truth table."""
    out = np.zeros(16, dtype=np.int8)
    for x in range(16):
        b = [(x >> i) & 1 for i in range(4)]
        val = coeffs[0]
        for i in range(4):
            val ^= coeffs[1 + i] * b[i]
        for (i, j) in deg2:
            val ^= b[i] * b[j]
        out[x] = val
    return out


bent = [evaluate([c0, c1, c2, c3, c4], [(0, 1), (2, 3)])
        for c0 in (0, 1) for c1 in (0, 1) for c2 in (0, 1) for c3 in (0, 1) for c4 in (0, 1)]
bent = np.array(bent)
w = bent.sum(1)
ok(sorted(set(w.tolist())) == [6, 10],
   "the coset x1x2 + x3x4 + RM(1,4) has 32 elements, all of weight 6 or 10 (bent)")
B = bent[w == 6]
ok(len(B) == 16, "   sixteen of them have weight 6 -- the blocks")
inter = B @ B.T
np.fill_diagonal(inter, -1)
ok(int(inter.max()) == 2 and int(inter[inter >= 0].min()) == 2,
   "   every two of those six-sets meet in EXACTLY 2 points: the 2-(16,6,2) biplane")
r = B.sum(0)
ok(sorted(set(r.tolist())) == [6], "   every point is in exactly 6 blocks, and 16*6 = 6*16")
print("   so A(16,8,6) = 16 and a tier of full parity classes is 32 * 16 = 512, equator or no")
print("   equator: dropping the 4320 to make room for more six-sets cannot help.")

print()
if all(ok_all):
    print("ALL CHECKS PASSED -- %d of %d" % (sum(ok_all), len(ok_all)))
else:
    print("FAILURES: %d of %d checks failed" % (len(ok_all) - sum(ok_all), len(ok_all)))
    sys.exit(1)
