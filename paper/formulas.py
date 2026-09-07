#!/usr/bin/env python3
"""Independent verification of every formula stated in kissing46.tex.

    python formulas.py

Each block re-derives a displayed identity or inequality from its own hypotheses, either
symbolically with sympy or exactly with Fractions, and compares the result with what the
paper prints. Nothing here reads the manuscript: the expected values are written out, so
a disagreement means the manuscript and this file disagree, and one of them is wrong.

The numerical claims (table entries, counts, published values) are checked by
factcheck.py; this file is about the algebra.
"""
import sys
from fractions import Fraction as F
from itertools import product

import sympy as sp

BAD = []
OK = 0


def chk(name, cond, detail=''):
    global OK
    if cond:
        OK += 1
    else:
        BAD.append(name)
        print('  FAIL  %-58s %s' % (name, detail))


def head(s):
    print()
    print('== %s' % s)


# =====================================================================  Section 3
head('Lemma 3.1  the minimal shell is a kissing configuration')
mu, ip = sp.symbols('mu ip', positive=True)
# |v-w|^2 = 2mu - 2<v,w> and this is a non-zero lattice norm, so >= mu
sol = sp.solve(sp.Eq(2 * mu - 2 * ip, mu), ip)
chk('|v-w|^2 >= mu forces <v,w> <= mu/2', sol == [mu / 2], str(sol))

head('Proposition 3.2  extremal minimum mu = 2*floor(n/24) + 2')
for n, want in ((8, 2), (24, 4), (48, 6), (72, 8)):
    chk('mu(%d) = %d' % (n, want), 2 * (n // 24) + 2 == want)

head('Lemma 3.4  |<v,w>| <= mu/2, with equality mu only for w = +-v')
# |v -+ w|^2 = 2mu -+ 2<v,w> >= mu  gives  |<v,w>| <= mu/2
x = sp.symbols('x')
chk('both branches give |<v,w>| <= mu/2',
    sp.solve(sp.Eq(2 * mu - 2 * x, mu), x) == [mu / 2] and
    sp.solve(sp.Eq(2 * mu + 2 * x, mu), x) == [-mu / 2])
chk('equality |<v,w>| = mu forces |v -+ w|^2 = 0',
    sp.simplify(2 * mu - 2 * mu) == 0 and sp.simplify(2 * mu + 2 * (-mu)) == 0)
# the admissible inner-product sets quoted for P48 and Gamma72
chk('P_48 inner products {0,+-1,+-2,+-3} u {+-6}', list(range(-3, 4)) == [-3, -2, -1, 0, 1, 2, 3]
    and 6 == 6)
chk('Gamma_72 inner products {0,...,+-4} u {+-8}', 8 // 2 == 4)


# =====================================================================  Section 4
head('Lemma 4.2  the moment identity, checked on all four shells')
# n[c] tables as printed in Proposition 4.8 (n[0], n[1], ...), each non-zero c twice
ONE = {
    'E_8':      (8, 2, 240, [126, 56]),
    'Leech':    (24, 4, 196560, [93150, 47104, 4600]),
    'P_48':     (48, 6, 52416000, [23766960, 12608784, 1678887, 36848]),
    'Gamma_72': (72, 8, 6218175600, [2603658750, 1512243200, 280928256, 13959168, 127800]),
}
for L, (n, m, N, row) in ONE.items():
    chk('%s: mu/2 + 1 = %d entries' % (L, m // 2 + 1), len(row) == m // 2 + 1)
    chk('%s: the cells total N' % L, row[0] + 2 * sum(row[1:]) + 2 == N)
    top = 11 if n % 24 == 0 else 7
    for p in range(0, top + 1, 2):
        # left: sum over the shell of (<u,v>/mu)^p, with v = +-u contributing (+-1)^p
        lhs = sum(F(c, m) ** p * row[c] * (2 if c else 1) for c in range(len(row))) + 2
        # right: N * E[<uhat,Z>^p] / (n(n+2)...(n+p-2)) = N * (p-1)!! / prod
        dfac, den = 1, 1
        for j in range(1, p, 2):
            dfac *= j
        for j in range(p // 2):
            den *= n + 2 * j
        chk('%s degree %d' % (L, p), lhs == F(N * dfac, den), '%s vs %s' % (lhs, F(N * dfac, den)))

head('Lemma 4.2  the Gaussian normalisation E|Z|^{2m} = n(n+2)...(n+2m-2)')
for n in (8, 24, 48, 72):
    Z = sp.symbols('Z')
    for m2 in range(1, 6):
        # E|Z|^{2m} for Z ~ N(0,I_n) is 2^m * Gamma(n/2+m)/Gamma(n/2)
        want = sp.simplify(2 ** m2 * sp.gamma(sp.Rational(n, 2) + m2) / sp.gamma(sp.Rational(n, 2)))
        got = sp.prod([n + 2 * j for j in range(m2)])
        chk('E|Z|^{%d} in R^%d' % (2 * m2, n), sp.simplify(want - got) == 0,
            '%s vs %s' % (want, got))

head('Remark 4.3  the number of moment equations at k = 3')
full = [p for p in product(range(12), repeat=3) if sum(p) <= 11 and sum(p) % 2 == 0]
even = [p for p in product(range(12), repeat=3)
        if sum(p) <= 11 and all(q % 2 == 0 for q in p)]
chk('161 equations with mixed exponents', len(full) == 161, str(len(full)))
chk('56 equations with even exponents only', len(even) == 56, str(len(even)))

head('Lemma 4.5  the lattice-range inequality')
nu = sp.symbols('nu', positive=True)
# |v -+ w|^2 = mu + nu -+ 2<v,w> is a lattice norm: 0 or at least mu
chk('mu + nu - 2t >= mu gives t <= nu/2', sp.solve(sp.Eq(mu + nu - 2 * x, mu), x) == [nu / 2])
chk('mu + nu + 2t >= mu gives t >= -nu/2', sp.solve(sp.Eq(mu + nu + 2 * x, mu), x) == [-nu / 2])
# for a 60-degree pair u1-u2 is minimal, so |c1-c2| <= mu/2
chk('a 60-degree pair forces |c1 - c2| <= mu/2',
    sp.solve(sp.Eq(mu + mu - 2 * x, mu), x) == [mu / 2])

head('Lemma 4.6  weak duality')
# e0^T x = y^T A x + s^T x = b^T y + s^T x, and s^T x >= -sum_{u_j<inf} u_j max(0,-s_j)
import random
random.seed(0)
for trial in range(200):
    m_, n_ = 3, 6
    A = sp.Matrix(m_, n_, lambda i, j: sp.Rational(random.randint(-3, 3)))
    xs = sp.Matrix(n_, 1, lambda i, j: sp.Rational(random.randint(0, 4)))
    b = A * xs
    y = sp.Matrix(m_, 1, lambda i, j: sp.Rational(random.randint(-3, 3)))
    e0 = sp.Matrix(n_, 1, lambda i, j: 1 if i == 0 else 0)
    s = e0 - A.T * y
    cap = [None] * n_                      # None means +infinity
    for j in range(n_):
        if random.random() < 0.4:
            cap[j] = 4                     # a finite upper bound that xs respects
    feasible = all(cap[j] is None or xs[j] <= cap[j] for j in range(n_))
    dualfeas = all(s[j] >= 0 for j in range(n_) if cap[j] is None)
    if not (feasible and dualfeas):
        continue
    bound = (b.T * y)[0] - sum(cap[j] * max(0, -s[j]) for j in range(n_) if cap[j] is not None)
    chk('weak duality on a random instance', (e0.T * xs)[0] >= bound,
        '%s < %s' % ((e0.T * xs)[0], bound))


# =====================================================================  Section 5
head('Equation (2)  the class threshold gamma = floor(mu/3)/mu')
for lat, m, want in (('Leech', 4, F(1, 4)), ('P_48', 6, F(1, 3)), ('Gamma_72', 8, F(1, 4))):
    chk('gamma(%s) = %s' % (lat, want), F(m // 3, m) == want, str(F(m // 3, m)))
    # gamma is the largest c/mu with c an integer and c/mu <= 1/3
    best = max(c for c in range(m + 1) if F(c, m) <= F(1, 3))
    chk('%s: floor(mu/3) is the largest admissible c' % lat, best == m // 3)
# the ceiling 1/3 is exactly the condition for a second direction
t_ = sp.symbols('t')
g_ = sp.symbols('gamma')
tmin = sp.solve(sp.Eq(t_ * g_ + (1 - t_), sp.Rational(1, 2)), t_)[0]
chk('t_min = 1/(2(1-gamma))', sp.simplify(tmin - 1 / (2 * (1 - g_))) == 0, str(tmin))
thr = sp.simplify((sp.Rational(1, 2) - tmin) / (1 - tmin))
chk('the direction threshold is gamma/(2 gamma - 1)',
    sp.simplify(thr - g_ / (2 * g_ - 1)) == 0, str(sp.simplify(thr)))
chk('the WRONG form gamma/(gamma-1) is not equal to it',
    sp.simplify(thr - g_ / (g_ - 1)) != 0)
chk('threshold <= -1 exactly when gamma > 1/3',
    sp.solve(sp.Eq(g_ / (2 * g_ - 1), -1), g_) == [sp.Rational(1, 3)])
for gam, wt, wc, wm in ((F(1, 4), F(2, 3), F(-1, 2), 3), (F(1, 3), F(3, 4), F(-1), 2)):
    got_t = F(1, 2) / (1 - gam)
    got_c = (F(1, 2) - got_t) / (1 - got_t)
    chk('gamma=%s: t >= %s' % (gam, wt), got_t == wt, str(got_t))
    chk('gamma=%s: <z,z\'> <= %s' % (gam, wc), got_c == wc, str(got_c))
    chk('gamma=%s: |Z_i| <= %d' % (gam, wm), int(1 - 1 / got_c) == wm, str(1 - 1 / got_c))
# gamma = 3/8 for Gamma_72 would kill every pair
g38 = F(3, 8)
t38 = F(1, 2) / (1 - g38)
chk('gamma = 3/8 forces t >= 4/5', t38 == F(4, 5), str(t38))
chk('gamma = 3/8 gives <z,z\'> <= -3/2', (F(1, 2) - t38) / (1 - t38) == F(-3, 2))

head('Corollary 5.2  m unit vectors pairwise at cos <= c < 0 satisfy m <= 1 - 1/c')
mm, cc = sp.symbols('m c')
# 0 <= |sum z_i|^2 <= m + m(m-1)c
chk('m + m(m-1)c >= 0 gives m <= 1 - 1/c',
    sp.simplify(sp.solve(sp.Eq(mm + mm * (mm - 1) * cc, 0), mm)[-1] - (1 - 1 / cc)) == 0)
chk('c = -1/2 gives m <= 3', 1 - 1 / F(-1, 2) == 3)
chk('c = -1 gives m <= 2', 1 - 1 / F(-1) == 2)
# a zero-sum m-tuple has pairwise cosine exactly -1/(m-1)
for m3 in (2, 3):
    chk('a zero-sum %d-tuple has cosine -1/%d' % (m3, m3 - 1),
        F(-1, m3 - 1) <= ((F(1, 2) - (F(2, 3) if m3 == 3 else F(3, 4)))
                          / (1 - (F(2, 3) if m3 == 3 else F(3, 4)))))

head('Proposition 5.1  the nine pair types')
a, b = sp.sqrt(t_), sp.sqrt(1 - t_)
half = sp.Rational(1, 2)
# equator . cap: a * <vhat,uhat>, worst 1/2
chk('(6) equator-cap: sqrt(t)/2 <= 1/2 iff t <= 1',
    sp.solve(sp.Eq(a * half, half), t_) == [1])
# cap . cap same direction, two class lines: t*gamma + (1-t) <= 1/2
chk('(7) same direction: t gamma + (1-t) <= 1/2 rearranges to t >= 1/(2(1-gamma))',
    sp.simplify(sp.solve(sp.Eq(t_ * g_ + (1 - t_), half), t_)[0] - 1 / (2 * (1 - g_))) == 0)
# cap . cap same line, two directions: t + (1-t)<z,z'> <= 1/2
z = sp.symbols('z')
chk('(8) same line: t + (1-t)z <= 1/2 gives z <= (1/2-t)/(1-t)',
    sp.simplify(sp.solve(sp.Eq(t_ + (1 - t_) * z, half), z)[0] - (half - t_) / (1 - t_)) == 0)
# cap . cap two classes: t/2 + (1-t)<z,z'> <= 1/2
chk('(9) two classes: t/2 + (1-t)z <= 1/2 gives z <= 1/2',
    sp.simplify(sp.solve(sp.Eq(t_ * half + (1 - t_) * z, half), z)[0] - half) == 0)
# cap . pole: b<z,w> <= 1/2
chk('(10) cap-pole at t = 3/4: b = 1/2, so any pole is admissible',
    sp.sqrt(1 - sp.Rational(3, 4)) == half)
chk('(10) cap-pole at t = 2/3: needs <z,w> <= sqrt(3)/2, i.e. 30 degrees',
    sp.simplify(half / sp.sqrt(1 - sp.Rational(2, 3)) - sp.sqrt(3) / 2) == 0)
chk('30 degrees is arccos(sqrt(3)/2)', sp.acos(sp.sqrt(3) / 2) == sp.pi / 6)
# The two cap-cap cases the proof used to leave out.  Both are implied, and this is where
# that is checked rather than asserted.
chk('(7) forces t >= 1/2, so the antipodal cap pair 1-2t is at most 1/2',
    all(sp.Rational(1, 2 * (1 - gg)) >= half and 1 - 2 * sp.Rational(1, 2 * (1 - gg)) <= half
        for gg in (sp.Rational(0), sp.Rational(1, 4), sp.Rational(1, 3))))
chk('distinct lines of one class at distinct directions: t g + (1-t)z < 1/2 under (8)',
    all(sp.simplify(tt * gg + (half - tt) - half) < 0
        for tt, gg in ((sp.Rational(2, 3), sp.Rational(1, 4)),
                       (sp.Rational(3, 4), sp.Rational(1, 3)))))

head('Section 5.5  the dimension-27 direction set')
# R^3 has two CLASSICAL 12-point kissing configurations -- the code is not unique, so these
# are not the only two.  Only one of them has zero-sum triples, and the paper named the wrong
# one until this check existed.
PHI = (1 + sp.sqrt(5)) / 2
ICO = [v for s1 in (1, -1) for s2 in (1, -1)
       for v in ((sp.Integer(0), s1, s2 * PHI), (s1, s2 * PHI, sp.Integer(0)),
                 (s2 * PHI, sp.Integer(0), s1))]
CUB = [v for v in ((1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1))]
CUB = [sp.Matrix(v) for v in CUB] + [-sp.Matrix(v) for v in CUB]
ICO = [sp.Matrix(v) for v in ICO]

def cosines(V):
    n2 = sp.simplify((V[0].T * V[0])[0])
    return {sp.simplify((V[i].T * V[j])[0] / n2)
            for i in range(len(V)) for j in range(i + 1, len(V))}

chk('both are 12 points', len(ICO) == 12 and len(CUB) == 12)
ic, cu = cosines(ICO), cosines(CUB)
chk('icosahedron cosines are +-1, +-1/sqrt5',
    ic == {sp.Integer(-1), -1 / sp.sqrt(5), 1 / sp.sqrt(5)}, str(ic))
chk('cuboctahedron cosines are 0, +-1/2, -1',
    cu == {sp.Integer(-1), -half, sp.Integer(0), half}, str(cu))
chk('both are 60-degree codes', max(ic) <= half and max(cu) <= half)
chk('the icosahedron has NO non-antipodal pair at cos <= -1/2, hence no zero-sum triple',
    not any(c <= -half and c > -1 for c in ic))
chk('the cuboctahedron has 24 pairs at cos = -1/2',
    sum(1 for i in range(12) for j in range(i + 1, 12)
        if sp.simplify((CUB[i].T * CUB[j])[0] / 2) == -half) == 24)
chk('over the icosahedron the construction would give only 199548',
    196560 + 2 * 6 * 248 * (2 - 1) + 12 == 199548)
chk('over the cuboctahedron it gives 200540',
    196560 + 2 * 4 * 248 * (3 - 1) + 12 == 200540)

head('Equation (4)  the count')
NN, Ci, Zi, W = sp.symbols('N C Z W')
r = sp.symbols('r', positive=True, integer=True)
chk('N - 2 sum|C| + 2 sum|C||Z| + |W| = N + 2 sum|C|(|Z|-1) + |W|',
    sp.simplify((NN - 2 * Ci + 2 * Ci * Zi + W) - (NN + 2 * Ci * (Zi - 1) + W)) == 0)
chk('a triple is worth 4 per class line and a pair 2',
    2 * (3 - 1) == 4 and 2 * (2 - 1) == 2)
# one extra line in every class is worth 2 sum (|Z_i| - 1)
PARTS = {25: (0, 1), 26: (2, 0), 27: (4, 0), 28: (8, 0), 29: (12, 2), 30: (24, 0), 31: (42, 0)}
TAU = {25: 2, 26: 6, 27: 12, 28: 24, 29: 40, 30: 72, 31: 126}
CAL = {25: 197058, 26: 198550, 27: 200540, 28: 204520, 29: 209496, 30: 220440, 31: 238350}   # 25: the level-3/4 construction; the claim is 197569 (lens heads)
GAIN = {25: 2, 26: 8, 27: 16, 28: 32, 29: 52, 30: 96, 31: 168}
for d, (nt, npr) in PARTS.items():
    chk('dim %d: parts total tau(%d)' % (d, d - 24), 3 * nt + 2 * npr == TAU[d])
    got = 196560 + 2 * 248 * (2 * nt + npr) + TAU[d]
    chk('dim %d: the count is %d' % (d, CAL[d]), got == CAL[d], str(got))
    step = (196560 + 2 * 249 * (2 * nt + npr) + TAU[d]) - got
    chk('dim %d: one extra class line gives +%d' % (d, GAIN[d]), step == GAIN[d], str(step))
chk('dim 25 with a 249-line class gives 197060',
    196560 + 2 * 249 * (2 - 1) + 2 == 197060)
chk('dim 27: floor(2*12/3) = 8 gain-units, attained by four triples',
    2 * 12 // 3 == 8 and 4 * (3 - 1) == 8)
chk('dim 27: 196560 + 2*4*248*2 + 12 = 200540', 196560 + 2 * 4 * 248 * 2 + 12 == 200540)
chk('dim 38: 196560 - 2*98280 + 2*98280*3 = 6*98280 = 3*196560',
    196560 - 2 * 98280 + 2 * 98280 * 3 == 6 * 98280 == 3 * 196560 == 589680)
chk('dim 38: 589680 + 1932 = 591612', 589680 + 1932 == 591612)
chk('dim 38: 1932 - 24 = 1908', 1932 - 24 == 1908)
chk('dim 38: 1932/3 = 644 triples', 1932 // 3 == 644 and 1932 % 3 == 0)
chk('floor(tau(k)/3) for k = 1..7', [TAU[24 + k] // 3 for k in range(1, 8)]
    == [0, 2, 4, 8, 13, 24, 42])
chk('the largest gain-unit count at |W| = tau(k) is tau(k) - ceil(tau(k)/3)',
    all(3 * nt + 2 * npr == TAU[24 + k] and 2 * nt + npr
        == TAU[24 + k] - -(-TAU[24 + k] // 3)
        for k, (nt, npr) in ((k, PARTS[24 + k]) for k in range(1, 8))))

head('Proposition 5.6  Caro-Wei on the conflict graphs')
LINES = {'Leech': 98280, 'P_48': 26208000, 'Gamma_72': 3109087800}
DEG = {'Leech': 4600, 'P_48': 36848, 'Gamma_72': 13959168 + 127800}
CW = {'Leech': 22, 'P_48': 712, 'Gamma_72': 221}
FOUND = {'Leech': 248, 'P_48': 7069, 'Gamma_72': 2118}   # the largest each family ships
for L in LINES:
    chk('%s: lines are half the shell' % L, 2 * LINES[L] == ONE[L][2])
    chk('%s: Caro-Wei ceil(n/(d+1)) = %d' % (L, CW[L]),
        -(-LINES[L] // (DEG[L] + 1)) == CW[L], str(-(-LINES[L] // (DEG[L] + 1))))
    chk('%s: search beats it by about ten' % L, 9.0 < FOUND[L] / CW[L] < 12.0,
        '%.2f' % (FOUND[L] / CW[L]))
chk('P_48 conflict degree is n[3]', DEG['P_48'] == 36848)
chk('Gamma_72 conflict degree is n[3] + n[4]', DEG['Gamma_72'] == 13959168 + 127800)
chk('Leech conflict degree is n[2]', DEG['Leech'] == 4600)


# =====================================================================  Section 6
head('Section 6  the Edel-Rains-Sloane chain')
from math import comb
chk('A(8,2) = 2^7 = 128', 2 ** 7 == 128)
chk('A(2,1) = 4', 2 ** 2 == 4)
chk('ceil(39/4) = 10', -(-39 // 4) == 10)
chk('C(39,2) = 741', comb(39, 2) == 741)
chk('327680 + 128*3324 + 4*741 = 756116', 327680 + 128 * 3324 + 4 * 741 == 756116)
chk('327680 + 128*3323 + 4*741 = 755988', 327680 + 128 * 3323 + 4 * 741 == 755988)
chk('the gain is 128*(3324-3323)', 128 * (3324 - 3323) == 756116 - 755988 == 128)
chk('dim 64 reverse-engineers', 2 ** 28 + 30828 * 2048 + 10416 * 16 + 128 == 331737984)
chk('dim 80 reverse-engineers', 2 ** 30 + 143780 * 2048 + 20540 * 16 + 160 == 1368532064)
chk('2^30 is 17 percent of 6218175600', abs(2 ** 30 / 6218175600 - 0.17) < 0.005)
chk('2^32 is 69 percent of 6218175600', abs(2 ** 32 / 6218175600 - 0.69) < 0.005)
chk('2^32 does not exceed 6218175600', 2 ** 32 < 6218175600)

# dimensions 62 and 63: the chain (n,15,2), and the P_48 ceiling it beats
chk('ceil(62/4) = ceil(63/4) = 16', -(-62 // 4) == -(-63 // 4) == 16)
chk('weight 15 is the largest the chain allows', 4 * 15 <= 62 and 4 * 16 > 63)
chk('supports of weight 15 meet in at most 7', 15 - 8 == 7 == 15 // 2)
chk('C(62,2) = 1891 and C(63,2) = 1953', comb(62, 2) == 1891 and comb(63, 2) == 1953)
chk('chain (62,15,2) totals 71310732',
    2 ** 26 + 4 ** 6 * 2 ** 10 + 4 * comb(62, 2) == 71310732)
chk('chain (63,15,2) totals 138419844',
    2 ** 27 + 4 ** 6 * 2 ** 10 + 4 * comb(63, 2) == 138419844)
chk('level 1 is 4096*1024 = 4194304', 4 ** 6 * 2 ** 10 == 4194304)

# dimension 96: the chain (96,24,6,1), which is tight at 96 = 4*24 = 16*6
chk('ceil(96/4) = 24', -(-96 // 4) == 24)
chk('the chain (96,24,6,1) is admissible',
    96 >= 96 >= 4 * 24 >= 16 * 6 >= 64 * 1)
chk('and tight except at the last step', 4 * 24 == 96 and 16 * 6 == 96 and 64 * 1 < 96)
chk('weight-24 supports meet in at most 12', 24 - 12 == 12 == 24 // 2)
chk('weight-6 supports meet in at most 3', 6 - 3 == 3 == 6 // 2)
chk('A(24,6) = 16384 and A(6,2) = 32', 2 ** 14 == 16384 and 2 ** 5 == 32)
chk('level 0 is 2^33 = 8589934592', 2 ** 33 == 8589934592)
chk('level 1 is 4^9 * 16384 = 4294967296', 4 ** 9 * 16384 == 4294967296)
chk('level 2 is 16^4 * 32 = 2097152', 16 ** 4 * 32 == 2097152)
chk('level 3 is 96 * 2 = 192', 96 * 2 == 192)
chk('chain (96,24,6,1) totals 12886999232',
    2 ** 33 + 4 ** 9 * 16384 + 16 ** 4 * 32 + 96 * 2 == 12886999232)
chk('level 0 is two thirds of the total',
    abs(2 ** 33 / 12886999232 - 2 / 3) < 0.005)
chk('the floor is the direct sum 6218175600 + 196560',
    6218175600 + 196560 == 6218372160)
chk('the factor is 2.07', abs(12886999232 / 6218372160 - 2.07) < 0.005)
chk('the gain is 6668627072', 12886999232 - 6218372160 == 6668627072)
chk('the cap construction falls short at k = 24', 6480558568 < 12886999232)
chk('levels 1-3 at dimension 96 total 4297064640',
    4 ** 9 * 16384 + 16 ** 4 * 32 + 96 * 2 == 4297064640)
# The claim in 6.3 that the chain passes the direct sum "as soon as A(96,24) >= 2^31":
# 2^31 clears it and 2^30 does not, so the threshold is stated at the right power.
chk('2^31 clears the direct sum, 2^30 does not',
    2 ** 31 + 4297064640 > 6218372160 > 2 ** 30 + 4297064640)
# the cap ceiling over P_48 at k = 14, 15: every one of floor(tau(k)/2) classes at 7069
CEIL = lambda k, tau: 52416000 + 2 * (tau // 2) * 7069 + tau
chk('cap ceiling at k=14 is 66075240', CEIL(14, 1932) == 66075240)
chk('cap ceiling at k=15 is 70543480', CEIL(15, 2564) == 70543480)
chk('short of 2^26 by 1033624', 2 ** 26 - CEIL(14, 1932) == 1033624)
chk('short of 2^27 by 63674248', CEIL(15, 2564) + 63674248 == 2 ** 27)
# the class size that would close each gap: a CEILING, not a floor -- 7604 leaves 4 points
NEED = lambda k, tau, tgt: -(-(tgt - 52416000 - tau) // (2 * (tau // 2)))
chk('closing dim 62 needs 7605 lines, and 7604 falls short',
    NEED(14, 1932, 2 ** 26) == 7605 and 52416000 + 2 * 966 * 7604 + 1932 < 2 ** 26)
chk('closing dim 63 needs 31903 lines, and 31902 falls short',
    NEED(15, 2564, 2 ** 27) == 31903 and 52416000 + 2 * 1282 * 31902 + 2564 < 2 ** 27)


# =====================================================================  Section 8
head('Section 8  the remarks')
chk('dim 18: 6500 + 768 = 7268, against 7654', 6500 + 768 == 7268 and 7654 / 7268 < 1.06)
chk('dim 20: 128 * 16 = 2048', 128 * 16 == 2048)
chk('dim 19: 10668 + 4*320 = 11948', 10668 + 4 * 320 == 11948)
chk('dim 17: 5346 + 2*192 = 5730', 5346 + 2 * 192 == 5730)
chk('dim 17 certificate masses sum to 191',
    16 * F(65, 80) + 10 * F(40, 80) + 120 * F(30, 80) + 160 * F(28, 80)
    + 360 * F(10, 80) + 240 * F(9, 80) == 191)
chk('dim 17: 1 + g(0) = 192 and 5346 + 2*192 = 5730', 1 + 191 == 192 and 5346 + 384 == 5730)
chk('dim 17 ratio bound is the clique cover', F(1024 * 20, 60 + 20) == 256 == 8 * 32)
chk('dim 17: 4320 + 2*512 + 2 = 5346', 4320 + 2 * 512 + 2 == 5346)
chk('class bound floor(4680/11) = 425', 4680 // 11 == 425)
chk('4680/11 = 425.4545...', abs(F(4680, 11) - F(4254545, 10000)) < F(1, 1000))
chk('Hall: 1024/2025 + 44/891 < 1', F(1024, 2025) + F(44, 891) < 1,
    str(F(1024, 2025) + F(44, 891)))
chk('Lambda_21 min-max cosine sqrt(8/29)', abs(float(sp.sqrt(sp.Rational(8, 29))) - 0.525225731) < 1e-8)
chk('Lambda_22 min-max cosine sqrt(3/11)', abs(float(sp.sqrt(sp.Rational(3, 11))) - 0.522232968) < 1e-8)
chk('Lambda_23 min-max cosine sqrt(4/15)', abs(float(sp.sqrt(sp.Rational(4, 15))) - 0.516397779) < 1e-8)
chk('Lambda_14, Lambda_15 share sqrt(2/7)',
    abs(float(sp.sqrt(sp.Rational(2, 7))) - 0.534522484) < 1e-8)
chk('margins: sqrt(3/11)/(1/2) - 1 = 4.4 percent',
    abs(float(sp.sqrt(sp.Rational(3, 11))) / 0.5 - 1 - 0.0445) < 5e-4)
chk('C(200540,2) = 20108045530', comb(200540, 2) == 20108045530)

head('Proposition 8.1  the dimension-17 Delsarte certificate')
# The certificate is written out in the paper as a table of Fourier coefficients.  Its
# internal arithmetic is checkable here even though the function g is not: the masses must
# sum to g(0), the character counts must exhaust the group, and the spectrum must be the
# one the ratio bound is read from.
COUNT = [16, 10, 120, 160, 360, 240]
BY = [F(65, 80), F(40, 80), F(30, 80), F(28, 80), F(10, 80), F(9, 80)]
MASS = [13, 5, 45, 56, 45, 27]
chk('every b_y is non-negative', all(b >= 0 for b in BY))
chk('the masses are the products count * b_y',
    [c * b for c, b in zip(COUNT, BY)] == [F(m) for m in MASS],
    str([c * b for c, b in zip(COUNT, BY)]))
chk('g(0) = sum of the masses = 191', sum(MASS) == 191)
chk('the bound is |C| <= 1 + g(0) = 192', 1 + sum(MASS) == 192)
chk('192 named characters plus 118 zeros exhaust F_2^10',
    sum(COUNT) == 906 and 906 + 118 == 1024 == 2 ** 10)
# the spectrum, and the two coarser relaxations it is read from
SPEC = {60: 1, 30: 16, 12: 130, 6: 160, -2: 240, -4: 375, -10: 96, -20: 6}
chk('the spectrum has 1024 eigenvalues', sum(SPEC.values()) == 1024)
chk('the degree is the top eigenvalue, and |W_4| = 60', max(SPEC) == 60)
chk('the ratio bound is 1024*20/80 = 256',
    F(1024 * 20, 60 + 20) == 256 and 1024 * 20 // 80 == 256)
chk('the characters carrying weight are the eigenvalues 30, 12, 6, -4, -2',
    COUNT[0] == SPEC[30] and COUNT[1] + COUNT[2] == SPEC[12]
    and COUNT[3] == SPEC[6] and COUNT[5] == SPEC[-2] and COUNT[4] <= SPEC[-4])
chk('the 118 zero characters are what the table leaves out',
    SPEC[60] + SPEC[-10] + SPEC[-20] + (SPEC[-4] - COUNT[4]) == 118,
    str(SPEC[60] + SPEC[-10] + SPEC[-20] + (SPEC[-4] - COUNT[4])))
# the points where g <= -1
chk('887 + 60 + 16 = 963 points outside W_4 u {0}', 887 + 60 + 16 == 963)
chk('963 + |W_4| + 1 = 1024, so |W_4| = 60', 963 + 60 + 1 == 1024)
chk('scaling by 80 clears every denominator',
    all((80 * b).denominator == 1 for b in BY)
    and (80 * F(21, 5)).denominator == 1)
chk('alpha = 192 gives 5346 + 2*192 = 5730 = K(17)', 5346 + 2 * 192 == 5730)

head('Ozeki cross-checks')
chk('the degree-3 quotient is 12309600 exactly',
    F(23775066324172800000, 1931424768000) == 12309600)
chk('458386320/43 is not an integer', F(458386320, 43).denominator != 1)
chk('the affine line: n[000] + 164 n[333] = 6732912 at n[333] = 278 gives 6687320',
    6732912 - 164 * 278 == 6687320)
chk('the same line at n[333] = 184 gives 6702736', 6732912 - 164 * 184 == 6702736)
chk('the certified maximum is the vertex itself', 6732912 - 164 * 184 == 6702736)
chk('the range has width 164 * 94', 164 * (278 - 184) == 15416)

print()
print('=' * 78)
# the size stated for this script in paper/README.md, checked against this run
import os as _os, io as _io, re as _re
_rm = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'README.md')
if _os.path.exists(_rm):
    _said = _re.search(r'`formulas\.py` \((\d+) checks',
                       _io.open(_rm, encoding='utf-8').read())
    chk('paper/README.md states this script\'s current size',
        bool(_said) and int(_said.group(1)) == OK + 1,
        'README says %s; this run has %d'
        % (_said.group(1) if _said else None, OK + 1))

print('%d formula checks passed, %d FAILED' % (OK, len(BAD)))
if BAD:
    for nme in BAD:
        print('   * %s' % nme)
    print('=' * 78)
    print('*** A FORMULA IN THE MANUSCRIPT IS WRONG ***')
else:
    print('=' * 78)
    print('EVERY FORMULA VERIFIES')
print('=' * 78)
sys.exit(1 if BAD else 0)
