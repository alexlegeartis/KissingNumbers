"""THE LAYER IDENTITY, dimensions 17 to 24 -- verified as a STRUCTURE, not fitted to a table.

    python layers.py

Fix an octad O of the binary Golay code and split the 24 Leech coordinates as
R^24 = R^16 (+) R^8, the R^8 being O's eight coordinates.  Write p(x) for the R^8 part of a
Leech minimal vector x (norm^2 = 32 throughout).  This script establishes, by exhaustive
enumeration in exact integer arithmetic:

  (L1)  |p(x)|^2 takes only the values 0, 8, 16, 32 -- never 24, and never anything else;
  (L2)  the set of values p(x) takes is exactly
            {0}                                        (1 position)
            the 240 roots of 2E8                       (|p|^2 = 8)
            the 2160 norm-4 vectors of 2E8             (|p|^2 = 16)
            the 240 DOUBLED roots of 2E8               (|p|^2 = 32)
        -- note the last line: of the 17520 vectors of 2E8 of norm 32, only the 240 that are
        twice a root occur;
  (L3)  the multiplicity is CONSTANT on each of those four sets: 4320, 512, 32, 1.

(L3) is the whole content.  It makes the kissing number of every coordinate section additive
over positions, so for EVERY subspace W of the R^8,

    tau( Leech (^) (R^16 (+) W) )  =  4320 + 512 a(W) + 32 b(W) + c(W)

with a, b, c the number of roots, norm-4 vectors and doubled roots of 2E8 lying in W.  When W
is spanned by roots and L = E8 (^) W is the root lattice they generate, c = a = N_2(L) and
b = N_4(L), which is the identity

    tau(Lambda_{16+k})  =  4320 + 513 N_2(L_k) + 32 N_4(L_k),
    L_1..L_8 = A1, A2, A3, D4, D5, E6, E7, E8.

The script then evaluates both sides for the standard chain A1 < A2 < A3 < D4 < D5 < E6 < E7
< E8 of sub-root-systems -- by the formula AND by counting the section directly -- and checks
both against the published kissing numbers of the laminated lattices.

The "513" is 512 + 1: a tier and a pole, which sit at DIFFERENT positions (a root and twice
that root) and are grouped only because both sets are indexed by the roots.  Writing it as one
coefficient is an accident of counting, not a fact about the geometry.
"""
import os, sys, random, itertools, collections
from fractions import Fraction

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


# ---------------------------------------------------------------- 2E8, from its definition
def two_e8_upto(maxnorm):
    """Every x in 2E8 with |x|^2 <= maxnorm.

    2E8 = 2*(D8 u (D8 + (1/2)^8)) = {x in Z^8 : x all even, sum(x)/2 even}
                                  u {x in Z^8 : x all odd,  sum(x) = 0 mod 4}.
    Enumerated by pruned depth-first search, so nothing about the answer is assumed.
    """
    out = []
    for parity in (0, 1):
        step = [v for v in range(-int(maxnorm ** 0.5) - 1, int(maxnorm ** 0.5) + 2)
                if v % 2 == parity or (v % 2 + 2) % 2 == parity]
        step = [v for v in step if abs(v) % 2 == parity and v * v <= maxnorm]

        def rec(d, rem, cur, tot):
            if d == 8:
                if (parity == 0 and (tot // 2) % 2 == 0) or (parity == 1 and tot % 4 == 0):
                    out.append(tuple(cur))
                return
            for v in step:
                if v * v <= rem:
                    cur.append(v)
                    rec(d + 1, rem - v * v, cur, tot + v)
                    cur.pop()

        rec(0, maxnorm, [], 0)
    return out


E8SHELL = collections.defaultdict(set)
for _v in two_e8_upto(32):
    E8SHELL[sum(t * t for t in _v)].add(_v)
print("2E8 built from its definition (all-even with half-sum even, or all-odd with sum = 0 mod 4)")
print("   shell norms <= 32 and their sizes: %s"
      % {k: len(E8SHELL[k]) for k in sorted(E8SHELL) if k})
ok(len(E8SHELL[8]) == 240, "2E8 has 240 roots (norm^2 = 8)")
ok(len(E8SHELL[16]) == 2160, "2E8 has 2160 vectors of norm^2 = 16")
ok(len(E8SHELL[32]) == 17520, "2E8 has 17520 vectors of norm^2 = 32")
DOUBLED = set(tuple(2 * t for t in r) for r in E8SHELL[8])
ok(len(DOUBLED) == 240 and DOUBLED <= E8SHELL[32],
   "240 of those 17520 are twice a root")

# ---------------------------------------------------------------- Leech, split at an octad
A = leech_load().astype(np.int64)
ok(A.shape == (196560, 24) and ((A * A).sum(1) == 32).all(),
   "Leech minimal vectors: 196560, all of norm^2 = 32")
G = golay24()
OCTADS = G[G.sum(1) == 8]
ok(OCTADS.shape[0] == 759, "the Golay code has 759 octads")

print("\n(L1)-(L3), over %d octads (the split is O's eight coordinates):" % 12)
PROFILE = {0: (1, 4320), 8: (240, 512), 16: (2160, 32), 32: (240, 1)}
rng = random.Random(7)
sample = [0] + rng.sample(range(1, 759), 11)
prof_ok, set_ok = True, True
for t in sample:
    idx = np.nonzero(OCTADS[t])[0]
    P = A[:, idx]
    q = (P * P).sum(1)
    norms = sorted(set(q.tolist()))
    cnt = collections.Counter(map(tuple, P.tolist()))
    seen = collections.defaultdict(collections.Counter)
    for v, m in cnt.items():
        seen[sum(x * x for x in v)][m] += 1
    if norms != [0, 8, 16, 32]:
        prof_ok = False
    for nq, (npos, mult) in PROFILE.items():
        if dict(seen[nq]) != {mult: npos}:
            prof_ok = False
    pos = {nq: set(v for v in cnt if sum(x * x for x in v) == nq) for nq in norms}
    if not (pos[0] == {(0,) * 8} and pos[8] == E8SHELL[8]
            and pos[16] == E8SHELL[16] and pos[32] == DOUBLED):
        set_ok = False
ok(prof_ok, "(L1)+(L3) |p|^2 in {0,8,16,32} only, with 1x4320, 240x512, 2160x32, 240x1")
ok(set_ok, "(L2) the four position sets are {0}, the 2E8 roots, its norm-16 shell, "
           "and the doubled roots")
print("   (24 never occurs: a (-+3,+-1^23) vector puts 8 or 16 on an octad, a (+-2^8) puts")
print("    4|B ^ O| in {0,8,16,32} since two octads meet in 0, 2, 4 or 8 points, and a")
print("    (+-4^2) puts 0, 16 or 32.)")

# the equator itself
idx0 = np.nonzero(OCTADS[0])[0]
rest = [i for i in range(24) if i not in set(idx0.tolist())]
EQ = A[(A[:, idx0] ** 2).sum(1) == 0]
shape4 = int((np.abs(EQ) == 4).sum(1).astype(int).__eq__(2).sum())
shape2 = int((np.abs(EQ) == 2).sum(1).astype(int).__eq__(8).sum())
ok(len(EQ) == 4320 and shape4 == 480 and shape2 == 3840,
   "the |p| = 0 layer is 4320 vectors: 480 of shape (+-4^2) and 3840 of shape (+-2^8)")
print("   that layer is the Barnes-Wall lattice Lambda_16 = BW16 (Conway-Sloane, SPLAG ch. 4:")
print("   Leech (^) (complement of an octad) = Lambda_16), and 4320 is its kissing number.")


# ---------------------------------------------------------------- sections by a subspace
def complement_basis(rows):
    """Integer basis of {x : rows . x = 0}, exactly, over Q.  Returns (rank, basis)."""
    M = [[Fraction(t) for t in r] for r in rows]
    piv, r = [], 0
    for c in range(8):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        inv = M[r][c]
        M[r] = [t / inv for t in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c)
        r += 1
    free = [c for c in range(8) if c not in piv]
    N = []
    for f in free:
        v = [Fraction(0)] * 8
        v[f] = Fraction(1)
        for i, c in enumerate(piv):
            v[c] = -M[i][f]
        den = 1
        for t in v:
            den = den * t.denominator // np.gcd(den, t.denominator)
        N.append([int(t * den) for t in v])
    return r, np.array(N, dtype=np.int64).reshape(len(N), 8)


R8 = np.array(sorted(E8SHELL[8]), dtype=np.int64)
R16 = np.array(sorted(E8SHELL[16]), dtype=np.int64)
R32 = np.array(sorted(DOUBLED), dtype=np.int64)


def abc(rows):
    rank, N = complement_basis(rows)
    if N.size == 0:
        return rank, len(R8), len(R16), len(R32), None
    inW = lambda S: (S @ N.T == 0).all(1)
    return rank, int(inW(R8).sum()), int(inW(R16).sum()), int(inW(R32).sum()), N


# Bourbaki's simple roots of E8, written in the 2E8 coordinates (twice the usual ones)
SIMPLE = {1: [1, -1, -1, -1, -1, -1, -1, 1], 2: [2, 2, 0, 0, 0, 0, 0, 0],
          3: [-2, 2, 0, 0, 0, 0, 0, 0], 4: [0, -2, 2, 0, 0, 0, 0, 0],
          5: [0, 0, -2, 2, 0, 0, 0, 0], 6: [0, 0, 0, -2, 2, 0, 0, 0],
          7: [0, 0, 0, 0, -2, 2, 0, 0], 8: [0, 0, 0, 0, 0, -2, 2, 0]}
ok(all(tuple(v) in E8SHELL[8] for v in SIMPLE.values()),
   "all eight Bourbaki simple roots are roots of 2E8")

# nested sub-diagrams: a1 ... a8 with a2 attached to a4 and the chain a3-a4-a5-a6-a7-a8
CHAIN = [(1, 'A1', [3]), (2, 'A2', [3, 4]), (3, 'A3', [3, 4, 5]), (4, 'D4', [3, 4, 5, 2]),
         (5, 'D5', [3, 4, 5, 2, 6]), (6, 'E6', [3, 4, 5, 2, 6, 1]),
         (7, 'E7', [3, 4, 5, 2, 6, 1, 7]), (8, 'E8', [3, 4, 5, 2, 6, 1, 7, 8])]
# Kissing numbers of the laminated lattices, Conway-Sloane, SPLAG Table 6.1 (Lambda_16..Lambda_24).
TAU_LAMBDA = {17: 5346, 18: 7398, 19: 10668, 20: 17400, 21: 27720,
              22: 49896, 23: 93150, 24: 196560}
# Published kissing-number records, for the difference in the last column.
RECORD = {17: 5730, 18: 8358, 19: 11948, 20: 19448, 21: 29768,
          22: 49896, 23: 93150, 24: 196560}

print("\nthe chain A1 < A2 < A3 < D4 < D5 < E6 < E7 < E8, by the formula AND by counting:")
print("   k dim type rank  N_2  N_4  doubled   4320+512a+32b+c   counted   tau(Lambda)  record-tau")
P0 = A[:, idx0]
form_ok, count_ok = True, True
for k, nm, gens in CHAIN:
    rank, a, b, c, N = abc([SIMPLE[g] for g in gens])
    val = 4320 + 512 * a + 32 * b + c
    inW = np.ones(len(A), bool) if N is None else (P0 @ N.T == 0).all(1)
    direct = int(inW.sum())
    if rank != k or val != TAU_LAMBDA[16 + k]:
        form_ok = False
    if direct != TAU_LAMBDA[16 + k]:
        count_ok = False
    print("   %d %3d %-4s  %d   %4d %5d   %4d   %14d  %8d   %9d   %6d"
          % (k, 16 + k, nm, rank, a, b, c, val, direct, TAU_LAMBDA[16 + k],
             RECORD[16 + k] - TAU_LAMBDA[16 + k]))
ok(form_ok, "the formula reproduces tau(Lambda_17) .. tau(Lambda_24), all eight")
ok(count_ok, "and so does a DIRECT count of the Leech vectors whose p lies in W")
ok(all(abc([SIMPLE[g] for g in gens])[3] == abc([SIMPLE[g] for g in gens])[1]
       for _, _, gens in CHAIN), "c = a on every one of them, which is why 512 + 1 = 513")

# ---------------------------------------------------------------- is the chain the best W?
print("\nis any other root-spanned section better?  greedy from a random root, 2000 restarts:")
best = {}
Rl = R8.tolist()
rng2 = random.Random(1)
for _ in range(2000):
    B = []
    for t in rng2.sample(range(240), 8):
        cand = B + [Rl[t]]
        rank, a, b, c, _N = abc(cand)
        if rank != len(cand):
            continue
        B = cand
        v = 4320 + 512 * a + 32 * b + c
        if v > best.get(rank, 0):
            best[rank] = v
beat = [k for k in range(1, 9) if best.get(k, 0) > TAU_LAMBDA[16 + k]]
for k in range(1, 9):
    print("   rank %d: best found %8d,  tau(Lambda_%d) = %8d   %s"
          % (k, best.get(k, 0), 16 + k, TAU_LAMBDA[16 + k],
             "equal" if best.get(k, 0) == TAU_LAMBDA[16 + k] else "DIFFERENT"))
ok(not beat, "no restart beat the laminated chain at any rank (a SEARCH, not a proof; "
             "that max #roots in a rank-k subspace of E8 is 2,6,12,24,40,72,126,240 is "
             "classical -- Dynkin 1952, Borel-de Siebenthal 1949)")

print()
if all(ok_all):
    print("ALL CHECKS PASSED -- %d of %d" % (sum(ok_all), len(ok_all)))
else:
    print("FAILURES: %d of %d checks failed" % (len(ok_all) - sum(ok_all), len(ok_all)))
    sys.exit(1)
