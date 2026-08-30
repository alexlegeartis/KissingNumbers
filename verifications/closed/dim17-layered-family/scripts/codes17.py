"""The dimension-17 tiers as BINARY CODES, and how much room the constraints really leave.

In the integer frame (equator norm^2 8, tier A norm^2 6, tier B norm^2 16) every layering
inequality becomes a small integer bound:

    e.e' <= 4      a.a' <= 2      b.b' <= 7      e.a <= 4      e.b <= 6      a.b <= 4

Tier B consists of all-(+-1) vectors, so b.b' = 16 - 2d with d the Hamming distance of the sign
patterns: b.b' <= 7 means distance >= 5, and since we will see the patterns live in an
even-distance coset, distance >= 6.  Cohn-Li use 192 of them.  A(16,6) = 256 (Nordstrom-
Robinson), so the question is whether the OTHER constraints -- against the equator and against
tier A -- are what cost the missing 64, or whether 192 was simply not optimised.

This script identifies:
  * which 8-subsets carry the equator's 1^8 vectors (they should be the 30 octads of RM(1,4)),
    and the sign parity class used on each;
  * tier B's sign patterns as a binary code: size, distance distribution, and whether it sits
    in a coset of RM(2,4) = the [16,11,4] extended Hamming code (which is what the equator
    constraint forces);
  * tier A's supports and per-support sign codes.
"""
import os, sys, collections, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

H = os.path.dirname(os.path.abspath(__file__))
E = np.load(os.path.join(H, '..', 'data', 'd17_equator.npy'))
import glob
tf = sorted(glob.glob(os.path.join(H, '..', 'data', 'd17_tier_*.npy')))
A = np.load([f for f in tf if '0p5' in f][0])
B = np.load([f for f in tf if '0p33' in f][0])
Ei = np.round(E * np.sqrt(8)).astype(np.int64)
Ai = np.round(A * np.sqrt(6)).astype(np.int64)
Bi = np.round(B * 4).astype(np.int64)
print(__doc__)
assert set((Ei * Ei).sum(1).tolist()) == {8} and set((Ai * Ai).sum(1).tolist()) == {6}
assert set((Bi * Bi).sum(1).tolist()) == {16}
print("integer frames verified: e.e=8, a.a=6, b.b=16")

# --- equator: the octads --------------------------------------------------
oct_rows = Ei[np.abs(Ei).sum(1) == 8]
supports = collections.Counter(tuple(np.nonzero(np.abs(r))[0].tolist()) for r in oct_rows)
print("\nequator 1^8 vectors: %d, on %d distinct 8-subsets, %s vectors each"
      % (len(oct_rows), len(supports), sorted(set(supports.values()))))
OCT = sorted(supports)
par = {}
for O in OCT:
    ps = {int(sum(1 for i in O if r[i] < 0) % 2) for r in oct_rows
          if tuple(np.nonzero(np.abs(r))[0].tolist()) == O}
    par[O] = ps
print("sign parities used per octad: %s" % sorted({tuple(sorted(v)) for v in par.values()}))

# --- tier B as a binary code ----------------------------------------------
beta = ((1 - Bi) // 2).astype(np.uint8)                 # +1 -> 0, -1 -> 1
words = sorted({int("".join(map(str, r)), 2) for r in beta})
print("\ntier B: %d vectors -> %d distinct sign patterns" % (len(Bi), len(words)))
dd = collections.Counter()
for i in range(len(words)):
    for j in range(i + 1, len(words)):
        dd[bin(words[i] ^ words[j]).count('1')] += 1
print("tier B distance distribution: %s" % dict(sorted(dd.items())))
print("tier B minimum distance: %d" % min(dd))

# RM(1,4) from the octads, then RM(2,4) = its dual
rm1 = {0}
for O in OCT:
    v = 0
    for i in O:
        v |= 1 << i
    rm1 |= {x ^ v for x in rm1}
print("\nspan of the octads: %d words (RM(1,4) has 32)" % len(rm1))
rm2 = [x for x in range(1 << 16)
       if all(bin(x & y).count('1') % 2 == 0 for y in rm1)]
print("its dual RM(2,4): %d words, min weight %d"
      % (len(rm2), min(bin(x).count('1') for x in rm2 if x)))
diffs = {words[i] ^ words[j] for i in range(len(words)) for j in range(len(words))}
print("tier B differences all inside RM(2,4)? %s" % diffs.issubset(set(rm2)))

# --- tier A ----------------------------------------------------------------
supA = collections.Counter(tuple(np.nonzero(r)[0].tolist()) for r in Ai)
print("\ntier A: %d vectors on %d distinct 6-subsets, %s vectors each"
      % (len(Ai), len(supA), sorted(set(supA.values()))))
mx = 0
for S in supA:
    for O in OCT:
        mx = max(mx, len(set(S) & set(O)))
print("maximum |support(a) n octad| over tier A: %d   (bound is 4)" % mx)
