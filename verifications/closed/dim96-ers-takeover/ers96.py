#!/usr/bin/env python3
"""Edel-Rains-Sloane at n = 93..96, with the code tables as of 2026-08-23.

    tau_n >= N(n) = sum_nu  A(n, n_nu, n_nu) * A(n_nu, ceil(n_nu/4))
    over an admissible chain  n >= n_0 >= 4 n_1 >= 16 n_2 >= 64 n_3 >= ...

The formula, its provenance and its reproduction of ERS's own published totals at n = 64 and
n = 80 are in ../dim32-44-ers-audit/, whose PIPELINE is imported here for the binary code
table A(n,d) (n <= 60).

WHICH WAY EACH INEQUALITY POINTS.  Everything this script computes is a LOWER bound on what
ERS gives, because every input is a construction.  That settles n = 96 -- a lower bound on
ERS above Gamma_72 proves ERS wins -- and it does NOT settle n = 93, 94, 95, where a lower
bound below Gamma_72 proves nothing.  For those three the script instead reports the
THRESHOLD: how large the level-1 support code would have to be for ERS to overtake this
repository's claim, and how much room the Johnson bound leaves above the best construction
in hand.  See the "exposure" block at the end; it is the honest statement of the risk and it
is why 93-95 are printed at all.

TWO INPUTS COME FROM OUTSIDE PIPELINE'S TABLE.

  1.  The level-0 sign code A(n_0, ceil(n_0/4)) for n_0 in 93..96, read off codetables.de
      (Grassl, best known LINEAR codes) on 2026-08-26.  Exactly the source, and exactly the
      quantity, that ERS themselves used: codetables gives [64,28,16], and 2^28 is ERS's
      published level-0 input at n = 64.

          n_0 = 93  [93,31,24]   A >= 2^31          n_0 = 95  [95,32,24]   A >= 2^32
          n_0 = 94  [94,31,24]   A >= 2^31          n_0 = 96  [96,33,24]   A >= 2^33

      **n_0 = 96 is where this doubles**, and that single step is what puts dimension 96
      out of reach of anything built on Gamma_72.

  2.  The level-nu SUPPORT codes A(n, w, w), which no table covers at n = 93..96.  Lower
      bounds here come from the q-ary grid map: split n = w * q with q a prime power, index
      the cells of a w x q grid, and let each support pick one cell per row.  Two supports
      meet in exactly (w - Hamming distance of their choice words), so a q-ary code of
      length w and minimum distance >= w/2 gives  A(n, w, w) >= q^k, which is what level nu
      needs (|S ^ S'| <= w/2).  The q-ary codes are codetables.de again:

          96 = 24 x 4    [24,9,12]_4    A(96,24,24) >= 4^9 = 262144
          96 =  6 x 16   RS [6,4,3]_16  A(96, 6, 6) >= 16^4 = 65536   (built and checked below)
          92 = 23 x 4    [23,8,12]_4    A(n,23,23)  >= 4^8 =  65536   for n = 93, 94, 95
"""
import os, sys
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'dim32-44-ers-audit'))
import PIPELINE as P

GAMMA72 = 6218175600
# this repository's own claims, the values ERS would actually have to beat
# The claims come from RESULTS.md, never from a copy here.  They were hardcoded until
# 2026-08-26 and went stale the moment the class families grew: this file still read
# 6402983248 at dimension 95 after the claim had moved to 6472821892, so the threshold
# and the exposure factor it printed were both computed against a superseded claim.
# ers_exposure.claims() already reads RESULTS.md; use it rather than keep a second copy.
from ers_exposure import claims as _claims, acw as _acw
CLAIMED = {n: v for n, v in _claims().items() if n in (93, 94, 95)}
assert set(CLAIMED) == {93, 94, 95}, CLAIMED

# --------------------------------------------------------------------- verify the grid map
# The one grid code that can be built and checked here from scratch: RS[6,4,3] over GF(16).
def gf16_tables():
    exp, log, x = [0] * 32, [0] * 16, 1
    for i in range(15):
        exp[i] = x
        log[x] = i
        x <<= 1
        if x & 16:
            x ^= 0b10011                      # x^4 + x + 1
    for i in range(15, 30):
        exp[i] = exp[i - 15]
    return exp, log

EXP, LOG = gf16_tables()
def mul(a, b):
    return 0 if a == 0 or b == 0 else EXP[LOG[a] + LOG[b]]

pts = list(range(1, 7))                        # 6 distinct evaluation points in GF(16)*
worst = 6
for m in range(16 ** 4):
    c = [(m >> (4 * i)) & 15 for i in range(4)]
    if not any(c):
        continue
    wt = 0
    for p in pts:                              # evaluate the degree-<4 polynomial at p
        v, pk = 0, 1
        for ci in c:
            v ^= mul(ci, pk)
            pk = mul(pk, p)
        wt += (v != 0)
    worst = min(worst, wt)
assert worst >= 3, worst
print("RS[6,4,3]_16 verified: 65536 codewords, minimum weight %d >= 3" % worst)
print("  => 65536 six-subsets of a 96-set pairwise meeting in <= 6 - 3 = 3 = w/2\n")

# --------------------------------------------------------------------- the inputs
SIGN0 = {93: 1 << 31, 94: 1 << 31, 95: 1 << 32, 96: 1 << 33}     # codetables.de 2026-08-26
GRID = {(96, 24): 4 ** 9, (96, 6): 16 ** 4,
        (95, 23): 4 ** 8, (94, 23): 4 ** 8, (93, 23): 4 ** 8}

# Acw returns max(GRID, derived), so an entry LARGER than the derivation would raise the
# level-1 charge -- and the dimension-96 total with it -- on this table's authority alone.
# Every entry currently equals what ers_exposure.acw derives, and this says so out loud:
# the table is a cache of the derivation, never a source of its own.  Should the derivation
# ever improve, the assertion is what tells you to delete the stale row rather than let
# max() quietly keep it.
for _nw, _v in sorted(GRID.items()):
    _derived = _acw(*_nw)[0]
    assert _v <= _derived, (
        "GRID%s = %d EXCEEDS the derived %d: a hand-written bound would enter the claim "
        "without a construction behind it" % (_nw, _v, _derived))

def Acw(n, w):
    """A LOWER bound on A(n,w,w): supports of weight w pairwise meeting in <= w/2."""
    if w == 1:
        return n                                # distinct singletons meet in 0
    if w == 2:
        return comb(n, 2)                       # distinct pairs meet in <= 1
    if w == 3:                                  # partial Steiner triple system; Hanani's
        return (n // 3) * ((n - 1) // 2)        # maximum is at least this for every n
    # Anything not hand-checked above comes from ers_exposure.acw, which DERIVES the
    # q-ary grid map for every (n, w) rather than listing the few needed here.  This
    # returned 1 until 2026-08-26, which silently zeroed the level-2 charge: at n = 95
    # it missed A(95,5,5) >= 19^3 (95 = 5 x 19, RS[5,3,3] over GF(19)) and A(95,4,4),
    # so the lower levels came to 206 instead of 536366 and the threshold printed here
    # was 265852 against the true 265786 (265798 since the axis layers were rebuilt on
    # 2026-08-27 and dimension 95 rose; the threshold tracks the claim).  An understated
    # rival OVERSTATES the
    # threshold, i.e. it flatters the claim, which is the direction that matters.
    return max(GRID.get((n, w), 1), _acw(n, w)[0])

def Abin(w, d):
    """A(w,d), or the trivial 1 where nothing is tabulated (w > 60)."""
    try:
        return P.Abin(w, d)
    except KeyError:
        return 1

def sd(w):
    return -((-w) // 4)                         # ceil(w/4)

def term(n, w):
    return SIGN0[n] if w == n else Acw(n, w) * Abin(w, sd(w))

def best(n):
    b, bc = -1, None
    for ch in P.chains(n):
        v = sum(term(n, w) for w in ch)
        if v > b:
            b, bc = v, ch
    return b, bc

# --------------------------------------------------------------------- 1. what ERS reaches
print("Lower bounds on the Edel-Rains-Sloane count, from constructions only:\n")
print("  n    best chain            N(n) >=           vs Gamma_72   vs this repo's claim")
RES = {}
CHAIN = {}
for n in (93, 94, 95, 96):
    v, ch = best(n)
    RES[n] = (v, ch)
    cl = CLAIMED.get(n)
    print("  %2d   %-20s  %-16d  x%-11.3f %s"
          % (n, str(ch), v, v / GAMMA72,
             ("x%.3f" % (v / cl)) if cl else "-- (nothing claimed at 96)"))

v96, ch96 = RES[96]
print("\n  n = 96, chain %s, term by term:" % (ch96,))
for w in ch96:
    print("     w = %2d   A(96,%d,%d) = %-8d x A(%d,%d) = %-11d  ->  %13d"
          % (w, w, w, 1 if w == 96 else Acw(96, w), w, sd(w),
             SIGN0[96] if w == 96 else Abin(w, sd(w)), term(96, w)))
print("     %sTOTAL  tau(96) >= %d" % (" " * 50, v96))
print("\n  A LOWER bound on ERS, ABOVE Gamma_72, so ERS wins dimension 96 outright:")
print("     level 0 alone, the [96,33,24] code                 %d" % SIGN0[96])
print("     the Gamma_72 cap construction at k = 24            6480558568")
print("     ERS chain (96,24,6,1)                              %d" % v96)

# --------------------------------------------------------------------- 2. the exposure
print("\n" + "=" * 92)
print("EXPOSURE: what would overturn dimensions 93, 94, 95")
print("=" * 92)
print("""
The three numbers above for n = 93..95 are lower bounds on ERS, so they do NOT prove ERS
loses there.  What follows is the threshold instead: hold the chain at (n, 23, w2, 1) with w2 chosen, take
the level-0 and lower levels as given, and ask how large the level-1 support code A(n,23,23)
would have to be for ERS to reach this repository's claim.  A(23,6) = %d, so each unit of
A(n,23,23) is worth that much.
""" % Abin(23, 6))
print("  n   claim           ERS level 0     lower levels   A(n,23,23) needed   have    Johnson <=")
for n in (93, 94, 95):
    # the level-2 weight is CHOSEN, not fixed: 16*w2 <= 4*23 admits w2 <= 5, and the
    # strongest rival gives the tightest threshold.  Fixing it at 5 understated the
    # rival, because w2 = 4 charges more (A(n,4,4) is a packing, A(n,5,5) a grid map).
    w2 = max((w for w in range(1, 6) if 16 * w <= 4 * 23), key=lambda w: term(n, w))
    CHAIN[n] = (n, 23, w2, 1)
    lower = term(n, w2) + term(n, 1)
    need = -((-(CLAIMED[n] - SIGN0[n] - lower)) // Abin(23, 6))   # ceiling division
    have = Acw(n, 23)
    jb = P.johnson_cw(n, 24, 23)
    print("  %2d  %-14d  %-14d  %-13d  %-18d  %-6d  %d"
          % (n, CLAIMED[n], SIGN0[n], lower, need, have, jb))
    RES[n] = RES[n] + (need, have, jb)

# ------------------------------------------------------------------ cross-check, two ways
# Two independent paths to the same three thresholds: ers_exposure builds the whole ranking
# from RESULTS.md and derives every constant-weight bound, while this file costs one chain
# per dimension from codetables values entered by hand.  They disagreed from the day both
# existed until 2026-08-26 and nothing caught it, because the only check between them ran
# ONE WAY: ers_exposure asserts its bounds are at least the ones hardcoded here, which stays
# true precisely when this file understates.  Assert equality instead -- a one-sided check
# between two implementations is blind to the side it does not test.
from ers_exposure import ranked as _ranked
_rk = {r[0]: r for r in _ranked()}
for _n in (93, 94, 95):
    _here, _there = RES[_n][2], _rk[_n][6]
    assert _here == _there, (
        'threshold at dimension %d: ers96 says %d on chain %s, ers_exposure says %d on %s'
        % (_n, _here, CHAIN[_n], _there, _rk[_n][3]))
print('cross-check: ers_exposure.ranked() agrees on all three thresholds '
      '(%d / %d / %d)' % tuple(RES[k][2] for k in (93, 94, 95)))

n = 95
need, have, jb = RES[95][2:]
print("""
So dimension %d is the exposed one: it needs A(95,23,23) to stay under %s, the grid map
gives %s, and the Johnson bound allows up to %s.  The margin is a factor %.2f on the best
construction in hand.  For calibration, ERS's own A(64,16,16) = 30828 beats the same grid map
(4^7 = 16384 at n = 64, since codetables' best quaternary length-16 code of distance 8 is
[16,7,8]) by a factor %.2f, so a factor %.2f is not obviously out of reach and
this is a genuine risk to dimension %d rather than a settled question.  What protects the
claim is that the value to beat throughout 73-95 is the DIRECT SUM, not an ERS total.
ERS did publish inside this range -- dimension 80, at 1 368 532 064 -- and the direct sum
6 218 175 840 beats it by a factor 4.5.  The exposure is to a referee with a better
constant-weight code, not to the literature as it stands.

None of it touches dimension 96, where the inequality points the other way.
""" % (n, "{:,}".format(need).replace(",", " "), "{:,}".format(have).replace(",", " "),
       "{:,}".format(jb).replace(",", " "), need / have, 30828 / 4 ** 7, need / have, n))
