#!/usr/bin/env python3
"""The Edel-Rains-Sloane LEVEL-0 floor, swept over every dimension this repository claims.

The dimension-96 work asked "how big is A(n_0, ceil(n_0/4)) at n_0 = 96?" and found 2^33
where 2^32 had been assumed.  The same question has an answer in EVERY dimension, and this
script asks it everywhere rather than at one point.

    tau(n) >= max over n_0 <= n of A(n_0, ceil(n_0/4))

is already a complete kissing configuration on its own: take the n_0 coordinates, put
(+-1/sqrt(n_0)) on them and 0 elsewhere, and two sign vectors at Hamming distance d meet at
1 - 2d/n_0 <= 1/2 exactly when d >= n_0/4.  It is level 0 of Edel-Rains-Sloane with a single
support, it needs nothing else, and it is a one-line consequence of any published code.

    python ers_sweep.py

The k values below are the largest with minimum distance >= ceil(m/4) in Grassl's table of
best known LINEAR codes, https://codetables.de, read 2026-08-26.  Within a group of four the
value is carried down by SHORTENING (a [m,k,d] code gives [m-1,k-1,>=d]), which is why only
the group tops have to be looked up; where the table was read directly at m < 4d it agreed.
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DELIV = os.path.join(HERE, '..', '..', '..')
sys.path.insert(0, os.path.join(DELIV, 'common'))
sys.path.insert(0, os.path.join(HERE, '..', 'dim32-44-ers-audit'))
import PIPELINE as P

# ---- k for the best known [m, k, ceil(m/4)] binary linear code, m = 4d ------------------
# codetables.de, 2026-08-26.  Read directly; the whole finding rests on these twelve numbers.
TOP_K = {48: 24, 52: 21, 56: 23, 60: 25, 64: 28, 68: 25, 72: 28,
         80: 30, 84: 28, 88: 30, 92: 30, 96: 33}
DIRECT = {                       # rows read directly rather than shortened, as a cross-check
    61: 25, 62: 26, 63: 27,      # [61,25,16] [62,26,16] [63,27,16]
    69: 25, 70: 26, 71: 27,      # [69,25,18] [70,26,18] [71,27,18]
    93: 31, 94: 31, 95: 32,      # [93,31,24] [94,31,24] [95,32,24]
}

def sign_code(m):
    """lower bound on A(m, ceil(m/4)) from the best known linear code, or from shortening

    TOP_K and DIRECT are LINEAR codes and cover 48..96 only, so both lookups miss below 48 and
    this returned 1 -- a silent identity, and the third of that shape found in this directory
    on 2026-08-26 (Acw's GRID.get((n,w), 1) and ers_exposure.sign0 were the others).

    IT IS LATENT HERE, and the measurement matters more than the fix: `sign_code` is wrong at
    every length 32..47, but `level0` only ever reports for n >= 49 (see the loops below), and
    for those the maximum is attained at some n_0 >= 48 where the tables do cover.  Checked:
    ZERO of the 43 swept dimensions change -- 49..63, 68..71 and 73..96, the claims at or
    above 49, which is all this script reports.  Were the sweep extended below 49 it would bite
    once -- level0(38) and level0(39) are 180 224 and 327 680, and read 1 without this.

    Fixed anyway, because a lookup that answers 1 when it means "I have no entry" is the shape
    that cost this directory two real defects the same day, and the next person to widen the
    range would inherit it.  P.Abin is the repository's binary table over n = 0..512 and admits
    NONLINEAR codes, which is legitimate here: level zero needs only a binary code of length m
    with d >= ceil(m/4), linear or not.
    """
    best = DIRECT.get(m, 0)
    d = -((-m) // 4)
    top = 4 * d                                  # the group top has the same required d
    if top in TOP_K:
        best = max(best, TOP_K[top] - (top - m))
    linear = 1 << best if best > 0 else 1
    try:
        return max(linear, P.Abin(m, d))
    except KeyError:                             # a table miss, and nothing else
        return linear

def level0(n):
    """max over n_0 <= n of A(n_0, ceil(n_0/4)); the whole configuration, on its own"""
    b, bn = 1, None
    for n0 in range(32, n + 1):
        v = sign_code(n0)
        if v > b:
            b, bn = v, n0
    return b, bn

# ---- this repository's claims, read from the generated RESULTS.md -----------------------
CLAIM = {}
for line in open(os.path.join(DELIV, 'RESULTS.md'), encoding='utf-8'):
    m = re.match(r'\| (\d+) \| ([\d\s  ]+) \| \*\*([\d\s  ]+)\*\* \|', line)
    if m:
        CLAIM[int(m.group(1))] = int(re.sub(r'\D', '', m.group(3)))

print(__doc__)
print("  dim   this repo claims   ERS level 0 alone   n_0   margin    verdict")
broken, tight = [], []
for n in sorted(d for d in CLAIM if d >= 49):
    v, n0 = level0(n)
    c = CLAIM[n]
    if v > c:
        verd, mark = "**BEATEN BY LEVEL 0 ALONE**", broken
    elif v == c:
        # the claim IS this construction: dimensions 62 and 63, after the cap construction
        # was shown unable to reach them with any class family at all
        verd, mark = "this IS the claim (dim62-63-ers-chain)", None
    elif c < 2 * v:
        verd, mark = "under a factor 2 of margin", tight
    else:
        verd, mark = "clear", None
    if mark is not None:
        mark.append(n)
    print("  %3d   %-16d   %-17d   %3s   x%-7.2f  %s"
          % (n, c, v, n0 if n0 else "-", c / v, verd))

print()
print("  BEATEN: %s" % (broken if broken else "none"))
print("  within a factor 2: %s" % (tight if tight else "none"))

if broken:
    print("""
  These are not "unpublished evaluations that a referee might do".  A(62,16) >= 2^26 and
  A(63,16) >= 2^27 are the codes [62,26,16] and [63,27,16], which are in the table, and
  ceil(62/4) = ceil(63/4) = 16, so the sign vectors of those codes ARE kissing configurations
  in dimensions 62 and 63 with nothing added.  By the standard this repository already
  applies in dimensions 32-44 -- where an evaluation of Edel-Rains-Sloane with the current
  tables is treated as the live record, and is what shows Cohn's table stale at 32, 33, 34
  and 37 -- these two claims are not improvements.

  published.py's dimension 49-63 note says the Edel-Rains-Sloane construction "does NOT reach
  52416000 anywhere in this range" because "the value 2^28 = A(64,16) that carries dimension
  64 is not available below n_0 = 64".  2^28 is indeed not available below 64.  But 2^27 is
  available at n_0 = 63 and 2^26 at n_0 = 62, and both are far above 52 416 000, so the
  conclusion does not follow from the premise.""")


# ---- level 1: for the dimensions level 0 does not beat, how big a target is it? ---------
# NOTE ON WHAT IS *NOT* HERE.  The clean thing would be an UPPER bound on the ERS count, and
# there is one: A(n,w,w) <= Johnson's nested-floor bound is a theorem.  It was computed and
# it is useless -- at n = 63 the chain (32,8,2) with Johnson support codes allows 1.3e11,
# 2000x the claim, because Johnson is astronomically above anything constructible (it gives
# 3.2e5 at n = 64, w = 16 where ERS's actual code has 30828).  Declaring every dimension "not
# safe" on that basis would be noise, not a risk assessment.  So no upper bound is offered.
# What follows is the size of the TARGET instead: how large the level-1 support code would
# have to be, next to the largest one anybody has published at a comparable length.
print("""
  Level 0 is only the first term.  For a chain topped at n_0 the next level has weight
  w <= n_0/4 and contributes A(n,w,w) * A(w, ceil(w/4)).  For every dimension level 0 does
  NOT already beat, here is the best available sign code at that level and how big the
  support code would have to be to close the rest of the gap.  For scale, the two largest
  such codes anyone has published are ERS's own A(64,16,16) = 30828 and A(80,16,16) = 143780.
""")
print("  dim   gap to close      w    A(w,ceil(w/4))   A(n,w,w) needed   vs 143780")
for n in sorted(d for d in CLAIM if d >= 49):
    v, n0 = level0(n)
    if v >= CLAIM[n] or n0 is None:
        continue
    gap = CLAIM[n] - v
    bw, bs = None, 0
    for w in range(2, n0 // 4 + 1):
        try:
            sc = P.Abin(w, -((-w) // 4))
        except KeyError:
            continue
        if sc > bs:
            bw, bs = w, sc
    need = -((-gap) // bs)
    flag = "  <-- inside a factor 2 at level 0" if n in (60, 61, 68, 95) else ""
    print("  %3d   %-15d  %3d   %-14d   %-15d   x%.2f%s"
          % (n, gap, bw, bs, need, need / 143780.0, flag))
print("""
  Everything from 49 to 59 needs a length-49-to-59 weight-12 code of 139000-155000 words,
  roughly the published length-80 value at less than three quarters of the length, and 73-95
  need millions.  The four flagged rows are the ones worth watching: 60 and 61 need about
  24000, 68 needs 45000 and 95 needs 257000 -- the first three of those are the same order as
  A(64,16,16) = 30828, which is why they are called exposed rather than safe.

  None of this touches 62 and 63.  Those need no level 1 at all, and as of 2026-08-23 they
  ARE the level-0 construction: verifications/improved/dim62-63-ers-chain/.  The cap
  construction over P_48 cannot reach them with any class family -- its ceiling, every
  class at the maximum 7069, is 66075240 and 70543480 -- so they were moved rather than
  patched; see ../../improved/dim49-63-p48-caps/scripts/ceiling.py.
""")
