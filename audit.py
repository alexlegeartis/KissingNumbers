#!/usr/bin/env python3
"""Consistency audit of every claim in this repository.

    python audit.py

This does NOT re-prove anything -- each package's own script does that. What it does is
check that the claims are mutually consistent, agree with the published state of the art,
and respect the elementary facts that tau is non-decreasing and bounded above by the records
in higher dimensions. A construction can be individually correct and collectively wrong; this
catches the second kind of error, which is the kind that survives a careful reading.

Where a claim is cheap to recompute, it is recomputed here rather than quoted:

  * dimension 71 comes from theta.py, exactly, in a fraction of a second;
  * dimensions 49-63, 73-80 and 81-95 are recomputed from the shipped class-size tables;
  * the rest are quoted, each with the script that establishes it.

CHECKS
  1. every claim strictly exceeds the best previously published value in its dimension,
     with that value taken from common/published.py and never typed by hand;
  2. every claim in a dimension Cohn's table covers cites exactly Cohn's figure;
  3. the claims are monotone in the dimension, and no claim exceeds the best value known in
     ANY higher dimension -- including this repository's own claims;
  4. no claim is made in a dimension that is attained by the published construction, withdrawn,
     short of its floor, or already in the literature.  The four sets are ATTAINED, WITHDRAWN,
     SHORT and RECOVERED below, and the run prints their contents rather than this docstring
     repeating them -- which is how "69" stayed here after 69 left SHORT and became a claim;
  4b. every claim exceeds what the FULL Edel-Rains-Sloane chain reaches from published codes,
     via ers_exposure.best_chain, IMPORTED rather than copied.  Level 0 alone is not enough:
     dimensions 62 and 63 were once claimed AT level 0 while the chain (n,15,2) reached more,
     and the sweep knew while this audit did not.  Every chain input is a construction, so the
     floor is a LOWER bound on ERS -- clearing it is necessary, not sufficient, and
     verifications/closed/dim96-ers-takeover/ers_exposure.py reports the remaining room.
     Four claims EQUAL the floor rather than exceed it, and do so by being the chain
     themselves: 39, 62, 63 and 96.
     THIS CHECK WAS HOLLOW BELOW DIMENSION 48 UNTIL 2026-08-26.  ers_exposure.sign0 read
     only tables of LINEAR codes covering 48..96 and returned 1 for anything else, which
     dropped the chain's largest term entirely; best_chain then returned 0.78, 0.81, 0.48
     and 0.29 of the totals ERS themselves published at 32, 36, 40 and 44.  The repository's
     binary code table was never the weak part -- sign0 simply was not consulting it.  It
     now falls back to that table, which also admits NONLINEAR codes, and best_chain returns
     1.255, 1.104, 1.073 and 1.000 at those four dimensions: at or above published ERS,
     which is what an evaluator with better tables must give.  Nothing in 49..96 moved, so
     the paper's section 8.2 margins are unaffected.  closed/dim32-44-ers-audit/ remains the
     dedicated low-range evaluator, with Brouwer's tables including Echols arXiv:2608.13906;
     it is what shows the chain reaching 570 236 at dimension 38 against a claim of 591 612,
     and 756 116 at dimension 39, which IS the claim there;
  4c. no claim exceeds a PROVEN upper bound.  Steps 1-4b ask whether a claim is an
     improvement; this asks whether it is possible at all, against the upper-bound column of
     Cohn's table -- the only place in this repository that column is read;
  5. every package directory named by a claim exists and contains a README;
  1b. no floor falls below an Edel-Rains-Sloane total published at or below that dimension.
     K is monotone, so a published value above our floor would mean the floor is wrong --
     dimension 80 is the one that bites, and three files here once denied it exists;
  4a. where two constructions reach the same dimension, the better one is the one claimed.
     Dimensions 62 and 63 are why: the cap construction reaches them and the ERS chain
     reaches further, and whichever call came second would otherwise have won;

  Then the documentation, checked the way the claims are.  Prose is where every defect of
  the last two days lived, and none of it moved a digit:

  5a. every script in every package is named in that package's README;
  5b. every package that ships a SHA-256 manifest still matches it;
  5c. the working directory's README accounts for every path beside it;
  5d. verifications/README.md counts the tiers and packages it actually has;
  5e. every path a .gitignore rule names still exists;
  5f. every pointer into an outside document says which document it numbers -- a preprint
     and its journal version need not number their equations alike;
  5g. the README's description of the paper still partitions the claims;
  5h. the root README names every script that reads Cohn's coordinate files;
  5i. update.py's improvement guide lists the right dimensions per status;
  5j. this docstring lists every check this script runs;
  5o. every rotated axis layer lies in the axis space R^k, not in the coordinates
      the direction set happens to be written in
  5p. the layer floor in verify_poles.py's own header matches the shipped layers
  5r. the record axis layers are sealed and match the README table describing them
  5q. the axis-layer figures that live only in prose -- the denominator digit-range
      and what the shell rule cost -- match the shipped layers
  5n. no .py in the tree carries a tau table disagreeing with common/published.py
  5s. the run instructions quote the number of scripts run_all.py has, and the
      root README the number of files and the size the published tree has.  Both
      went stale by a fifth without failing anything: nothing measured them;
  5m. no markdown table outside KNOWLEDGE.md quotes a superseded value of a claim -- the
      check is relational, because the separator is sometimes U+202F and a grep for the
      old digits misses it;
  5l. no LaTeX escape in a .md or .tex has been halved into a control character -- writing
      "\tau" through a non-raw Python string leaves a TAB, and three files had it;
  5k. .gitattributes states the line-ending mix the tree actually has, `* -text` is
     still there, and no source or documentation file mixes endings within itself.  The
     scan opens ten extensions and that is the scope of the claim -- one archived HTML
     copy of Cohn's page does mix, deliberately, and .gitattributes says so.  The mix is harmless --
     `* -text` preserves bytes, which is what the SHA-256 manifest needs -- but a
     maintainer who thought the tree uniform might remove `* -text` as redundant;
  6. every claim says what would have to change for it to improve, and the claim COUNT
     agrees everywhere it is stated.

Nothing here needs anything outside this directory.  Where a check is about the author's
working folder -- the root README, the paper's PDF -- it says so and stands down in a clone.
"""
import fnmatch
import io
import hashlib, os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'common'))
from published import (COHN, floor_for, upper_for,       # noqa: E402
                       ERS_PUBLISHED)
from status import STATUS, LEGEND, status_for, runs             # noqa: E402
from published import TAU as PUBTAU
from theta import onepoint                                # noqa: E402

IMP = os.path.join(HERE, 'verifications', 'improved')

fail = 0


def bad(msg):
    global fail
    fail += 1
    print("   FAIL: " + msg)


# --------------------------------------------------------------------------- the claims
# (dimension, value, package, how it is established here)
claims = []


SUPERSEDED = []


def add(d, v, pkg, how):
    """record a claim, keeping the BEST construction when two of them reach the same dimension

    Dimensions 62 and 63 are the reason this is not just an append: the P_48 cap construction
    reaches them, but a single Edel-Rains-Sloane sign code reaches further, and whichever call
    happened to come second would otherwise have won.  A displaced claim is remembered so the
    report can say what it was."""
    for t, (d0, v0, pkg0, how0) in enumerate(claims):
        if d0 != d:
            continue
        if v > v0:
            SUPERSEDED.append((d, v0, pkg0, v, pkg))
            claims[t] = (d, v, pkg, how)
        else:
            SUPERSEDED.append((d, v, pkg, v0, pkg0))
        return
    claims.append((d, v, pkg, how))


# --- recomputed exactly, right now, from Venkov's theorem
sol48, ok48, _ = onepoint(48, 6, 52416000, 11)
sol72, ok72, _ = onepoint(72, 8, 6218175600, 11)
if not (ok48 and ok72):
    bad("the one-point surplus equations do not hold")
# P_48's n[0] is NOT a claim of this project: 23766960 is due to Boyvalenkov-
# Cherkashin, Results Math. 80 (2025), no. 1, Paper No. 3, and is equation (3) of the
# preprint arXiv:2312.05121, whose same paragraph observes that those vectors form a
# 47-dimensional kissing configuration.  The equation number belongs to the preprint:
# it is the version this repository read, and two versions need not number alike.
# We recompute it because it is the input to the two-point system and because agreeing
# with a published value to the digit is the strongest check available here.
RECOVERED = {47: 23766960, 46: 12309600}
if int(sol48[0]) != RECOVERED[47]:
    bad("P_48 one-point n[0] = %d, expected the published 23766960" % int(sol48[0]))
add(71, int(sol72[0]), 'dim70-71-gamma72-cross-sections', 'recomputed by theta.py')

# --- recomputed from the shipped class-size tables
TAU = {k: COHN[k] for k in range(1, 24)}
# tau_lat(k), the best known LATTICE kissing number
TAULAT = {1: 2, 2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240, 9: 272, 10: 336,
          11: 438, 12: 756, 13: 918, 14: 1422, 15: 2340, 16: 4320, 17: 5346, 18: 7398}
ANTIP = TAULAT
# lambda(k): the dimension 49-63 scheme sits at t = 3/4, where a part is an ANTIPODAL PAIR, so
# what it needs is the largest antipodally symmetric 60-degree code -- not tau_lat(k).  Cohn's
# record configurations are antipodally symmetric for k = 9, 10, 11, 13, 14, 15 (verified by
# dim49-63-p48-caps/scripts/verify_capdirs.py), giving lambda(k) = tau(k)/2 there.  k = 12 is
# excluded because tau(12) = 841 is odd.  For k <= 8 the records are the root systems.
SYMK = {9, 10, 11, 13, 14, 15}
LAMBDA = {k: (TAU[k] // 2 if k in SYMK else TAULAT[k] // 2) for k in range(1, 16)}
# root-system triple/pair partitions verified by dim73-95/scripts/triples.py
ROOTPARTS = {1: (0, 1), 2: (2, 0), 3: (4, 0), 4: (8, 0), 5: (12, 2), 6: (24, 0),
             7: (42, 0), 8: (80, 0)}

# ONE family of 1400 pairwise-disjoint P_48 classes, regenerated and verified completely by
# dim49-63-p48-caps/scripts/regenerate.py from three small files and a seed.  It replaced a
# 1170-class family plus a 130-class extension whose base no shipped script reproduces; running
# the same construction out to 1400, with the first 400 images chosen to overlap least, beats
# that pair on all fifteen dimensions.  lambda(15) = 1282, so the family must reach 1282.
S48 = np.sort(np.load(os.path.join(
    IMP, 'dim49-63-p48-caps', 'data', 'class_sizes_p48.npy')).ravel())[::-1].astype(np.int64)
c48 = np.concatenate([[0], np.cumsum(S48)])
for k in range(1, 16):
    lam = LAMBDA[k]
    if lam > len(S48):
        bad('k = %d needs %d classes, the family has %d' % (k, lam, len(S48)))
        continue
    Sk = int(c48[lam])
    add(48 + k, 52416000 + 2 * Sk + TAU[k], 'dim49-63-p48-caps',
        'lambda(%d) = %d, S = %d from class_sizes_p48.npy' % (k, lam, Sk))

Z = np.load(os.path.join(IMP, 'dim73-95-gamma72-caps', 'data', 'disjoint_classes.npz'))
S72 = np.sort(np.array([len(Z[k]) for k in Z.files], dtype=np.int64))[::-1]
c72 = np.concatenate([[0], np.cumsum(S72)])
for k in range(1, 9):
    T, P = ROOTPARTS[k]
    gain = 4 * int(c72[T]) + 2 * int(c72[T + P] - c72[T])
    poles = 2 if k == 1 else 0            # k = 1 sits at t = 3/4, where poles are automatic
    # the axis layer at t = 2/3, rotated where that beats the lattice shell (k >= 5)
    _rf = os.path.join(IMP, 'dim73-95-gamma72-caps', 'data', 'poles_rot_k%d.npz' % k)
    _pf = os.path.join(IMP, 'dim73-95-gamma72-caps', 'data', 'poles_k%d.npy' % k)
    if os.path.exists(_rf):
        poles += int(len(np.load(_rf)['keep']))
    elif os.path.exists(_pf):
        poles += 2 * len(np.load(_pf))
    add(72 + k, 6218175600 + gain + poles, 'dim73-95-gamma72-caps',
        'recomputed from disjoint_classes.npz')

# TWO class families are shipped and each dimension takes whichever gives it more, by the
# rule stated in final81-95.py; they are never mixed, since disjointness ACROSS families is
# not claimed.  The audit has to apply the same rule -- reading only the first of them made
# it recompute dimensions 86-94 several million low while RESULTS.md and the paper carried
# the right values, and audit.py is what writes RESULTS.md.
DATA = os.path.join(IMP, 'dim73-95-gamma72-caps', 'data')
FAMS = []
for _fn in ('class_sizes_32000.npy', 'class_sizes_gpu.npy'):
    _fp = os.path.join(DATA, _fn)
    if os.path.exists(_fp):
        _s = np.sort(np.load(_fp).ravel())[::-1].astype(np.int64)
        FAMS.append((_fn, _s, np.concatenate([[0], np.cumsum(_s)])))
if not FAMS:
    bad('no Gamma_72 class-size table found')
D80 = [c for c in claims if c[0] == 80][0][1]
# k = 9..23 read their zero-sum triple partitions from the shipped CERTIFICATES, choosing per
# k whichever of the best lattice and the record kissing configuration gives more units; both
# are re-checked exactly by scripts/verify_parts.py.
run = D80
for k in range(9, 24):
    best = None
    cands = [('lattice', 'lambda9_W.npy' if k == 9 else 'lam%d_W.npy' % k, 'parts_%d.npz' % k),
             ('record', 'rec_%d_W.npy' % k, 'rec_%d_parts.npz' % k),
             # since 2026-08-30: the Kappa sections K_12 and K_13, which are larger than
             # Lambda_12 and Lambda_13 and split perfectly into zero-sum triples
             ('kappa', 'kap%d_W.npy' % k, 'kap%d_parts.npz' % k)]
    if k == 11:
        cands.append(('record(Z[sqrt2])', None, 'rec11_parts.npz'))
    for src, wf, pf in cands:
        f = os.path.join(DATA, pf)
        if not os.path.exists(f):
            continue
        z = np.load(f)
        T, P = len(z['triples']), len(z['pairs'])
        if best is None or 2 * T + P > 2 * best[0] + best[1]:
            best = (T, P, src)
    if best is None:
        bad("no cap-direction certificate for k = %d" % k)
        continue
    T, P, src = best
    cands = [4 * int(_c[T]) + 2 * int(_c[T + P] - _c[T])
             for _nm, _s, _c in FAMS if len(_s) >= T + P]
    if not cands:
        bad('no class family has %d classes, needed for k = %d' % (T + P, k))
        continue
    gain = max(cands)
    # The axis layer, in whichever form the package ships for this k.  FOUR are possible and
    # the largest is the one used, which is how final81-95.py chooses:
    #
    #   poles_rec_<k>.npz    a rotated copy of the RECORD configuration of R^k, which is not
    #                        the cap set and is not bounded by |Z| (k = 19, 20, 21);
    #   poles_sqrt2_<k>.npz  the same idea over Z[sqrt2] (k = 11);
    #   poles_rot_<k>.npz    a rotated copy of the cap directions;
    #   poles_<k>.npy        the older lattice-shell rule, antipodal pairs.
    #
    # Reading only the rotated file recomputed dimensions 83, 91, 92 and 93 at their previous
    # values -- and then agreed with a RESULTS.md that was equally stale, which is exactly
    # what this recomputation exists to prevent.  scripts/verify_record.py, scripts/verify11.py
    # and scripts/verify_poles.py check the four kinds exactly.
    # A LAYER BELONGS TO A CAP SET.  Every kind below is indices into W, or a rotation of
    # W, so a layer built on Lambda_k says nothing about K_k -- taking the max over all of
    # them regardless of src would let a retired layer of one configuration be counted for
    # another.  Where the chosen source is the Kappa section, only its own layer counts.
    _lay = []
    if src == 'kappa':
        _names = ('poles_kap_%d.npz' % k, 'poles_kaprec_%d.npz' % k)
    else:
        _names = ('poles_rec_%d.npz' % k, 'poles_sqrt2_%d.npz' % k, 'poles_rot_%d.npz' % k)
    for _nm in _names:
        _p = os.path.join(DATA, _nm)
        if os.path.exists(_p):
            _lay.append(int(len(np.load(_p)['keep'])))
    pf = os.path.join(DATA, 'poles_%d.npy' % k)
    if os.path.exists(pf) and src != 'kappa':
        _lay.append(2 * len(np.load(pf)))
    if _lay:
        gain += max(_lay)
    val = max(6218175600 + gain, run)
    run = val
    add(72 + k, val, 'dim73-95-gamma72-caps',
        'recomputed from the %s partition at k = %d and class_sizes_32000.npy' % (src, k))

# --- quoted, each verified by the named script
add(25, 197569, 'dim25-lens-heads', 'verify.py (exact arithmetic; fullcheck.py is the independent net)')
add(26, 199632, 'dim26-27-iota-triangles', 'verify26.py (exact: integers, and sympy for the hexagon)')
add(27, 201010, 'dim26-27-iota-triangles', 'verify27.py (exact: integers, and sympy for the rotated cuboctahedron)')
add(38, 591612, 'dim38-leech-large-codimension', 'scripts/verify.py')
add(39, 756116, 'dim39-ers-constant-weight', 'scripts/verify.py')
# The Edel-Rains-Sloane chain (n, 15, 2).  Level 0 is a single sign code: [62,26,16] and
# [63,27,16] (codetables.de, read 2026-08-23) have d = 16 = ceil(62/4) = ceil(63/4), so their
# 2^26 and 2^27 sign vectors are kissing configurations outright.  Levels 1 and 2 are the rest
# of the construction and are BUILT, not cited: a cyclic [15,6,8]_4 turns into 4096 supports of
# weight 15 by the grid map, each carrying the 1024 sign patterns of a [15,10,4], and level 2 is
# every pair of coordinates with all four signs.  Both totals exceed the CEILING of the P_48 cap
# construction in those dimensions -- 66075240 and 70543480, every class at the maximum 7069 --
# so the cap route cannot reach them with any class family at all.
add(62, (1 << 26) + 4096 * 1024 + 1891 * 4, 'dim62-63-ers-chain',
    'scripts/verify_chain.py (chain (62,15,2): 2^26 + 4096*1024 + 1891*4)')
add(63, (1 << 27) + 4096 * 1024 + 1953 * 4, 'dim62-63-ers-chain',
    'scripts/verify_chain.py (chain (63,15,2): 2^27 + 4096*1024 + 1953*4)')
add(70, 1249778250, 'dim70-71-gamma72-cross-sections', 'derive.py (exact dual certificate)')
add(69, 627822180, 'dim68-69-gamma72-cross-sections',
    'verify69.py (A_3 at the Hermite floor, det 256, exhibited inside the D_4 four-tuple)')
add(68, 361275480, 'dim68-69-gamma72-cross-sections',
    'verify68.py (the exact optimum of the programme, D_4 Gram at the Hermite floor)')
# Dimension 96, the Edel-Rains-Sloane chain (96,24,6,1), on exactly the footing of 62 and 63:
# their construction evaluated at a length they never evaluated.  Their section 3 gives seven
# worked examples -- 32, 36, 40, 44, 64, 80, 128 -- and 96 is not among them (arXiv:math/
# 0207291, read 2026-08-26), while Cohn's table stops at 48 and 72 and the Nebe-Sloane table it
# replaced listed 1-40, 42, 44, 48, 64, 72, 80, 128.  Substituting current code tables into
# their formula is sanctioned by their own Remark (1).  The four levels are
#     A(96,24) >= 2^33          [96,33,24] on codetables.de       8 589 934 592
#     A(96,24,24) >= 4^9        grid map, [24,9,12]_4             x A(24,6) = 16384
#     A(96,6,6)  >= 16^4        grid map, RS[6,4,3]_16, BUILT     x A(6,2)  = 32
#     A(96,1,1)  =  96          singletons                        x A(1,1)  = 2
# This was the FLOOR here until 2026-08-26, which hid it: a baseline set to one's own new
# result turns the claim into a non-result.  See common/published.py's floor_for docstring.
add(96, (1 << 33) + 4 ** 9 * 16384 + 16 ** 4 * 32 + 96 * 2, 'dim96-ers-takeover',
    'ers96.py (chain (96,24,6,1): 2^33 + 4^9*16384 + 16^4*32 + 96*2)')

claims.sort()
print(__doc__)
print("=" * 92)
print("%d claims, dimensions %d..%d" % (len(claims), claims[0][0], claims[-1][0]))
print("=" * 92)
print("  dim   previously published        this work            factor   package")
for d, v, pkg, how in claims:
    prev, _ = floor_for(d)
    print("  %3d   %-24s    %-20s %6.2f   %s" % (d, "{:,}".format(prev),
                                                 "{:,}".format(v), v / prev, pkg))

print()
print("1. every claim strictly exceeds the best previously published value")
for d, v, pkg, how in claims:
    prev, why = floor_for(d)
    if v <= prev:
        bad("dim %d: %d does not exceed %d (%s)" % (d, v, prev, why))
print("   %d claims checked against common/published.py" % len(claims))

print()
print("1b. no floor falls below a published Edel-Rains-Sloane total")
# K is monotone, so a floor in dimension d must be at least every ERS total published at or
# below d.  This is the check that was missing when three files claimed nobody had published
# ERS in 73-95: dimension 80 is in that range.  A published value ABOVE our floor would mean
# the improvement is measured against the wrong baseline and overstated.
_ers_worst = 0
for _d, _v, _pkg, _how in claims:
    _prev, _ = floor_for(_d)
    _pub = [(_m, _w) for _m, _w in ERS_PUBLISHED.items() if _m <= _d]
    if not _pub:
        continue
    _m, _w = max(_pub, key=lambda mw: mw[1])
    if _prev < _w:
        bad("dim %d: floor %d is below ERS's published dimension-%d total %d"
            % (_d, _prev, _m, _w))
    else:
        _ers_worst = max(_ers_worst, _w / _prev)
print("   every floor clears it; the tightest is %.3f of its floor" % _ers_worst)

print("2. where Cohn's table has an entry, the cited previous value is exactly Cohn's")
n = sum(1 for d, v, p, h in claims if d in COHN)
for d, v, pkg, how in claims:
    if d in COHN and floor_for(d)[0] != COHN[d]:
        bad("dim %d: floor_for gives %d, Cohn's table says %d"
            % (d, floor_for(d)[0], COHN[d]))
print("   %d of %d claims fall inside Cohn's table (%s), all consistent"
      % (n, len(claims), runs(sorted(COHN))))

print()
print("3. monotonicity, and no claim exceeding a known record above it")
prev_d, prev_v = None, 0
for d, v, pkg, how in claims:
    if v < prev_v:
        bad("dim %d: %d is below dimension %d's %d (tau is non-decreasing)"
            % (d, v, prev_d, prev_v))
    prev_d, prev_v = d, v
mine = {d: v for d, v, p, h in claims}
OTHER = {64: 331737984,
         96: 12886999232,      # Edel-Rains-Sloane, closed/dim96-ers-takeover/
         128: 218044170240}
best_above = {}
for kd in sorted(set(list(COHN) + list(OTHER) + list(mine))):
    best_above[kd] = max(COHN.get(kd, 0), OTHER.get(kd, 0), mine.get(kd, 0))
for d, v, pkg, how in claims:
    for kd, kv in best_above.items():
        if kd > d and kv and v > kv:
            bad("dim %d: %d exceeds the best known %d in dimension %d" % (d, v, kv, kd))
print("   monotone across all %d claims; every claim respects every record above it"
      % len(claims))

print()
print("4. no claim in a dimension that is attained, withdrawn, or short of its floor")
ATTAINED = {28, 29, 30, 31}           # the cap model reproduces these exactly (26 was here
                                      # until dim26-27-iota-triangles beat it, 2026-09-08)
WITHDRAWN = {44, 45}                  # beaten by Sun-Wang, arXiv:2607.20359v3
SHORT = set()                         # 69 was here until its Gram was exhibited (see the
                                      # dim68-69 package); nothing is short of its floor now
for d, v, pkg, how in claims:
    if d in ATTAINED:
        bad("dim %d is already attained by the published construction" % d)
    if d in WITHDRAWN:
        bad("dim %d was withdrawn; see verifications/superseded/" % d)
    if d in SHORT:
        bad("dim %d does not reach its floor; see verifications/closed/" % d)
    if d in RECOVERED:
        bad("dim %d was already in the literature; see verifications/recovered/" % d)
# Printed from the sets above, not typed out beside them: the hand-written form of this
# line is exactly where "69" survived in the docstring after 69 stopped being short.
print("   %s absent (attained); %s absent (withdrawn)"
      % (runs(sorted(ATTAINED)), runs(sorted(WITHDRAWN))))
print("   %s absent (already in the literature -- recovered, not claimed)"
      % runs(sorted(RECOVERED)))
if SHORT:
    print("   %s absent (short of the floor)" % runs(sorted(SHORT)))

print()
print("4a. where two constructions reach the same dimension, the better one is claimed")
if SUPERSEDED:
    for d, vlo, plo, vhi, phi in sorted(SUPERSEDED):
        print("   dim %d: %s reaches %d, superseded by %s at %d (+%d)"
              % (d, plo, vlo, phi, vhi, vhi - vlo))
else:
    print("   no dimension is reached by more than one construction")

print()
print("4b. every claim beats what the FULL Edel-Rains-Sloane chain reaches")
# Level 0, a single sign code, is only the first term.  Gating on it is not enough: dimensions
# 62 and 63 were once claimed AT level 0 while the chain (n,15,2) reached 6% and 3% more.  The
# floor here is ers_exposure.best_chain, IMPORTED rather than copied -- a second copy of the
# code tables is what had drifted elsewhere, and a floor disagreeing with the sweep would be
# worse than none.  DIRECTION: every chain input is a construction, so this is a LOWER bound on
# ERS.  Below it a claim is definitely not an improvement; above it a claim is not thereby safe.
# n_0 -> largest k with d >= ceil(n_0/4) among best known binary linear codes,
# codetables.de, read 2026-08-23.  Kept only as a CROSS-CHECK on the chain: level 0 is
# one term of the chain, so best_chain must never come out below it.
ERS0 = {48: 24, 52: 21, 56: 23, 60: 25, 61: 25, 62: 26, 63: 27, 64: 28, 68: 25,
        72: 28, 76: 26, 79: 29, 80: 30, 84: 27, 88: 28, 92: 30, 93: 31, 94: 31, 95: 32}
_ers = os.path.join(HERE, 'verifications', 'closed', 'dim96-ers-takeover')
sys.path.insert(0, _ers)
try:
    from ers_exposure import best_chain          # noqa: E402
except Exception as _e:                          # pragma: no cover - reported, never silent
    best_chain = None
    bad("cannot import ers_exposure.best_chain (%s); the ERS floor is NOT being checked" % _e)


def _d95_threshold():
    """The dimension-95 exposure, DERIVED, for the '96 and above' row of RESULTS.md.

    This used to be three transcribed numbers in the row text, and they went stale twice in
    one night: once when the class family moved the claim, and once more because ers96.py was
    deriving its own copy from a hardcoded CLAIMED dict.  The threshold depends on the claim,
    so anything that moves a claim moves it too.  Note it is a fact about ONE CHAIN through
    the dimension, not about the dimension -- a bare "A(95,23,23) < N" reads as the latter --
    so the chain is named in the text.
    """
    from ers_exposure import ranked               # noqa: E402
    row = next((r for r in ranked() if r[0] == 95), None)
    assert row is not None, "dimension 95 is not in the ERS ranking"
    _n, _c, _v, ch, w1, have, need, _ratio, _why, _isers = row

    def sp(v):
        return "{:,}".format(int(v)).replace(",", " ")

    return ("A(%d,%d,%d) < %s" % (_n, w1, w1, sp(need)), sp(have),
            "(%s)" % ", ".join(str(x) for x in ch))
if best_chain is not None:
    beaten, at_floor = [], []
    for d, v, pkg, how in claims:
        floor = int(best_chain(d)[0])
        lvl0 = max([1 << k for n0, k in ERS0.items() if n0 <= d] or [0])
        if floor < lvl0:                         # the chain contains level 0, so this cannot be
            bad("dim %d: the chain floor %d is below the level-0 floor %d -- the two "
                "computations disagree" % (d, floor, lvl0))
        if floor > v:
            beaten.append((d, v, floor))
        elif floor == v:
            at_floor.append(d)
    for d, v, floor in beaten:
        bad("dim %d: claimed %d but the ERS chain already reaches %d from published codes "
            "-- NOT an improvement" % (d, v, floor))
    if not beaten:
        print("   all %d claims exceed it, or equal it by BEING it" % len(claims))
        if at_floor:
            print("   %s equal it: %s"
                  % (runs(at_floor), ', '.join(sorted({p for d, _v, p, _h in claims
                                                       if d in at_floor}))))
        print("   this is a LOWER bound on ERS, so clearing it is necessary, not sufficient;")
        print("   ers_exposure.py reports how much room each claim has (least: dimension 68)")

print()
print("5. every package exists and is documented")
for d, v, pkg, how in claims:
    # Usually improved/, but dim96-ers-takeover sits under closed/: it was written to settle
    # dimension 96 AGAINST the cap construction, and it still does that -- the claim there is
    # Edel-Rains-Sloane's, not this paper's.  Look in both rather than move it.
    p = os.path.join(IMP, pkg)
    if not os.path.isdir(p):
        p = os.path.join(HERE, 'verifications', 'closed', pkg)
    if not os.path.isdir(p):
        bad("dim %d: package %s does not exist" % (d, pkg))
    elif not os.path.exists(os.path.join(p, 'README.md')):
        bad("dim %d: package %s has no README.md" % (d, pkg))

# The tiers are whatever verifications/ holds -- closed, improved, recovered, superseded
# today.  Naming them here instead cost check 5a and 5b the superseded tier entirely.
_VDIR = os.path.join(HERE, 'verifications')
_TIERS = sorted(_t for _t in os.listdir(_VDIR)
                if os.path.isdir(os.path.join(_VDIR, _t))
                and not _t.startswith(('.', '__')))

def _named(fname, prose):
    """Is fname mentioned in prose in its own right?

    A plain `fname in prose` counts capdirs.py as documented when only verify_capdirs.py is
    mentioned, k4.py when only scan_k4.py is, triples.py when only lamtriples.py is -- and
    the second of those was hiding a real gap.  Require the character before the name not to
    be one that could continue a filename.
    """
    return re.search(r'(?<![A-Za-z0-9_])' + re.escape(fname), prose) is not None

print()
print("5a. every script in every package is named in that package's README")
# A reviewer reads one README and then lists the directory.  Anything in the listing the
# README never mentions is code with no stated role -- and the two jobs run_all.py reports
# as skipped were exactly that, so the run said "skipped: ... (needs scripts/classes)" and
# the README named neither the script nor what would let it run.
_undoc = []
for _tier in _TIERS:
    _tdir = os.path.join(HERE, 'verifications', _tier)
    for _pkg in sorted(os.listdir(_tdir)):
        _pdir = os.path.join(_tdir, _pkg)
        _rd = os.path.join(_pdir, 'README.md')
        if not os.path.isdir(_pdir) or not os.path.exists(_rd):
            continue
        _prose = io.open(_rd, encoding='utf-8', errors='replace').read()
        for _dp, _dn, _fn in os.walk(_pdir):
            _dn[:] = [_x for _x in _dn if _x != '__pycache__']
            for _f in sorted(_fn):
                if _f.endswith('.py') and not _named(_f, _prose):
                    _undoc.append('%s/%s: %s'
                                  % (_tier, _pkg,
                                     os.path.relpath(os.path.join(_dp, _f), _pdir)
                                     .replace(os.sep, '/')))
# The top-level scripts too.  update.py -- the maintainer's entry point, which runs every
# driver, the audit, the table generators, the paper's checkers and the whole suite -- was
# named in no markdown file at all, because this check looked only inside packages.
_rtext = io.open(os.path.join(HERE, 'README.md'), encoding='utf-8', errors='replace').read()
for _f in sorted(_f for _f in os.listdir(HERE) if _f.endswith('.py')):
    if not _named(_f, _rtext):
        _undoc.append('(top level): %s' % _f)

if _undoc:
    for _u in _undoc[:12]:
        bad("no README names this script: %s" % _u)
    if len(_undoc) > 12:
        bad("... and %d more" % (len(_undoc) - 12))
else:
    print("   every .py in every package, and beside them, is named in a README")

print()
print("5b. every package that ships a SHA-256 manifest still matches it")
# A documentation edit is enough to break one: adding verify_manifest.py to dim27's README
# did it, and audit.py passed, because only run_all.py runs the package's own checker.
_manifests = _manbad = 0
for _tier in _TIERS:
    _tdir = os.path.join(HERE, 'verifications', _tier)
    for _pkg in sorted(os.listdir(_tdir)):
        _man = os.path.join(_tdir, _pkg, 'data', 'SHA256SUMS.txt')
        if not os.path.exists(_man):
            continue
        _manifests += 1
        _root = os.path.join(_tdir, _pkg)
        for _line in io.open(_man, encoding='utf-8').read().split('\n'):
            if not _line.strip():
                continue
            _want, _name = _line.split('  ', 1)
            _f = os.path.join(_root, _name)
            if not os.path.exists(_f):
                _manbad += 1
                bad("%s: manifest names a file that is gone: %s" % (_pkg, _name))
            elif hashlib.sha256(io.open(_f, 'rb').read()).hexdigest() != _want:
                _manbad += 1
                bad("%s: %s no longer matches the manifest" % (_pkg, _name))
if not _manbad:
    print("   %d manifest%s, every entry matching"
          % (_manifests, '' if _manifests == 1 else 's'))

print()
print("5c. the working directory's README accounts for every path beside it")
# It named five and there were seven: colab/ and gpu/ had appeared since it was written,
# and were missing from .gitignore too.  Skipped when this package travels on its own.
_root = os.path.normpath(os.path.join(HERE, os.pardir))
_rootmd = os.path.join(_root, 'README.md')
if not os.path.exists(_rootmd):
    print("   (no working directory beside this package; nothing to check)")
else:
    _rtext = io.open(_rootmd, encoding='utf-8', errors='replace').read()
    _unnamed = [_e for _e in sorted(os.listdir(_root))
                if not _e.startswith('.') and _e != 'README.md' and _e not in _rtext]
    if _unnamed:
        for _e in _unnamed:
            bad("%s is beside this package and its README never names it" % _e)
    else:
        print("   every path beside this package is named in its README")

print()
print("5d. verifications/README.md counts the tiers and packages it actually has")
# It said five tiers, seven packages in improved/, eight in closed/, and 75 sections of
# KNOWLEDGE.md, against four, nine, ten and 96 on disk -- and its improved/ table had no row
# for dim62-63-ers-chain or dim68-69-gamma72-cross-sections at all.  Each of those numbers
# counts something, so each can be read back and compared.
_WORDS = {'one': 1, 'two': 2, 'three': 3, 'four': 4, 'five': 5, 'six': 6, 'seven': 7,
          'eight': 8, 'nine': 9, 'ten': 10, 'eleven': 11, 'twelve': 12}
_vr = os.path.join(_VDIR, 'README.md')
if not os.path.exists(_vr):
    bad("verifications/README.md is missing")
else:
    _vt = io.open(_vr, encoding='utf-8', errors='replace').read()
    _vrbad = 0
    _flat = ' '.join(_vt.split())

    def _stated(pat):
        """The number the README states in a phrase like 'Nine packages', word or digit."""
        _m = re.search(pat, _flat, re.I)
        if not _m:
            return None
        _tok = _m.group(1)
        return int(_tok) if _tok.isdigit() else _WORDS.get(_tok.lower())

    for _pat, _have, _what in (
            (r'([A-Za-z]+|\d+) tiers', len(_TIERS), 'tiers'),
            (r'\(improved/\)[^|]*?\n?.{0,80}?([A-Za-z]+|\d+) packages',
             len(os.listdir(os.path.join(_VDIR, 'improved'))), 'packages in improved/'),
            (r'\(closed/\)[^|]*?\n?.{0,120}?([A-Za-z]+|\d+) packages',
             len(os.listdir(os.path.join(_VDIR, 'closed'))), 'packages in closed/')):
        _said = _stated(_pat)
        if _said is None:
            _vrbad += 1
            bad("verifications/README.md no longer states a count of %s" % _what)
        elif _said != _have:
            _vrbad += 1
            bad("verifications/README.md says %d %s; there are %d" % (_said, _what, _have))

    # every package directory must have a row in that file
    _norow = []
    for _tier in _TIERS:
        for _pkg in sorted(os.listdir(os.path.join(_VDIR, _tier))):
            if os.path.isdir(os.path.join(_VDIR, _tier, _pkg)) and _pkg not in _vt:
                _norow.append('%s/%s' % (_tier, _pkg))
    if _norow:
        _vrbad += len(_norow)
        for _n in _norow:
            bad("verifications/README.md has no row for %s" % _n)

    # and the KNOWLEDGE.md section count it quotes
    _kp = os.path.join(HERE, 'KNOWLEDGE.md')
    _ksaid = re.search(r'KNOWLEDGE\.md\) . (\d+) sections', _flat)
    if os.path.exists(_kp) and _ksaid:
        _khave = len(re.findall(r'^## ', io.open(_kp, encoding='utf-8',
                                                 errors='replace').read(), re.M))
        if int(_ksaid.group(1)) != _khave:
            _vrbad += 1
            bad("verifications/README.md says KNOWLEDGE.md has %s sections; it has %d"
                % (_ksaid.group(1), _khave))
    # KNOWLEDGE.md's OWN first paragraph counts itself, and nothing was checking that one:
    # it said 75 sections and about 4200 lines while the file stood at 96 and 5799.  A number
    # in a file's own opening sentence is the kind no other file has a reason to look at.
    if os.path.exists(_kp):
        _ktext = io.open(_kp, encoding='utf-8', errors='replace').read()
        _khave = len(re.findall(r'^## ', _ktext, re.M))
        _klines = len(_ktext.rstrip(chr(10)).split(chr(10)))
        _kself = re.search(r'\*\*(\d+) sections, about (\d+) lines\*\*', _ktext)
        if not _kself:
            _vrbad += 1
            bad("KNOWLEDGE.md no longer states its own size in the form the audit checks")
        else:
            if int(_kself.group(1)) != _khave:
                _vrbad += 1
                bad("KNOWLEDGE.md says it has %s sections; it has %d"
                    % (_kself.group(1), _khave))
            if abs(int(_kself.group(2)) - _klines) > 100:
                _vrbad += 1
                bad("KNOWLEDGE.md says about %s lines; it has %d"
                    % (_kself.group(2), _klines))
    if not _vrbad:
        print("   tiers, package counts, every package row, and KNOWLEDGE.md all agree")
        print("   (KNOWLEDGE.md's own header counts itself correctly too)")

print()
print("5e. every path a .gitignore rule names still exists")
# verifications/open/dim17-layered-family/... survived the package's move to closed/, so the
# cache it excluded was shipping again.  A rule that matches nothing is indistinguishable
# from a rule with nothing to match, which is why this is checked rather than noticed.
_gi = os.path.join(HERE, '.gitignore')
_gone = []
if os.path.exists(_gi):
    for _rule in io.open(_gi, encoding='utf-8'):
        _rule = _rule.strip()
        if not _rule or _rule.startswith(('#', '!')) or '/' not in _rule.rstrip('/'):
            continue                      # a bare pattern names no directory
        _dirpart = _rule.rstrip('/').rsplit('/', 1)[0]
        if any(_c in _dirpart for _c in '*?['):
            continue                      # wildcards in the directory part: not a fixed path
        if not os.path.exists(os.path.join(HERE, _dirpart)):
            _gone.append((_rule, _dirpart))
if _gone:
    for _r, _d in _gone:
        bad(".gitignore rule %r names %s, which does not exist" % (_r, _d))
else:
    print("   every rooted rule points at a directory that is there")

# Files outside this directory that state the claim count, named one at a time.  Not the
# whole repository root: research/ and notes/ are deliberately historical and the root
# README says as much, so sweeping them would fail on text that is correct as history.
def _paper_dir():
    """the paper's directory, inside this package or beside it -- whichever holds it

    It shipped INSIDE on 2026-08-30 so that its checkers run in a clone; before that it was
    a sibling.  Assuming one layout is how this file broke a clone twice, so both are tried
    and neither is required."""
    for _c in (os.path.join(HERE, 'paper'), os.path.join(HERE, os.pardir, 'paper')):
        if os.path.exists(os.path.join(_c, 'kissing46.tex')):
            return _c
    return None


_PAPERDIR = _paper_dir()
_ALSO = [os.path.join(HERE, os.pardir, 'README.md')] + (
    [os.path.join(_PAPERDIR, 'kissing46.tex')] if _PAPERDIR else [])
_walked = [(os.path.dirname(_q), [], [os.path.basename(_q)])
           for _q in _ALSO if os.path.exists(_q)]

# The text sweeps of 5l, 5m and 5n used to walk os.path.dirname(HERE) -- the PARENT of this
# package.  In the author's working folder that is right, because paper/ is there.  In a
# CLONE it is whatever directory the clone was made in, so the audit swept the reviewer's
# unrelated files: a scratch script holding a set literal named TAU crashed 5n outright, on
# the first command the README gives.  Scan this package, and the paper beside it only when
# that paper is THIS one.  Checks 5f and 6 above already had it right.
_ROOTS = [HERE] if _PAPERDIR is None or _PAPERDIR.startswith(HERE) else [HERE, _PAPERDIR]


def _sweep():
    """os.walk over this package and, when it is beside it, the paper -- and nothing else"""
    for _r in _ROOTS:
        for _tup in os.walk(_r):
            yield _tup

print()
print("5f. every pointer into an outside document says which document")
# "Equation (3)" of tau(47) = 23766960 was the preprint's in two files and the journal
# version's in four.  Both were plausible, neither was checked, and two versions need not
# number alike -- so a referee following the manuscript could land somewhere else.  Nothing
# numeric can catch it: the value is recomputed here and agrees either way.  What is
# checkable is that a numbered pointer names its document in the same passage.
_POINTERS = (
    # (phrase, the identifiers that fix which document its numbering belongs to)
    ("equation (3)", ("2312.05121",)),
    ("Table 3 of Sun", ("2607.20359",)),
)
_loose = []
for _dp, _dn, _fn in list(os.walk(HERE)) + _walked:
    _dn[:] = [d for d in _dn if d not in ('__pycache__', '.git', 'classes', 'state')]
    for _f in _fn:
        if not _f.endswith(('.md', '.py', '.tex')) or _f == 'audit.py':
            continue                       # audit.py quotes the phrases it looks for
        _q = os.path.join(_dp, _f)
        try:
            _t = io.open(_q, encoding='utf-8').read()
        except (IOError, UnicodeDecodeError):
            continue
        # LaTeX writes equation~(3) with a tie, and prose wraps mid-phrase, so the
        # tie and every run of whitespace collapse to one space before matching.
        # Without this the check could not fire on the manuscript at all.
        _flat = re.sub(r'[~\s]+', ' ', _t)
        _low = _flat.lower()
        for _phrase, _ids in _POINTERS:
            _i = _low.find(_phrase.lower())
            while _i >= 0:
                _win = _flat[max(0, _i - 300):_i + 300]
                if not any(_id in _win for _id in _ids):
                    _loose.append((os.path.relpath(_q, HERE), _phrase))
                _i = _low.find(_phrase.lower(), _i + 1)
if _loose:
    for _q, _ph in sorted(set(_loose)):
        bad("%s says %r without naming the document it numbers" % (_q, _ph))
else:
    print("   every numbered pointer into an outside work names its document")

print()
print("5g. the README's description of the paper still partitions the claims")
# Three stale facts in one sentence -- a page count, and two method-to-dimension lists that
# had not followed 62, 63, 68 and 69 when they moved.  Each was true when written.
# The page count in the same sentence is checked by paper/factcheck.py, over both READMEs
# that state it, so the PDF is opened in one place only.  This check needs no PDF -- it used
# to branch on one, and so skipped itself in exactly the clean checkout where it matters.
_rm = os.path.join(HERE, 'README.md')
_wr = io.open(_rm, encoding='utf-8').read()
_m = re.search(r'is a (\d+)-page account of(.{0,600}?)\n\n', _wr, re.S)
if not _m:
    bad("README.md has no 'is a N-page account of' sentence to check")
else:
    # every dimension the sentence names, once each, minus the two it calls recovered
    _body = re.sub(u'\u00a7\\d+(\\.\\d+)?', ' ', _m.group(2))     # drop section numbers
    _seen, _dupe = [], []
    for _a, _b, _c in re.findall(u'(\\d+)\\s*[-\u2013]\\s*(\\d+)|(\\d+)', _body):
        _seen.extend([int(_c)] if _c else list(range(int(_a), int(_b) + 1)))
    _dupe = sorted(set(d for d in _seen if _seen.count(d) > 1))
    _named = set(_seen) - {46, 47}
    _claimed = set(d for d, v, pkg, how in claims)
    if _dupe:
        bad("README.md's description lists %s under more than one method" % _dupe)
    elif _named != _claimed:
        bad("README.md's description names %s and omits %s"
            % (sorted(_named - _claimed) or 'nothing extra',
               sorted(_claimed - _named) or 'nothing'))
    else:
        print("   its dimension lists partition the %d claims" % len(_claimed))

print()
print("5h. the root README names every script that reads Cohn's coordinate files")
# It named two.  Four read them, and one of the omitted two is frame17.py, which run_all.py
# itself labels "needs Cohn's data set".  A mention in a comment is not a read, so the scan
# requires the name AND a path resolution beside it.
_ROOT = os.path.join(HERE, os.pardir)
_SETS = ('dimensions1-24.txt', 'dimensions25-31.txt', 'dimensions32-47.txt')
_readers = set()
for _dp, _dn, _fn in os.walk(HERE):
    _dn[:] = [d for d in _dn if d not in ('__pycache__', '.git', 'classes', 'state')]
    for _f in _fn:
        if not _f.endswith('.py') or _f == 'audit.py':
            continue                       # audit.py quotes the names this check looks for
        _q = os.path.join(_dp, _f)
        try:
            _t = io.open(_q, encoding='utf-8').read()
        except (IOError, UnicodeDecodeError):
            continue
        if not any(_s in _t for _s in _SETS):
            continue
        # a reader resolves a path: argv, an environment variable, or os.path.join to one
        if re.search(r"(sys\.)?argv|environ|os\.path\.join\([^)]*dimensions", _t):
            # the PACKAGE, not the file: two of them are read through a shared extract.py
            # behind a registered entry point, and it is packages the README names.
            _rel = os.path.relpath(_q, os.path.join(HERE, 'verifications'))
            _readers.add(_rel.replace(chr(92), '/').split('/')[1])
_rrp = os.path.join(_ROOT, 'README.md')
if not os.path.exists(_rrp):
    # A clone of this directory alone: the outer working folder is not part of the
    # deliverable, so there is no root README to check against.  Say so rather than
    # crash, which is what this did -- and audit.py is the FIRST job of run_all.py, so
    # the crash was the first thing a reviewer of a clean checkout would have seen.
    print("   the outer working folder is not here, so there is no root README to check; "
          "%d reader(s) found" % len(_readers))
    _rr, _readers = None, set()
else:
    _rr = io.open(_rrp, encoding='utf-8').read()
_missing = sorted(f for f in _readers if _rr is not None and f not in _rr)
if _rr is None:
    pass                               # no root README here; already reported above
elif _missing:
    bad("the root README does not name the package(s) %s, which read Cohn's "
        "data sets" % ', '.join(_missing))
else:
    print("   all %d of them are named there" % len(_readers))
_said = re.search(r'(\w+) scripts in the deliverable read them', _rr) if _rr else None
if _rr is None:
    pass                               # no root README here; already reported above
elif not _said:
    bad("the root README no longer states how many scripts read the data sets")
else:
    _WORDS = {'one': 1, 'two': 2, 'three': 3, 'four': 4, 'five': 5, 'six': 6, 'seven': 7}
    if _WORDS.get(_said.group(1).lower()) != len(_readers):
        bad("the root README says %s scripts read Cohn's data sets; %d packages do"
            % (_said.group(1), len(_readers)))
    else:
        print("   and the count it states, %s, is the number found" % _said.group(1))

print()
print("5i. update.py's improvement guide lists the right dimensions per status")
# It named 39, 62 and 63 as `external` and omitted 96, which has been external since
# 2026-08-26 and is the largest claim here.  status.py decides the groups; read the guide
# against it.  Only the clause after the word "dimensions" is parsed -- the `classes` entry
# runs on for eight more lines of unrelated numbers.
_up = os.path.join(HERE, 'update.py')
if not os.path.exists(_up):
    bad("update.py is missing")
else:
    _ut = io.open(_up, encoding='utf-8').read()
    _groups = {}
    for _d, _v, _pkg, _how in claims:
        _groups.setdefault(status_for(_d)[0], set()).add(_d)
    _gbad = 0
    for _st in sorted(_groups):
        _m = re.search(re.escape(_st) + r'\s+[^\n]*?dimensions? ([-0-9,\s]+(?:and [-0-9]+)?)\.',
                       _ut)
        if not _m:
            bad("update.py's guide does not list the dimensions of `%s`" % _st)
            _gbad += 1
            continue
        _got = set()
        for _a, _b, _c in re.findall(r'(\d+)\s*-\s*(\d+)|(\d+)', _m.group(1)):
            _got |= set(range(int(_a), int(_b) + 1)) if _a else {int(_c)}
        if _got != _groups[_st]:
            bad("update.py's guide gives `%s` as %s; the claims give %s"
                % (_st, runs(sorted(_got)), runs(sorted(_groups[_st]))))
            _gbad += 1
    if not _gbad:
        print("   all %d status groups match the claim table" % len(_groups))

print()
print("5j. this script's docstring lists every check this script runs")
# The docstring stopped at check 5 while 5a-5i and 6 also ran, so a reader starting at the
# top was told about less than half of it.  Every heading printed below is now required to
# appear in __doc__ -- headings are what a reader navigates by.
_printed = set(re.findall(r'^print\("(\d+[a-z]?)\.', io.open(__file__, encoding='utf-8').read(),
                          re.M))
_undoc = sorted(_h for _h in _printed
                if not re.search(r'^\s*%s[.-]' % re.escape(_h), __doc__ or '', re.M)
                and not re.search(r'%s' % re.escape(_h), __doc__ or ''))
if _undoc:
    bad("the docstring does not mention check(s) %s" % ', '.join(_undoc))
else:
    print("   all %d printed headings are in the docstring" % len(_printed))

print()
print("5k. .gitattributes states the line-ending mix it actually has")
# It claimed the tree was uniformly LF and that an earlier mix had been normalised.  It had
# not been.  Nothing rests on uniformity -- `* -text` preserves bytes, which is what the
# dimension-27 manifest needs -- but a maintainer who believed it might drop `* -text`.
_ga = os.path.join(HERE, '.gitattributes')
_CR, _LF = chr(13) + chr(10), chr(10)
_TEXT = ('.py', '.md', '.txt', '.cff', '.json', '.tex', '.bib', '.yml', '.cfg', '.out')
_nlf = _ncrlf = 0
_mix = []
for _dp, _dn, _fn in os.walk(HERE):
    _dn[:] = [_d for _d in _dn if _d not in ('__pycache__', '.git')]
    for _f in _fn:
        if not _f.endswith(_TEXT):
            continue
        try:
            _t = io.open(os.path.join(_dp, _f), encoding='utf-8', newline='').read()
        except (IOError, UnicodeDecodeError):
            continue
        _c = _t.count(_CR)
        if _c and _t.count(_LF) - _c:
            _mix.append(os.path.relpath(os.path.join(_dp, _f), HERE))
        elif _c:
            _ncrlf += 1
        else:
            _nlf += 1
if _mix:
    bad("these files mix line endings within themselves: %s" % ', '.join(sorted(_mix)))
if not os.path.exists(_ga):
    bad(".gitattributes is missing; checkout would rewrite line endings")
else:
    _gt = io.open(_ga, encoding='utf-8').read()
    if not re.search(r'^\*\s+-text\s*$', _gt, re.M):
        bad(".gitattributes no longer disables end-of-line conversion (`* -text`)")
    _said = re.search(r'(\d+) text files use LF and (\d+) use CRLF', _gt)
    if not _said:
        bad(".gitattributes no longer states the line-ending mix")
    elif (int(_said.group(1)), int(_said.group(2))) != (_nlf, _ncrlf):
        bad(".gitattributes says %s LF and %s CRLF; the tree has %d and %d"
            % (_said.group(1), _said.group(2), _nlf, _ncrlf))
    else:
        print("   %d LF, %d CRLF, none of those mixed, over the %d extensions this "
              "scans; `* -text` preserves every byte" % (_nlf, _ncrlf, len(_TEXT)))

# ---- 5l. LaTeX escapes that were halved into control characters ------------------------
# "\tau" written through a non-raw Python string is TAB + "au", "\frac" is FF + "rac".
# The tex is scanned for this by gpu/propagate.py; the markdown was not, and three files
# carried it -- both "\tau"s in the dimension-38 README, one in the dimension-39 README,
# and the sentence in paper/README.md that describes the fault.  The signature is narrow
# on purpose: a control character immediately followed by two or more letters.
print()
print("5l. no LaTeX escape has been halved, and no Cyrillic look-alike letter")
# The second half was added on 2026-08-30, after `\u0413` (Cyrillic Ghe) was found
# standing for Gamma in the submission notes.  It renders exactly like a Greek Gamma,
# so a proofread cannot see it and a grep for Gamma does not find the sentence.  Runs
# of Cyrillic are legitimate -- KNOWLEDGE.md quotes a Russian journal title -- so only
# an ISOLATED Cyrillic letter counts.
_CYR = re.compile('[' + chr(0x0400) + '-' + chr(0x04FF) + ']+')
_CTRL = {chr(7): 'a', chr(8): 'b', chr(9): 't', chr(11): 'v', chr(12): 'f'}
_halved = []
for _dp, _dn, _fn in _sweep():
    _dn[:] = [_d for _d in _dn if _d not in ('__pycache__', '.git', 'research', 'colab')]
    for _f in _fn:
        if not _f.endswith(('.md', '.tex')):
            continue
        _p = os.path.join(_dp, _f)
        try:
            _t = io.open(_p, encoding='utf-8').read()
        except (IOError, UnicodeDecodeError):
            continue
        for _i, _L in enumerate(_t.split(chr(10)), 1):
            for _c, _n in _CTRL.items():
                for _m in re.finditer(re.escape(_c) + '[A-Za-z]{2,}', _L):
                    _halved.append('%s:%d looks like a lost backslash before %s%s'
                                   % (os.path.relpath(_p, os.path.dirname(HERE)), _i,
                                      _n, _m.group(0)[1:]))
        for _m in _CYR.finditer(_t):
            if len(_m.group()) > 1:
                continue                       # a Russian word quoted on purpose
            _halved.append('%s:%d has a lone Cyrillic U+%04X where a Greek or Latin letter '
                           'is meant' % (os.path.relpath(_p, os.path.dirname(HERE)),
                                         _t.count(chr(10), 0, _m.start()) + 1,
                                         ord(_m.group())))
if _halved:
    for _h in _halved[:8]:
        bad(_h)
else:
    print("   every .md and .tex outside research/ is free of both")

# ---- 5m. no table anywhere quotes a superseded value of a claim ------------------------
# A textual search for the old number is not enough: the package README separates thousands
# with U+202F, so grepping "361 275 383" missed it while it was on the page.  This is
# relational instead -- any row naming a claimed dimension and quoting a number within 3%
# of that claim must quote the claim exactly.  A floor or an ERS total is never that close.
print()
print("5m. no markdown table quotes a superseded value of a claim")
# In --write-results mode the claims are about to change and the documents have not been
# brought along yet, so this is a regeneration step, not a defect: run it only in the
# checking pass.  Otherwise a claim could never be raised -- the write would be blocked by
# documents that cannot be fixed until the write has happened.
if '--write-results' in sys.argv:
    print("   skipped: this run REGENERATES the claims the documents are compared against")
    _stale = []
else:
    _SEPS = (chr(0x202f), chr(0xa0), chr(0x2009), ',', chr(92) + ',')


    def _cell(_x):
        _x = _x.strip().replace('**', '').replace('`', '')
        for _sep in _SEPS:
            _x = _x.replace(_sep, ' ')
        _x = _x.replace(' ', '')
        # not .isdigit(): that is True for superscripts like 2³⁹, which int() then rejects
        return int(_x) if re.match(r'^[0-9]+$', _x) else None


    _canon = dict((_d, _v) for _d, _v, _p, _h in claims)
    _stale = []
    for _dp, _dn, _fn in _sweep():
        _dn[:] = [_x for _x in _dn if _x not in ('__pycache__', '.git', 'research', 'colab')]
        for _f in _fn:
            # RESULTS.md is the generated authority these claims are read FROM, and in
            # --write-results mode it is the output; comparing it with itself would deadlock a
            # regeneration on its own staleness.  Step 5 already compares it, claim by claim.
            if not _f.endswith('.md') or _f in ('KNOWLEDGE.md', 'RESULTS.md'):
                continue
            _p = os.path.join(_dp, _f)
            try:
                _lines = io.open(_p, encoding='utf-8').read().split(chr(10))
            except (IOError, UnicodeDecodeError):
                continue
            for _i, _L in enumerate(_lines, 1):
                if _L.count('|') < 3:
                    continue
                _parts = [_x.strip() for _x in _L.split('|')]
                _d = _cell(_parts[1]) if len(_parts) > 1 else None
                if _d not in _canon:
                    continue
                _vals = [_cell(_x) for _x in _parts[2:]]
                _vals = [_v for _v in _vals if _v is not None and _v >= 1000]
                if _canon[_d] in _vals:
                    continue
                for _v in _vals:
                    if _canon[_d] * 0.97 <= _v <= _canon[_d] * 1.03:
                        _stale.append('%s:%d quotes %d for dimension %d; the claim is %d'
                                      % (os.path.relpath(_p, os.path.dirname(HERE)), _i,
                                         _v, _d, _canon[_d]))
                        break
if _stale:
    for _t in _stale[:8]:
        bad(_t)
else:
    print("   every dimension row outside KNOWLEDGE.md and research/ is current")

# ---- 5n. no private tau table disagrees with common/published.py -----------------------
# axis_rotate.py and verify_poles.py both carried tau(20) = 17400 and tau(21) = 27720 until
# 2026-08-27.  Those are Lambda_20 and Lambda_21; the best known configurations in 20 and 21
# are not lattices, and both entries were 2048 low.  Nothing computed used them -- they set a
# printed column and an "AT THE CEILING" note -- so nothing failed, and the repository simply
# stated tau(20) two ways.  Retyping the two numbers would fix this instance; reading
# published.py fixes the class, and this check keeps the next copy from being retyped.
print()
print("5n. no script carries a tau table that disagrees with common/published.py")
_drift = []
for _dp, _dn, _fn in _sweep():
    _dn[:] = [_d for _d in _dn if _d not in ('__pycache__', '.git', 'research', 'colab')]
    for _f in _fn:
        if not _f.endswith('.py'):
            continue
        _p = os.path.join(_dp, _f)
        try:
            _t = io.open(_p, encoding='utf-8').read()
        except (IOError, UnicodeDecodeError):
            continue
        for _m in re.finditer(r'TAU\s*=\s*\{([^{}]*)\}', _t):
            _body = _m.group(1)
            if '"' in _body or chr(39) in _body:      # keyed by root-system name, not by k
                continue
            try:
                _d = eval('{' + _body + '}', {'__builtins__': {}}, {})
            except Exception:
                continue
            # eval of a brace body gives a SET for {1, 2, 3} and a dict for {1: 2}, and
            # .items() was called before anything checked which.  A set literal named TAU
            # anywhere in the swept tree crashed this, clone or not.
            if not isinstance(_d, dict) or not _d or not all(
                    isinstance(_k, int) and isinstance(_v, int) for _k, _v in _d.items()):
                continue
            for _k in sorted(_d):
                if _k in PUBTAU and _d[_k] != PUBTAU[_k]:
                    _drift.append('%s says tau(%d) = %d; common/published.py says %d'
                                  % (os.path.relpath(_p, os.path.dirname(HERE)),
                                     _k, _d[_k], PUBTAU[_k]))
if _drift:
    for _x in _drift[:8]:
        bad(_x)
else:
    print("   every tau table in the tree agrees with the published one, entry for entry")

# ---- 5o. the rotated axis layers lie in the axis space --------------------------------
# Dimension 72+k puts the cap directions AND the poles in one R^k.  Several direction sets
# are stored in more coordinates than they span -- Lambda_21's cross-section is 27 720
# vectors in Lambda_24's 24 coordinates spanning 21 -- and the first rotations built for
# those layers were isometries of the STORED coordinates, which carried the poles out of
# R^k.  Every inner product was right; the configuration needed R^96 instead of R^93.  The
# rotation is linear, so checking a basis of the span settles it, which is why this costs
# nothing and can sit in the fast audit.
print()
print("5o. every rotated axis layer lies in the axis space")
_pkg = os.path.join(HERE, 'verifications', 'improved', 'dim73-95-gamma72-caps')
if os.path.exists(os.path.join(_pkg, 'scripts', 'axis_block.py')):
    sys.path.insert(0, os.path.join(_pkg, 'scripts'))
    from axis_block import directions as _dirs, preserves_span as _pspan, rows as _rows
    _R = _rows()
    _leaks = []
    for _k in [_x for _x in range(9, 24) if _x != 11]:
        _f = os.path.join(_pkg, 'data', 'poles_rot_%d.npz' % _k)
        if not os.path.exists(_f):
            continue
        _W, _c, _m = _dirs(_k, _R)
        if not _pspan(_W, np.load(_f, allow_pickle=True)['N']):
            _leaks.append('k=%d: the rotation carries the poles out of R^%d, so dimension '
                          '%d does not hold' % (_k, _k, 72 + _k))
    if _leaks:
        for _x in _leaks:
            bad(_x)
    else:
        print("   all %d rotated layers are inside R^k, checked on a basis of the span"
              % len([1 for _k in range(9, 24)
                     if os.path.exists(os.path.join(_pkg, 'data', 'poles_rot_%d.npz' % _k))]))
else:
    bad("scripts/axis_block.py is missing: the axis-space check cannot run")

# ---- 5p. verify_poles.py's docstring states the size of the layers; is it still true? --
# The same claim in the PAPER is guarded by paper/factcheck.py, which reads the percentage
# out of the text and brackets it with the shipped layers.  This one sat in a script header
# with nothing reading it, and script headers are where a measured number goes stale in
# silence: it claimed "between 86% and 100% ... above 98% once k >= 12" while the layers
# moved twice underneath it.  Same treatment: derive, then parse.
print()
print("5p. the layer percentages in verify_poles.py's header are still the measured ones")
_vp = os.path.join(_pkg, 'scripts', 'verify_poles.py')
if os.path.exists(_vp):
    # Self-contained: 5o's imports live inside its own branch, and a check that runs only
    # when the check above it did is a check that can disappear without saying so.
    sys.path.insert(0, os.path.join(_pkg, 'scripts'))
    from axis_block import directions as _dirs5p, rows as _rows5p
    _R5p = _rows5p()
    _pcs = {}
    for _k in range(3, 24):
        # The layer file follows the CAP SET, not the dimension.  At k = 12 and 13 the
        # directions are the Kappa sections K_12 and K_13 and the layer is
        # poles_kap_<k>.npz; poles_rot_<k>.npz there would be a layer of Lambda_k, whose
        # indices point at vectors this W does not contain -- a percentage computed from
        # the pair would be arithmetic performed on two unrelated objects.
        _names = (('poles_kap_%d.npz' % _k,) if _R5p.get(_k, {}).get('W') == 'kappa'
                  else ('poles_rot_%d.npz' % _k, 'poles_rot_k%d.npz' % _k))
        _f = [_x for _x in _names if os.path.exists(os.path.join(_pkg, 'data', _x))]
        if not _f:
            continue
        _lay = len(np.load(os.path.join(_pkg, 'data', _f[0]))['keep'])
        if _k <= 8:
            # the direction set IS the root system, and |roots| = tau(k) for A3 .. E8
            _nz = PUBTAU[_k]
        else:
            _W5p, _c5p, _m5p = _dirs5p(_k, _R5p)
            _nz = len(_W5p)
        _pcs[_k] = 100.0 * _lay / _nz
    _flat5p = ' '.join(io.open(_vp, encoding='utf-8').read().split())
    _mm = re.search(r'keep between (\d+)% and 100% of the direction set; the caps are '
                    r'small and the smallest layer, at k = (\d+), is ([0-9.]+)% of it',
                    _flat5p)
    if not _mm:
        bad("verify_poles.py no longer states the layer floor in the shape 5p reads")
    else:
        _said, _saidk, _saidpc = (int(_mm.group(1)), int(_mm.group(2)),
                                  float(_mm.group(3)))
        _worst = min(_pcs, key=lambda _x: _pcs[_x])
        if not (_said <= _pcs[_worst] < _said + 1):
            bad("verify_poles.py says the layers keep at least %d%%; the smallest is %.2f%%"
                % (_said, _pcs[_worst]))
        elif _saidk != _worst or abs(_saidpc - _pcs[_worst]) > 0.005:
            bad("verify_poles.py names k = %d at %.2f%% as the smallest layer; it is "
                "k = %d at %.2f%%" % (_saidk, _saidpc, _worst, _pcs[_worst]))
        elif max(_pcs.values()) != 100.0:
            bad("verify_poles.py says 100%% is attained; the largest layer is %.2f%%"
                % max(_pcs.values()))
        else:
            print("   %d layers, %.2f%% at k = %d up to 100%% -- the header agrees"
                  % (len(_pcs), _pcs[_worst], _worst))
else:
    bad("scripts/verify_poles.py is missing")

# ---- 5r. the record-configuration axis layers, and the table that describes them --------
print()
print("5r. the record axis layers agree with the README table that states them")
_recs = {}
for _k in range(9, 24):
    for _nm in ('poles_rec_%d.npz' % _k, 'poles_sqrt2_%d.npz' % _k):
        _p = os.path.join(_pkg, 'data', _nm)
        if not os.path.exists(_p):
            continue
        _z = np.load(_p)
        _recs[_k] = int(len(_z['keep']))
        if _nm.startswith('poles_rec') and 'sha' not in _z.files:
            bad("data/%s is unsealed: run scripts/seal_record.py" % _nm)
        if _recs[_k] > PUBTAU[_k]:
            bad("the layer at k = %d is %d, above tau(k) = %d" % (_k, _recs[_k], PUBTAU[_k]))
if not _recs:
    print("   no record-configuration layer is shipped")
else:
    _rd = io.open(os.path.join(_pkg, 'README.md'), encoding='utf-8').read()
    _sec = _rd.find('The poles do not have to be a copy of the caps')
    if _sec < 0:
        bad("the package README no longer describes the record-configuration layers")
    else:
        # bound the window at the NEXT table: the per-k layer table below this one has
        # the same first three columns, and without the bound 5r reads fifteen rows out
        # of it against four shipped layers.
        _end = _rd.find('| dim | k | shell layer |', _sec)
        _win = _rd[_sec:_end if _end > _sec else _sec + 6000]
        _rows = re.findall(r'^\| (\d+) \| (\d+) \| ([0-9 \u00a0\u202f]+) \| '
                           r'([0-9 \u00a0\u202f]+) \| ([0-9 \u00a0\u202f]+) \|',
                           _win, re.M)
        _said = {}
        for _d, _kk, _before, _now, _tau in _rows:
            _said[int(_kk)] = (int(re.sub(r'[^0-9]', '', _now)),
                               int(re.sub(r'[^0-9]', '', _tau)))
        if set(_said) != set(_recs):
            bad("the README table describes k = %s; the shipped layers are k = %s"
                % (sorted(_said), sorted(_recs)))
        for _k in sorted(set(_said) & set(_recs)):
            if _said[_k][0] != _recs[_k]:
                bad("the README says dimension %d has %d poles; the shipped layer has %d"
                    % (72 + _k, _said[_k][0], _recs[_k]))
            elif _said[_k][1] != PUBTAU[_k]:
                bad("the README quotes tau(%d) = %d; the table says %d"
                    % (_k, _said[_k][1], PUBTAU[_k]))
        # the two axis spaces where the discriminant forbids a frame: 18 and 22, recomputed
        _noframe = []
        for _k in (18, 22):
            _W = np.load(os.path.join(_pkg, 'data', 'lam%d_W.npy' % _k)).astype(np.int64)
            _W = _W[:, np.abs(_W).sum(0) > 0]
            _r = int(np.linalg.matrix_rank(_W.astype(np.float64)))
            if _r % 2 == 0:
                _noframe.append(_k)
        if not all(('k = %d' % _k) in _rd[_sec:_sec + 6000] or
                   ('| **3** |' in _rd[_sec:_sec + 6000]) for _k in _noframe):
            bad("the README no longer names the axis spaces with no frame")
        print("   %d layers, %d described in the README table, tau(k) quoted correctly"
              % (len(_recs), len(_said)))

# The list of dimensions whose axis layer reaches tau(k) exactly is stated in prose in two
# documents and guarded by nothing.  It has just gained dimension 83, which had no layer at
# all.  Derive it and parse it: a list like this is precisely what goes stale in silence.
_ceil = []
for _k in range(1, 24):
    _best = 0
    for _nm in ('poles_rec_%d.npz' % _k, 'poles_sqrt2_%d.npz' % _k, 'poles_rot_%d.npz' % _k,
                'poles_rot_k%d.npz' % _k):
        _p = os.path.join(_pkg, 'data', _nm)
        if os.path.exists(_p):
            _best = max(_best, int(len(np.load(_p)['keep'])))
    for _nm in ('poles_%d.npy' % _k, 'poles_k%d.npy' % _k):
        _p = os.path.join(_pkg, 'data', _nm)
        if os.path.exists(_p):
            _best = max(_best, 2 * len(np.load(_p)))
    if _best and _k in PUBTAU and _best == PUBTAU[_k]:
        _ceil.append(72 + _k)


def _spanlist(_xs):
    """'74-82, 86, 88' from a sorted list of integers, en-dash for the ranges"""
    _out, _i = [], 0
    while _i < len(_xs):
        _j = _i
        while _j + 1 < len(_xs) and _xs[_j + 1] == _xs[_j] + 1:
            _j += 1
        _out.append(str(_xs[_i]) if _j == _i else
                    '%d%s%d' % (_xs[_i], chr(0x2013), _xs[_j]))
        _i = _j + 1
    return ', '.join(_out)


_want = _spanlist(sorted(_ceil))
_seen, _here = 0, 0
for _rel in ('../notes/PAPER-WORKLOG.md', '../notes/SUBMISSION-all-dimensions.md'):
    _fp = os.path.join(HERE, _rel)
    if not os.path.exists(_fp):
        continue
    _here += 1
    _flat = ' '.join(io.open(_fp, encoding='utf-8').read().split())
    for _m in re.finditer(r'reach(?:es|ing) \u03c4\(k\) exactly in (?:dimensions )?'
                          r'([0-9,\u2013 ]+?)\. ', _flat + ' '):   # to the PERIOD: the list is
                          # "74-83, 86, 88", and stopping at the first comma read
                          # only the first range and reported its own bug

        _seen += 1
        _got = _m.group(1).strip().rstrip('.').strip()
        if _got != _want:
            bad("%s says the layer reaches tau(k) exactly in %s; it is %s"
                % (_rel, _got, _want))
if not _here:
    # Both documents are the author's working folder, not the deliverable.  In a clone of
    # this package alone the sentence does not exist, so there is nothing to be wrong --
    # and reporting a failure there made `python audit.py`, the first command the README
    # gives a reviewer, end in "1 FAILURES" on a clean checkout.
    print("   the outer working folder is not here, so the sentence it states is not "
          "checked; the layer is at tau(k) in %s" % _want)
elif not _seen:
    bad("no document states which dimensions reach tau(k) exactly, in the shape 5r reads")
else:
    print("   the layer is at tau(k) in %s; %d document(s) agree" % (_want, _seen))


# ---- 5q. the denominator digit-ranges quoted in prose ---------------------------------
# paper/README.md and the package README each say how many digits the rotations'
# denominators run to.  Both said "three to ..." while the shipped layers start at two.  Both
# were right when written; the search only began preferring a small denominator on
# 2026-08-29, and D = 19 at k = 7 falsified them silently.
print()
print("5q. the axis-layer figures quoted only in prose match the shipped layers")
_NUM = {'one': 1, 'two': 2, 'three': 3, 'four': 4, 'five': 5, 'six': 6, 'seven': 7,
        'eight': 8,
        'nine': 9, 'ten': 10, 'eleven': 11, 'twelve': 12, 'thirteen': 13, 'fourteen': 14,
        'fifteen': 15, 'sixteen': 16, 'seventeen': 17, 'eighteen': 18, 'nineteen': 19}
_dg = {}
for _k in list(range(3, 9)) + [_x for _x in range(9, 24) if _x != 11]:
    for _f in ('poles_rot_%d.npz' % _k, 'poles_rot_k%d.npz' % _k):
        _p = os.path.join(_pkg, 'data', _f)
        if os.path.exists(_p):
            _dg[_k] = len(str(int(np.load(_p)['D'][0])))
if not _dg:
    bad("no rotated axis layer to read a denominator from")
else:
    _lo, _hi = min(_dg.values()), max(_dg.values())
    _said = []
    for _rel in ('../notes/PAPER-WORKLOG.md',
                 'verifications/improved/dim73-95-gamma72-caps/README.md'):
        _fp = os.path.join(HERE, _rel)
        if not os.path.exists(_fp):
            continue
        _flat = ' '.join(io.open(_fp, encoding='utf-8').read().split())
        for _m in re.finditer(r'denominator[^.]{0,60}?with (\w+) to (\w+) digits', _flat):
            _a, _b = _NUM.get(_m.group(1)), _NUM.get(_m.group(2))
            if _a is None or _b is None:
                # A number-word this table does not know reads as None and then compares
                # unequal, which is the right verdict by accident.  Say what happened.
                bad("5q cannot read the number-word '%s' or '%s'"
                    % (_m.group(1), _m.group(2)))
                continue
            _said.append((os.path.basename(os.path.dirname(_fp)) + '/README.md', _a, _b))
    if not _said:
        bad("no document states a denominator digit range in the shape 5q reads")
    else:
        _wrong = [(w, a, b) for w, a, b in _said if (a, b) != (_lo, _hi)]
        if _wrong:
            for _w, _a, _b in _wrong:
                bad("%s says the denominators run to %s..%s digits; they run %d..%d"
                    % (_w, _a, _b, _lo, _hi))
        else:
            print("   %d layers, denominators of %d to %d digits; %d document(s) agree"
                  % (len(_dg), _lo, _hi, len(_said)))

    # The other axis-layer figure that lives only in prose: what the shell rule cost, i.e.
    # what the rotated layers gain over it across 73-95.  Three documents quote it; two are
    # generated, this is the check that keeps all three honest.
    _sh = _ro = 0
    for _k in range(1, 24):
        _s1 = 0
        for _f in ('poles_%d.npy' % _k, 'poles_k%d.npy' % _k):
            _p = os.path.join(_pkg, 'data', _f)
            if os.path.exists(_p):
                _s1 = 2 * len(np.load(_p))
        _r1 = _s1
        # the CURRENT layer, whichever of the four kinds it is: the documents this compares
        # against say in so many words that their counts are the current ones.
        for _f in ('poles_rot_%d.npz' % _k, 'poles_rot_k%d.npz' % _k,
                   'poles_rec_%d.npz' % _k, 'poles_sqrt2_%d.npz' % _k,
                   # k = 12 and 13 since 2026-08-30: a rotated copy of the Kappa cap set
                   'poles_kap_%d.npz' % _k):
            _p = os.path.join(_pkg, 'data', _f)
            if os.path.exists(_p):
                _r1 = max(_r1, len(np.load(_p)['keep']))
        _sh += _s1
        _ro += _r1
    _tot = _ro - _sh
    _quotes = []
    for _rel, _pat in (('README.md',
                        r'became the search space, and cost ([0-9\u202f ]+) spheres'),
                       ('../notes/PAPER-WORKLOG.md',
                        r'\*\*\+([0-9\u202f ]+) spheres across \w+ dimensions\*\*'),
                       ('../notes/SUBMISSION-all-dimensions.md',
                        r'are worth ([0-9\u202f ]+) spheres across')):
        _fp = os.path.join(HERE, _rel)
        if not os.path.exists(_fp):
            continue
        _m = re.search(_pat, ' '.join(io.open(_fp, encoding='utf-8').read().split()))
        if _m:
            _quotes.append((_rel, int(re.sub(r'[^0-9]', '', _m.group(1)))))
    if not _quotes:
        bad("no document states what the shell rule cost, in the shape 5q reads")
    else:
        _off = [(w, v) for w, v in _quotes if v != _tot]
        if _off:
            for _w, _v in _off:
                bad("%s says the axis layers are worth %d spheres over the shell rule; "
                    "the shipped layers give %d" % (_w, _v, _tot))
        else:
            print("   the axis layers are worth %d over the shell rule; %d document(s) agree"
                  % (_tot, len(_quotes)))

print()
print("5s. the READMEs quote the size run_all.py and this tree actually have")
# "19 scripts, about 7 minutes" against a table of 53 taking 28, and "~300 files / 37 MB"
# against 359 -- a fifth out, in the first command a reader runs and in the line that says
# how big the download is.  Both are measured here instead of read.
sys.path.insert(0, HERE)
try:
    from run_all import JOBS as _JOBS                                        # noqa: E402
except Exception as _e:                                                      # noqa: BLE001
    bad("run_all.py does not import: %s" % _e)
    _JOBS = None
if _JOBS is not None:
    _fast = sum(1 for _j in _JOBS if _j[4])
    _rmp = os.path.join(HERE, 'README.md')
    _rm = io.open(_rmp, encoding='utf-8').read()
    _q = re.findall(r'python run_all\.py\s*(--full)?\s*#\s*(\d+) scripts', _rm)
    if len(_q) != 2:
        bad("README.md no longer states the script counts of both run_all.py modes")
    else:
        _want = {'': _fast, '--full': len(_JOBS)}
        _off = [(f or 'fast', int(v)) for f, v in _q if int(v) != _want[f]]
        if _off:
            for _f, _v in _off:
                bad("README.md says the %s set is %d scripts; run_all.py has %d"
                    % (_f, _v, _want['' if _f == 'fast' else '--full']))
        else:
            print("   README.md gives %d fast and %d in all, which is the table"
                  % (_fast, len(_JOBS)))
    # the tree's own size, after its own .gitignore -- the number the root README quotes
    _pats = [_l.split('#')[0].strip().rstrip('/')
             for _l in io.open(os.path.join(HERE, '.gitignore'), encoding='utf-8')
             if _l.split('#')[0].strip()]
    _n = _b = 0
    for _dp, _dn, _fn in os.walk(HERE):
        _dn[:] = [_d for _d in _dn if _d not in ('__pycache__', '.git', 'selftest-tmp')]
        for _f in _fn:
            _p = os.path.join(_dp, _f)
            _rel = os.path.relpath(_p, HERE).replace(chr(92), '/')
            if any(fnmatch.fnmatch(_rel, _pt) or fnmatch.fnmatch(_rel, '*/' + _pt)
                   or fnmatch.fnmatch(_f, _pt) for _pt in _pats):
                continue
            _n += 1
            _b += os.path.getsize(_p)
    _mb = _b / 1e6
    if _rr is None:
        print("   the outer working folder is not here, so its size line is not checked; "
              "the tree is %d files, %.0f MB" % (_n, _mb))
    else:
        _m = re.search(r'the deliverable\*\*, ~?([0-9]+) files / ([0-9]+) MB', _rr)
        if not _m:
            bad("the root README no longer states the deliverable's file count and size")
        else:
            _sf, _sm = int(_m.group(1)), int(_m.group(2))
            _e1 = abs(_sf - _n) / float(_n)
            _e2 = abs(_sm - _mb) / _mb
            if _e1 > 0.05 or _e2 > 0.05:
                bad("the root README says %d files / %d MB; the tracked tree is %d / %.0f "
                    "(%.0f%% and %.0f%% out)" % (_sf, _sm, _n, _mb, 100 * _e1, 100 * _e2))
            else:
                print("   the root README's %d files / %d MB is the tracked %d / %.0f MB, "
                      "within 5%%" % (_sf, _sm, _n, _mb))

print()
print("6. every claim says what would have to change for it to improve")
_unclassified = [d for d, v, pkg, how in claims if status_for(d)[0] is None]
if _unclassified:
    bad("no status for dimension(s) %s -- classify them in common/status.py"
        % sorted(_unclassified))
else:
    _by = {}
    for d, v, pkg, how in claims:
        _by.setdefault(status_for(d)[0], []).append(d)
    for _st in sorted(_by):
        _ds = _by[_st]
        _rng = runs(_ds)
        print("   %-10s %2d dims  %-12s %s" % (_st, len(_ds), _rng, LEGEND[_st]))

# The counts quoted in the front matter drift when a dimension is added: CITATION.cff still
# said 44 and omitted 68-69 long after they were claimed.  Tie them to the claim list.
_cff = os.path.join(HERE, 'CITATION.cff')
if os.path.exists(_cff):
    _t = io.open(_cff, encoding='utf-8').read()
    _t = _t.replace('\r\n', '\n')
    if ('in %d\n  dimensions' % len(claims)) not in _t:
        bad("CITATION.cff does not say '%d dimensions'" % len(claims))
    else:
        print("   CITATION.cff agrees that there are %d" % len(claims))
    # ... and the COUNT being guarded is why the LIST drifted: the abstract still ended at
    # 73-95 after dimension 96 was claimed, because nothing compared the enumeration.  A
    # guarded number beside an unguarded list is the scoped-claim trap.
    _abs = ' '.join(_t.split())
    _abs = _abs[_abs.index('dimensions:'):] if 'dimensions:' in _abs else ''
    _abs = _abs[:_abs.index('.')] if '.' in _abs else _abs
    _listed = set()
    for _tok in re.findall(r'(\d+)\s*-\s*(\d+)|(\d+)', _abs):
        if _tok[2]:
            _listed.add(int(_tok[2]))
        else:
            _listed |= set(range(int(_tok[0]), int(_tok[1]) + 1))
    _want = {d for d, _v, _p, _h in claims}
    if _listed != _want:
        bad("CITATION.cff enumerates %s, the claims are %s (missing %s, extra %s)"
            % (runs(sorted(_listed)), runs(sorted(_want)),
               runs(sorted(_want - _listed)) or 'none',
               runs(sorted(_listed - _want)) or 'none'))
    else:
        print("   CITATION.cff enumerates exactly those %d dimensions" % len(_listed))
_rm = os.path.join(HERE, 'README.md')
if os.path.exists(_rm):
    _t = io.open(_rm, encoding='utf-8').read()
    if ('**%d dimensions**' % len(claims)) not in _t:
        bad("README.md does not say '**%d dimensions**'" % len(claims))
    else:
        print("   README.md agrees that there are %d" % len(claims))
# ... and the same count is written in prose all over the tree, where nothing tied it to
# anything.  Going from 46 to 47 needed fifteen edits in three separate sweeps, each found by
# grepping a different phrasing, which is the duplicated-constant problem at documentation
# scale.  Sweep for a claim-count that is not the claim count.  HISTORICAL lines are allowed,
# but only by being named here with a reason -- an allowlist that must be read, not a pattern
# that silently forgives.
# "claims" is also a VERB: "Dimension 83 claims none" and "dimension 73-80 claims take no
# poles" are not counts.  So require a determiner in front, or an unambiguous noun after.
_COUNT_RE = re.compile(
    r'(?:all|the|among the|of the)\s+(\d{2})\s+claims\b'
    r'|(?<![-\u2013\u2014\d.])(\d{2})\s+claims,'
    r'|counted among the\s+(\d{2})\b'
    r'|(?<![-\u2013\u2014\d.])(\d{2})\s+(?:improvements|dimensions improved)\b',
    re.I)
_COUNT_OK = {
    # dated entries in the knowledge base, describing the state at the time of writing
    ('KNOWLEDGE.md', '**`kissing_verifications` brought in line.**  46 claims'),
    ('KNOWLEDGE.md', '`recheck_all.py` recomputes all 46'),
}
_count_bad = []
for _dp, _dn, _fn in list(os.walk(HERE)) + _walked:
    _dn[:] = [d for d in _dn if d not in ('__pycache__', '.git', 'classes', 'state')]
    for _f in _fn:
        if not _f.endswith(('.md', '.py', '.cff', '.tex')) or _f == 'audit.py':
            continue                       # audit.py quotes the strings it is looking for
        _p = os.path.join(_dp, _f)
        _rel = os.path.relpath(_p, HERE).replace('\\', '/')
        for _i, _line in enumerate(io.open(_p, encoding='utf-8', errors='replace'), 1):
            for _m in _COUNT_RE.finditer(_line):
                _n = next(g for g in _m.groups() if g)
                if int(_n) == len(claims):
                    continue
                if any(_k in _line for _fk, _k in _COUNT_OK if _fk == _f):
                    continue
                _count_bad.append((_rel, _i, _line.strip()[:78]))
if _count_bad:
    for _rel, _i, _line in _count_bad:
        bad("%s:%d says a claim count that is not %d -- %s" % (_rel, _i, len(claims), _line))
else:
    print("   no file states a claim count other than %d" % len(claims))

byp = {}
for d, v, pkg, how in claims:
    byp.setdefault(pkg, []).append(d)
for pkg in sorted(byp):
    ds = byp[pkg]
    print("   %-38s %2d dims: %s" % (pkg, len(ds),
          runs(ds) if len(ds) > 2 else ds))

# RESULTS.md is WRITTEN from these claims (--write), so a divergence between the two would
# otherwise be repaired by overwriting the file rather than reported.  Compare them.
_rp = os.path.join(HERE, 'RESULTS.md')
# In --write-results mode RESULTS.md is the OUTPUT, so comparing against it would only
# report the differences that the write is about to remove -- and the write is guarded on
# there being no failures, so the comparison would prevent its own fix.
if os.path.exists(_rp) and '--write-results' not in sys.argv:
    _txt = io.open(_rp, encoding='utf-8').read()
    _pat = re.compile(u'^\| (\d+) \| ([0-9  ]+) \| \*\*([0-9  ]+)\*\*', re.M)
    _res = {int(m.group(1)): int(re.sub(r'[^0-9]', '', m.group(3))) for m in _pat.finditer(_txt)}
    _mine = {d: v for d, v, _p, _h in claims}
    _agree = True
    if set(_res) != set(_mine):
        bad("RESULTS.md lists dimensions %s, the audit %s"
            % (sorted(set(_res) - set(_mine)), sorted(set(_mine) - set(_res))))
        _agree = False
    for _d in sorted(set(_res) & set(_mine)):
        if _res[_d] != _mine[_d]:
            bad("dim %d: RESULTS.md says %d, the audit recomputes %d"
                % (_d, _res[_d], _mine[_d]))
            _agree = False
    # The summary line used to print whatever the comparison found, so a reader skimming the
    # tail of a failing run was told RESULTS.md agreed on the line after being told which
    # dimensions it did not agree on.
    if _agree:
        print("   RESULTS.md agrees with every recomputed claim (%d dimensions)" % len(_mine))
    else:
        print("   RESULTS.md does NOT agree with the recomputation; see the failures above")

print()
print("=" * 92)
print("TOTAL %d dimensions improved.  %s"
      % (len(claims), "ALL CHECKS PASSED" if not fail else "%d FAILURES" % fail))
print("=" * 92)

# ------------------------------------------------------------------ RESULTS.md
if '--write-results' in sys.argv and not fail:
    def gr(x):
        return "{:,}".format(x).replace(",", "\u202f")

    PKGDOC = {
        'dim25-lens-heads': '1006 heads in the lens of a minimal vector, no removal shared, plus one non-lattice equator point',
        'dim26-27-iota-triangles': 'the coset triangle of three norm-6 vectors on every triangle of directions, the side chosen per head; in 26 the second side is the involution image of the first',
        'dim38-leech-large-codimension': 'the Leech cap construction at codimension 14',
        'dim39-ers-constant-weight': 'Edel-Rains-Sloane with the 2026 constant-weight codes',
        'dim49-63-p48-caps': 'the cap construction over P_48, with explicit classes',
        'dim62-63-ers-chain': 'the Edel-Rains-Sloane chain (n, 15, 2), levels 1-2 built here',
        'dim68-69-gamma72-cross-sections': 'k-point moment LP with the Gram exhibited in Gamma_72',
        'dim70-71-gamma72-cross-sections': 'k-point moment LP on cross-sections of Gamma_72',
        'dim73-95-gamma72-caps': 'the cap construction over Gamma_72, with zero-sum triples',
        'dim96-ers-takeover': 'the Edel-Rains-Sloane chain (96, 24, 6, 1)',
    }
    # Every package but one lives under improved/.  dim96-ers-takeover is under closed/,
    # because it was written to close dimension 96 against the cap construction -- which it
    # still does; the claim it carries is Edel-Rains-Sloane's, not this paper's.
    PKGTREE = {'dim96-ers-takeover': 'closed'}
    out = []
    out.append("# Results")
    out.append("")
    out.append("**Generated by `audit.py --write-results`. Do not edit by hand.**")
    out.append("")
    out.append("%d dimensions, all of them lower bounds on the kissing number tau(n) that improve\non the best previously published value."
               % len(claims))
    out.append("")
    out.append("| dim | previously published | **this work** | factor | status | construction |")
    out.append("|---|---|---|---|---|---|")
    for d, v, pkg, how in claims:
        prev, _ = floor_for(d)
        st, _reason = status_for(d)
        out.append("| %d | %s | **%s** | %.2f | `%s` | [%s](verifications/%s/%s/) |"
                   % (d, gr(prev), gr(v), v / prev, st or '?', PKGDOC[pkg],
                      PKGTREE.get(pkg, 'improved'), pkg))
    out.append("")
    out.append("## Status: what would have to change")
    out.append("")
    out.append("A status describes THIS PROJECT'S MECHANISM in that dimension, not tau(n).  No"
               " matching upper\nbound is known in any dimension here, so no value below is"
               " known to be optimal --\n`exhausted` means the idea is finished, not the"
               " dimension.")
    out.append("")
    out.append("| status | dims | what it means | why |")
    out.append("|---|---|---|---|")
    _by = {}
    for d, v, pkg, how in claims:
        _by.setdefault(status_for(d), []).append(d)      # group by (status, REASON)
    for (st, reason) in sorted(_by, key=lambda k: (k[0], min(_by[k]))):
        ds = _by[(st, reason)]
        rng = runs(ds)
        out.append("| `%s` | %s (%d) | %s | %s |" % (st, rng, len(ds), LEGEND[st], reason))
    out.append("")
    # These two counts moved when dimensions 62 and 63 left the cap construction; they are
    # derived from the grouping above rather than written out, so they cannot go stale.
    _ncls = sum(len(v) for k, v in _by.items() if k[0] == 'classes')
    out.append("So the %d `classes` dimensions all improve together when a better family of"
               " pairwise-disjoint\nclasses is found -- which is where compute pays -- and the"
               " other %d need a new idea or a\nrefreshed table, not more search."
               % (_ncls, len(claims) - _ncls))
    out.append("")
    out.append("## How each number is established")
    out.append("")
    out.append("| dim(s) | package | in this audit |")
    out.append("|---|---|---|")
    seen = []
    for d, v, pkg, how in claims:
        key = (pkg, how)
        if key in seen:
            continue
        seen.append(key)
        ds = sorted(x[0] for x in claims if x[2] == pkg and x[3] == how)
        rng = ("%d-%d" % (ds[0], ds[-1])) if len(ds) > 2 else ", ".join(map(str, ds))
        out.append("| %s | `%s` | %s |" % (rng, pkg, how))
    out.append("")
    out.append("## Where the previously published figures come from")
    out.append("")
    out.append("Every figure in the second column is produced by `common/published.py`, never")
    out.append("typed by hand. Three regimes:")
    out.append("")
    out.append("| range | source |")
    out.append("|---|---|")
    out.append("| 25-39 | Cohn's table directly |")
    out.append("| 49-63 | 52 416 000 + tau(k): the direct sum P_48 (+) best k-dimensional "
               "lattice. Cohn's table has no entry between 48 and 72. |")
    out.append("| 70-71 | 331 737 984, inherited by monotonicity from dimension 64 "
               "(Edel-Rains-Sloane 1998). No table has any entry for 65-71. |")
    out.append("| 73-95 | 6 218 175 600 + tau(k): the direct sum Gamma_72 (+) best "
               "k-dimensional lattice. |")
    out.append("")
    out.append("## Not claimed, and why")
    out.append("")
    out.append("| dim(s) | why |")
    out.append("|---|---|")
    out.append("| 28, 29, 30, 31 | the cap model reproduces Cohn's table *exactly* "
               "there; no slack. See `verifications/improved/dim73-95-gamma72-caps/scripts/"
               "calibrate.py`. (Dimension 26 was here until 2026-09-08, when the two-triangle "
               "configuration beat the table by 1082.) |")
    out.append("| 32-37, 40-43 | exactly at the Edel-Rains-Sloane value. "
               "`verifications/closed/dim32-44-ers-audit/` |")
    out.append("| 44, 45 | withdrawn: beaten by Sun-Wang, arXiv:2607.20359v3. "
               "`verifications/superseded/dim44-45-p48-cross-sections/` |")
    out.append("| 46, 47 | already in the literature. %s is Boyvalenkov-Cherkashin, "
               "Results Math. 80 (2025), Paper No. 3, equation (3) of the preprint "
               "arXiv:2312.05121; %s is "
               "Sun-Wang, arXiv:2607.20359, Table 3. Both are recovered here by a "
               "different route and neither is claimed. "
               "`verifications/recovered/dim46-47-p48-cross-sections/` |"
               % (gr(23766960), gr(12309600)))
    out.append("| 48 | 52 416 000 is P_48 itself. |")
    out.append("| 64-67 | the Gamma_72 cross-sections fall below the dimension-64 floor by "
               "k = 3. `verifications/closed/dim69-gamma72-three-point/`  (68 and 69 ARE "
               "claimed, by the k-point LP with the Gram exhibited, not by that bound.) |")
    out.append("| 72 | 6 218 175 600 is Gamma_72 itself. |")
    out.append("| 97 and above | **not evaluated here.** Monotonicity gives at least "
               "dimension 96's 12 886 999 232, and nothing in this repository claims "
               "anything there. The Edel-Rains-Sloane chain is the thing to evaluate, as it "
               "was at 96. |")
    out.append("| 9-19 | every published record is checked MAXIMAL -- no sphere can be "
               "added to any of them, and the margin of each is recorded. "
               "`verifications/closed/dim09-19-record-maximality/` |")
    out.append("| 17-21 | the Cohn-Li mechanism is at its exact ceiling: provably exhausted "
               "in 18, 20, 21, and dimension 19 closed at exactly Ho's 11948. "
               "`verifications/closed/dim17-23-cohn-li-mechanism/` |")
    _EXPOSURE_NOTE = (
        "**An exposure, at dimensions 93-95.** These three ARE claimed, above; this is "
        "about how much room they have. Every Edel-Rains-Sloane input is a construction, "
        "so its totals are LOWER bounds: one above Gamma_72 settles a dimension, one below "
        "it settles nothing. At 93-95 the chain reaches only 2.68e9, 2.68e9 and 4.83e9, "
        "which proves nothing by being small, so the package reports the threshold instead. "
        "Dimension 95 stands as long as %s, against %s from the best construction in hand, "
        "on the chain %s. `verifications/closed/dim96-ers-takeover/`" % _d95_threshold())
    out.append("| 22, 23 | Lambda_22 and Lambda_23 are maximal spherical codes, and Leech "
               "level sets cap at exactly 93150. "
               "`verifications/closed/dim22-23-maximal-cross-sections/` |")
    out.append("| 24 | 196560 is the Leech lattice, and it is optimal (Odlyzko-Sloane, "
               "Levenshtein). |")
    out.append("")
    out.append(_EXPOSURE_NOTE)
    out.append("")
    # No dimension in the table above may be one this repository CLAIMS.  It listed 64-69
    # while claiming 68 and 69, and briefly listed 93-95 while claiming all three, because
    # the exposure note was written as a row.  Parse the rows back and check.
    _nc_start = out.index("| dim(s) | why |")
    _nc = set()
    for _row in out[_nc_start + 2:]:
        if not _row.startswith('|'):
            break
        _spec = _row.split('|')[1]
        for _a, _b, _c in re.findall(r'(\d+)\s*-\s*(\d+)|(\d+)', _spec):
            if _c:
                _nc.add(int(_c))
            else:
                _nc |= set(range(int(_a), int(_b) + 1))
    _clash = sorted(_nc & {d for d, _v, _p, _h in claims})
    if _clash:
        bad("RESULTS.md's 'Not claimed' table lists dimensions that ARE claimed: %s"
            % runs(_clash))
    else:
        print("   'Not claimed' lists no dimension that is claimed")
    io.open(os.path.join(HERE, 'RESULTS.md'), 'w',
            encoding='utf-8', newline='\n').write("\n".join(out) + "\n")
    print("wrote RESULTS.md")

sys.exit(1 if fail else 0)
