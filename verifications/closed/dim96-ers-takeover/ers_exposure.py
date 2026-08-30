#!/usr/bin/env python3
"""The FULL Edel-Rains-Sloane chain, evaluated in EVERY dimension this repository claims,
and the threshold that would overturn each claim.

    python ers_exposure.py

ers_sweep.py asks only about LEVEL 0, A(n_0, ceil(n_0/4)), which is a complete kissing
configuration on its own.  ers96.py asks the full question, but only at n = 93..96.  This
script asks it everywhere, because level 0 is not the whole construction:

    N(n) = sum_nu A(n, n_nu, n_nu) * A(n_nu, ceil(n_nu/4)),  n >= n_0 >= 4 n_1 >= 16 n_2 ...

and the level-1 term is what carries dimension 64 from 2^28 = 268 435 456 to ERS's published
331 737 984.  A claim that clears level 0 has NOT thereby cleared ERS.

WHICH WAY EACH INEQUALITY POINTS -- the trap this file exists to avoid.  Every ERS input is a
construction, so everything computed here is a LOWER bound on ERS.  A lower bound BELOW a
claim proves nothing.  So for each claimed dimension the script prints two things:

  reaches   the best ERS total assemblable from PUBLISHED codes.  A claim must exceed this
            or it is not an improvement at all.
  exposure  hold the chain fixed and ask how large the level-1 SUPPORT code A(n,w,w) would
            have to be for ERS to overtake the claim.  Divided by the best such code in hand,
            that ratio is the factor by which somebody would have to beat the best published
            construction to sink the claim.  Small ratio = fragile claim.

WHERE THE LEVEL-1 SUPPORT CODES COME FROM.  No table covers A(n,w,w) for n > 64 at these
weights.  Two sources, both published:

  ERS's own       A(64,16,16) >= 30828  (extended cyclic code)   } Edel-Rains-Sloane 1998,
                  A(80,16,16) >= 143780 (four orbits, L_2(79))   } and monotone in n
  the grid map    n = w*q, index the cells of a w x q grid and let each support take one cell
                  per row; two supports meet in w - d_H, so a q-ary [w,k,>=ceil(w/2)]_q code
                  gives A(w*q,w,w) >= q^k, and A(n,w,w) is monotone in n.  The quaternary
                  codes are codetables.de; only ones already cited in this package, or
                  obtained from them by SHORTENING ([m,k,d] -> [m-1,k-1,d]), are used.

RESULT.  Every claim clears what ERS reaches.  The exposure ratios are NOT uniform, and that
is the finding: dimension 68 sits far below every other claim, because its margin over ERS is
only 9 per cent and its level-1 weight 16 is the one weight at which ERS themselves built a
strong code.  It is the most fragile claim in the repository and this is where that is
written down.
"""
import os, sys, re
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
DELIV = os.path.join(HERE, '..', '..', '..')
sys.path.insert(0, os.path.join(HERE, '..', 'dim32-44-ers-audit'))
sys.path.insert(0, os.path.join(DELIV, 'common'))
import PIPELINE as P

# ---- level 0: A(n_0, ceil(n_0/4)) from codetables.de, read 2026-08-23 -------------------
TOP_K = {48: 24, 52: 21, 56: 23, 60: 25, 64: 28, 68: 25, 72: 28,
         80: 30, 84: 28, 88: 30, 92: 30, 96: 33}
DIRECT = {61: 25, 62: 26, 63: 27, 69: 25, 70: 26, 71: 27, 93: 31, 94: 31, 95: 32}

def sign0(m):
    """A(m, ceil(m/4)): the level-zero sign code, the largest term of the chain.

    TOP_K and DIRECT are LINEAR codes from codetables.de and cover only 48..96, so this
    returned 1 for every m below 48 -- silently dropping the chain's largest term.  Measured
    2026-08-26, that made best_chain return 0.78, 0.81, 0.48 and 0.29 of the totals Edel-
    Rains-Sloane themselves published at 32, 36, 40 and 44, so audit 4b was gating the low
    range against a floor missing its dominant term.  An understated rival flatters every
    claim measured against it, which is the direction that matters.

    So fall back to the repository's binary code table, which covers n = 0..512 and admits
    NONLINEAR codes -- level-zero sign vectors need only a binary code of length n_0 with
    d >= ceil(n_0/4), so nonlinear ones count, and at m = 65, 67, 68 they beat the best
    linear code by 4x, 1.5x and 1.2x.

    This can only RAISE the rival.  Measured before changing it: nothing in 49..96 moves,
    because in every claimed dimension there a better chain top is already available -- the
    paper's section 8.2 figures (dimension 68 at 1.47, dimension 95 at 4.06) are unchanged.
    Below 48 it is the whole difference, and best_chain now meets or exceeds all six of the
    published ERS totals it can reach (1.255, 1.104, 1.073, 1.000 at 32, 36, 40, 44; 1.000
    at 64 and 80), which is what a faithful evaluator with better tables must do.
    """
    best = DIRECT.get(m, 0)
    d = -((-m) // 4)
    top = 4 * d
    if top in TOP_K:
        best = max(best, TOP_K[top] - (top - m))
    linear = (1 << best) if best > 0 else 1
    try:
        return max(linear, P.Abin(m, d))
    except KeyError:
        return linear

# ---- the level-nu support codes A(n,w,w): supports pairwise meeting in <= w/2 -----------
Q4 = {24: 9, 23: 8, 22: 7, 21: 6, 20: 5, 19: 4, 18: 3, 17: 2,      # from [24,9,12]_4
      16: 7, 15: 6, 14: 5, 13: 4, 12: 3, 11: 2, 10: 1}             # from [16,7,8]_4
ERS_CW = {(64, 16): 30828, (80, 16): 143780}

def packing3(n):
    """maximum number of triples on n points meeting pairwise in <= 1 (Schonheim/Hanani)"""
    v = (n * ((n - 1) // 2)) // 3
    return v - 1 if n % 6 == 5 else v


def packing4(n):
    """maximum number of 4-subsets of an n-set meeting pairwise in <= 2.

    C(n,3)/4 exactly when an S(3,4,n) exists, i.e. n = 2 or 4 mod 6 (Hanani 1960); otherwise
    a LOWER bound by monotonicity from the largest such n.  Never the Johnson bound, which
    points the wrong way for what this script needs."""
    if n % 6 in (2, 4):
        return comb(n, 3) // 4, 'S(3,4,%d), Hanani' % n
    m = max(k for k in range(4, n + 1) if k % 6 in (2, 4))
    return comb(m, 3) // 4, 'S(3,4,%d) lengthened' % m


def _is_prime_power(q):
    """exactly one distinct prime divides q"""
    ps = set()
    m = q
    p = 2
    while p * p <= m:
        while m % p == 0:
            ps.add(p); m //= p
        p += 1
    if m > 1:
        ps.add(m)
    return len(ps) == 1


PRIMEPOWERS = [q for q in range(2, 256) if _is_prime_power(q)]


def gridmap(n, w):
    """the q-ary grid map at ANY weight, not just the quaternary table above.

    Split n >= w*q into a w x q grid and let each support take one cell per row.  Two
    supports meet in w minus the Hamming distance of their choice words, so a q-ary code
    of length w and distance ceil(w/2) gives supports meeting in <= floor(w/2), which is
    what A(n,w,w) asks.  Reed-Solomon supplies an MDS code whenever w <= q+1.

    Without this, acw fell back on published.Acw below weight 10 and returned 41664 for
    A(96,6,6) where the grid map gives 16^4 = 65536 -- the 763904 by which best_chain(96)
    disagreed with ers96.py, which had this map hardcoded.  An understated rival makes
    every margin in the paper's section 8.2 look wider than it is, so the error ran in
    the flattering direction and had to be removed."""
    best, why = 0, ''
    d = -((-w) // 2)                       # the distance the overlap condition needs
    for q in PRIMEPOWERS:
        if w * q > n or w > q + 1:
            continue
        v = q ** (w - d + 1)               # MDS: k = w - d + 1
        if v > best:
            best, why = v, 'grid map, RS [%d,%d,%d]_%d' % (w, w - d + 1, d, q)
    return best, why

def acw(n, w):
    if w == 1:
        return n, 'singletons'
    if w == 2:
        return comb(n, 2), 'pairs'
    if w == 3:
        return packing3(n), 'maximum partial STS'
    if w == 4:
        return packing4(n)
    best, why = 1, 'nothing published'
    if w in Q4 and 4 * w <= n:
        best, why = 4 ** Q4[w], 'grid map [%d,%d,%d]_4' % (w, Q4[w], -((-w) // 2))
    for (m, ww), v in ERS_CW.items():
        if ww == w and m <= n and v > best:
            best, why = v, 'ERS 1998 A(%d,%d,%d)' % (m, w, w)
    try:
        v = P.Acw(n, w)
        if v > best:
            best, why = v, P.Acw_src(n, w)[:34]
    except KeyError:
        pass                    # not tabulated at these parameters; the grid map below covers it
    g, gw = gridmap(n, w)
    if g > best:
        best, why = g, gw
    return best, why

def sd(w):
    return -((-w) // 4)

def abin(w, d):
    """A(w,d) from the published tables, or the trivial 1.

    The 1 is never reached: measured over ranked(), 115 051 calls and no fallback,
    because every level below the top has w <= n_0/4 and those lengths are all
    tabulated.  Every KeyError P.Abin raises in such a run comes from sign0, which
    answers with the linear bound rather than with 1.  Kept because A(w,d) >= 1 holds
    unconditionally, so should the tables ever shrink this UNDERSTATES the rival --
    the direction that flatters a claim, not the one that invents it.
    """
    try:
        return P.Abin(w, d)
    except KeyError:
        return 1

def term(n, w, n0):
    """the contribution of one level.  The TOP of the chain carries the sign code
    A(n_0, ceil(n_0/4)); its support code A(n, n_0, n_0) is 1 whenever two weight-n_0 words in
    length n are forced to meet in more than n_0/2, i.e. whenever 3*n_0 > 2*n, and is counted
    otherwise -- omitting it would understate the rival."""
    if w == n0:
        mult = 1 if 3 * n0 > 2 * n else acw(n, n0)[0]
        return mult * sign0(n0)
    return acw(n, w)[0] * abin(w, sd(w))

def best_chain(n):
    bv, bc = -1, None
    for ch in P.chains(n):
        if ch[0] > n:
            continue
        v = sum(term(n, w, ch[0]) for w in ch)
        if v > bv:
            bv, bc = v, ch
    return bv, bc


def claims():
    """{dimension: claimed value} straight out of RESULTS.md, the source of truth"""
    out = {}
    for line in open(os.path.join(DELIV, 'RESULTS.md'), encoding='utf-8'):
        m = re.match(r'\| (\d+) \| ([\d\s  ]+) \| \*\*([\d\s  ]+)\*\* \|', line)
        if m:
            out[int(m.group(1))] = int(re.sub(r'\D', '', m.group(3)))
    return out


def ranking():
    """every claimed dimension, as (dim, claim, reaches, chain, w1, have, need, factor, why,
    is_ers).

    `reaches` is what the full chain assembles from published codes -- a LOWER bound on ERS,
    since every input is a construction.  `need` is how large the level-1 support code would
    have to be for ERS to overtake the claim, and `factor` is that against the best code in
    hand.  `is_ers` marks the dimensions that ARE this construction, whose exposure is 1.00 by
    definition and says nothing.

    main() prints from this and paper/factcheck.py checks the paper against it, so the ranking
    has ONE implementation and the guard cannot disagree with the thing it guards."""
    CLAIM = claims()
    out = []
    for n in sorted(d for d in CLAIM if d >= 49):
        v, ch = best_chain(n)
        c, n0 = CLAIM[n], ch[0]
        w1 = ch[1] if len(ch) > 1 else None
        have, why = acw(n, w1) if w1 else (0, '')
        if w1 and abin(w1, sd(w1)) > 1:
            rest = sum(term(n, w, n0) for w in ch if w != w1)
            need = -((-(c - rest)) // abin(w1, sd(w1)))
            ratio = need / float(have) if have else float('inf')
        else:
            need, ratio = None, float('inf')
        # 62, 63 and 96 ARE the chain, so 'how much room does it have over the chain'
        # is not a question about them -- their ratio is 1.00 by construction and would
        # otherwise rank them as the most exposed claims in the repository.
        out.append((n, c, v, ch, w1, have, need, ratio, why, n in (62, 63, 96)))
    return out


def ranked(rows=None):
    """the rankable dimensions, most exposed first.

    Dimensions that ARE this construction are dropped here, in the one place that decision is
    made, as are those with no level-1 term to threaten."""
    rows = ranking() if rows is None else rows
    return sorted((r for r in rows if r[6] is not None and not r[9]), key=lambda r: r[7])


def main():
    """the report, and the gate: exits 1 if ERS reaches or beats any claim"""
    # ers96.py hardcodes the grid map it needs; this module derives it.  Two
    # implementations of one table drift, and this one did: before gridmap() existed
    # acw returned 41664 for A(96,6,6) against ers96.py's 65536, so best_chain(96) came
    # out 763904 below the value the paper quotes -- understating the rival, which is the
    # direction that flatters every margin below.  Tie them together.
    ERS96_GRID = {(96, 24): 4 ** 9, (96, 6): 16 ** 4,
                  (95, 23): 4 ** 8, (94, 23): 4 ** 8, (93, 23): 4 ** 8}
    for (nn, ww), vv in sorted(ERS96_GRID.items()):
        got = acw(nn, ww)[0]
        assert got >= vv, ('acw(%d,%d) = %d is below the value ers96.py uses, %d'
                           % (nn, ww, got, vv))
    assert best_chain(96)[0] == 12886999232, best_chain(96)

    CLAIM = claims()

    print(__doc__)
    print("=" * 106)
    print("  dim   this repo claims    ERS reaches      factor  chain             level-1 support code")
    print("=" * 106)
    rows = ranking()
    fail = [r[0] for r in rows if r[2] > r[1]]   # a rival construction ABOVE the claim
    for (n, c, v, ch, w1, have, need, ratio, why, IS_ERS) in rows:
        tail = ("A(n,%d,%d) >= %-8d %s" % (w1, w1, have, why)) if w1 else ""
        print("  %3d   %-17d  %-14d  %6.3f  %-17s %s" % (n, c, v, c / v, str(ch), tail))

    print()
    print("=" * 106)
    print("EXPOSURE: how far the best published level-1 support code would have to be beaten")
    print("=" * 106)
    print("  dim   claim            level-1   A(n,w,w) needed   best in hand   FACTOR   Johnson allows")
    for (n, c, v, ch, w1, have, need, ratio, why, ise) in sorted(rows, key=lambda r: r[7]):
        if need is None:
            continue
        jb = P.johnson_cw(n, 2 * ((w1 + 1) // 2), w1)
        print("  %3d   %-15d  w = %-4d  %-16d  %-13d  %6.2f   %-12d %s"
              % (n, c, w1, need, have, ratio, jb,
                 "(this claim IS the ERS chain)" if ise else ""))

    print()
    print("For calibration: ERS's own A(64,16,16) = 30828 beats the grid map's 4^7 = 16384 at the")
    print("same parameters by a factor 1.88.  A claim whose FACTOR above is near that is one good")
    print("constant-weight code away from being overtaken.")
    order = ranked(rows)
    if order:
        worst = order[0]
        print()
        print("MOST EXPOSED: dimension %d, factor %.2f" % (worst[0], worst[7]))
    print()
    if fail:
        print("*** ERS ALREADY REACHES OR BEATS THE CLAIM IN DIMENSIONS %s ***" % fail)
        raise SystemExit(1)
    print("No claim is below what the Edel-Rains-Sloane construction reaches from published codes;")
    print("dimensions 62, 63 and 96 EQUAL it, being that construction.  ALL CHECKS PASSED")


if __name__ == '__main__':
    main()
