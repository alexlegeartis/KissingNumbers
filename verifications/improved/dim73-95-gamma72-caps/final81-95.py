#!/usr/bin/env python3
"""Dimensions 81-95: the cap construction over Gamma_72 with explicit classes.

Two inputs, both produced and verified in this directory:

  * 32 000 pairwise-disjoint classes of Gamma_72 minimal lines (47 150 230 lines), built as
    automorphism images of a 2106-line class and verified in exact integer arithmetic --
    every vector of norm 8, every class with maximum |off-diagonal inner product| <= 2, all
    lines distinct.  Caro-Wei could only guarantee 221 lines per class.
  * zero-sum triple partitions of kissing configurations of R^k (lamtriples.py).

The gain is  2 * sum_i |C_i| * (|Z_i| - 1), so a triple part is worth 4 per line and a pair
part 2 (KNOWLEDGE.md section 50).  Per k we take whichever is better:

  k = 19..23: triple partitions of the A_j cross-sections of the Leech, which are smaller
              than Lambda_k but whose triples still beat pairs on the full Lambda_k
              (factor 61756 against 46575 at k = 23);
  k =  9..18: antipodal pairs on Lambda_k, since there the A_j cross-section has shrunk too
              far -- the honest fallback until triple partitions of the true Lambda_k are
              computed (section 53).

    python final.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'scripts'))
import capsha                                                      # noqa: E402
N = 6218175600
# Dimension 80 is READ from final73-80.py's own output, not repeated here: as a constant it
# had drifted 98 low, missing the axis layer that verify_poles.py certifies.  Run
# final73-80.py first; the fallback is only so this script still runs alone.
# The fallback is a FLOOR, not a copy: it is allowed to be behind, never ahead, and the
# assert below checks that on every run where the real value is available -- which is the
# only way a fallback that is almost never used can be kept honest.  It sat 32 low for
# some time with nothing comparing the two.
D80_FLOOR = 6218849196
_d80 = os.path.join(HERE, 'rows73-80.json')
if os.path.exists(_d80):
    import json as _json
    D80 = int(_json.load(open(_d80))['dim80'])
    assert D80 >= D80_FLOOR, (
        "the dimension-80 fallback %d is ABOVE final73-80.py's %d, so running this script "
        "without rows73-80.json would over-claim" % (D80_FLOOR, D80))
else:
    print("rows73-80.json is missing: using the floor %d for dimension 80, which may be "
          "behind. Run final73-80.py first." % D80_FLOOR)
    D80 = D80_FLOOR

# best known kissing numbers, Cohn's table verified 2026-08-20 (poles)
TAU = {9: 306, 10: 510, 11: 604, 12: 841, 13: 1154, 14: 1932, 15: 2564, 16: 4320,
       17: 5730, 18: 7654, 19: 11948, 20: 19448, 21: 29768, 22: 49896, 23: 93150}
# antipodal (lattice) kissing numbers, kept only to check that the triple scheme really does
# beat antipodal pairs on the full Lambda_k
ANTIP_FULL = {9: 272, 10: 336, 11: 438, 12: 756, 13: 918, 14: 1422, 15: 2340, 16: 4320,
              17: 5346, 18: 7398, 19: 10668, 20: 17400, 21: 27720, 22: 49896, 23: 93150}
# Zero-sum triple partitions are READ FROM THE SHIPPED CERTIFICATES, not hardcoded.  For each
# k the better of two configurations is taken: the best LATTICE (parts_<k>.npz over
# lam<k>_W.npy) and the RECORD kissing configuration of R^k (rec_<k>_parts.npz over
# rec_<k>_W.npy, from Cohn's dimensions1-24.txt).  The record wins exactly where it splits
# perfectly -- k = 9 (306 = 3*102), k = 10 (510 = 3*170), k = 14 (1932 = 3*644) -- and loses
# where it is a layered non-lattice packing with a thin triple structure.  Everything is
# re-checked exactly by scripts/verify_parts.py.
DATA = os.path.join(HERE, 'data')
TRI, SRC, CANDS = {}, {}, {}
for _k in range(9, 24):
    _best = None
    _cands = [('lattice', 'parts_%d.npz' % _k), ('record', 'rec_%d_parts.npz' % _k),
              ('kappa', 'kap%d_parts.npz' % _k)]
    if _k == 11:
        _cands.append(('record(Z[sqrt2])', 'rec11_parts.npz'))
    for _src, _pf in _cands:
        _f = os.path.join(DATA, _pf)
        if not os.path.exists(_f):
            continue
        _z = np.load(_f)
        if _src.startswith('record(Z'):
            _n = len(np.load(os.path.join(DATA, 'rec11_exact.npz'))['A'])
        else:
            _wf = {'lattice': 'lambda9_W.npy' if _k == 9 else 'lam%d_W.npy' % _k,
                   'record': 'rec_%d_W.npy' % _k,
                   'kappa': 'kap%d_W.npy' % _k}[_src]
            _n = len(np.load(os.path.join(DATA, _wf)))
        _T, _P = len(_z['triples']), len(_z['pairs'])
        CANDS.setdefault(_k, []).append((_T, _P, _n - 3 * _T - 2 * _P, _src))
        if _best is None or 2 * _T + _P > 2 * _best[0] + _best[1]:
            _best = (_T, _P, _n - 3 * _T - 2 * _P, _src)
    if _best is None:
        raise SystemExit("no certificate for k = %d" % _k)
    TRI[_k] = _best[:3]
    SRC[_k] = _best[3]

# ---------------------------------------------------------------------------------------
# CLASS FAMILIES.  A family is a set of PAIRWISE-DISJOINT classes of Gamma_72 minimal lines,
# a class being a set of lines with pairwise |<u,v>| <= 2.  Dimension 72+k consumes the T+P
# largest members of ONE family, and different dimensions may use different families -- each
# dimension's bound is an independent construction -- so the rule applied here is
#
#     gain(k) = max over VERIFIED families F with |F| >= T+P of
#                   4 * S_F(T) + 2 * (S_F(T+P) - S_F(T)),   S_F(j) = sum of F's j largest,
#
# stated once and applied uniformly, exactly as for the two P_48 families in dimensions
# 49-63.  A family with fewer than T+P classes cannot serve that dimension and is skipped;
# families are never mixed, since disjointness ACROSS families is not claimed.
FAMILIES = []
for _nm, _fn in (('images(Aut Gamma_72)', 'class_sizes_32000.npy'),
                 ('greedy(spread pool)', 'class_sizes_gpu.npy')):
    _fp = os.path.join(os.path.join(HERE, 'data'), _fn)
    if os.path.exists(_fp):
        _s = np.sort(np.asarray(np.load(_fp)).ravel())[::-1].astype(np.int64)
        FAMILIES.append((_nm, _s, np.concatenate([[0], np.cumsum(_s)])))
assert FAMILIES, "no class-size certificate found"
NCLS = max(len(_s) for _nm, _s, _c in FAMILIES)


def _poles_of(k, src):
    """the axis layer belonging to THIS cap source, or 0"""
    if src == 'kappa':
        names = ('poles_kap_%d.npz' % k, 'poles_kaprec_%d.npz' % k)
    else:
        names = ('poles_rec_%d.npz' % k, 'poles_sqrt2_%d.npz' % k, 'poles_rot_%d.npz' % k)
    got = []
    for fn in names:
        p = os.path.join(DATA, fn)
        if os.path.exists(p):
            got.append(int(len(np.load(p)['keep'])))
    pf = os.path.join(DATA, 'poles_%d.npy' % k)
    if os.path.exists(pf) and src != 'kappa':
        got.append(2 * len(np.load(pf)))
    return max(got) if got else 0


def best_gain(T, P):
    """the rule above; returns (gain, name of the family that attains it)"""
    cands = [(4 * int(_c[T]) + 2 * int(_c[T + P] - _c[T]), _nm)
             for _nm, _s, _c in FAMILIES if len(_s) >= T + P]
    assert cands, "no verified family has %d classes" % (T + P)
    return max(cands)


print(__doc__)
for _nm, _s, _c in FAMILIES:
    print("family %-22s %6d classes, sizes %4d..%4d, %d lines in total"
          % (_nm, len(_s), _s.min(), _s.max(), _s.sum()))

print()
print("  dim  k  scheme                          parts   gain          this work         previous")
def _capset(k, src):
    """(W, c) for the cap set a source names, in the coordinates a layer indexes"""
    if src == 'kappa':
        wf, cf = 'kap%d_W.npy' % k, 'kap%d_c.npy' % k
    elif src.startswith('record(Z'):
        return None, None                      # k = 11 lives in Z[sqrt2]; verify11.py
    elif src.startswith('record'):
        wf, cf = 'rec_%d_W.npy' % k, 'rec_%d_c.npy' % k
    else:
        wf, cf = ('lambda9_W.npy' if k == 9 else 'lam%d_W.npy' % k), None
    p = os.path.join(DATA, wf)
    if not os.path.exists(p):
        return None, None
    W = np.load(p).astype(np.int64)
    c = (np.load(os.path.join(DATA, cf)).astype(np.int64) if cf
         else np.ones(W.shape[1], np.int64))
    keep = (np.abs(W).sum(0) > 0) | (c != 1)
    return W[:, keep], c[keep]


rows = []
TABLE = []
POLESRC = {}                   # which of the three layer kinds each k uses
run = D80                      # tau is non-decreasing: a dimension whose own construction
                               # falls short is carried by the one below it
for k in range(9, 24):
    if k in TRI:
        # THE SOURCE IS CHOSEN ON THE TOTAL, not on the partition alone: 2T + P and the
        # pole count are different quantities and a source can win one while losing the
        # other.  At k = 13 the Kappa section wins the partition by 8 units and loses 188
        # poles; it wins the total by 33 268, which is what makes it the right choice --
        # and a rule that compared only the partition could not have said so.
        # ties broken by name, so the choice does not depend on the order the
        # candidate files happened to be listed in
        _pick = max(CANDS[k], key=lambda c: (best_gain(c[0], c[1])[0]
                                             + _poles_of(k, c[3]), c[3]))
        if _pick[3] != SRC[k]:
            TRI[k], SRC[k] = _pick[:3], _pick[3]
        T, P, S = TRI[k]
        need = T + P + S
        assert need <= NCLS, (k, need, NCLS)
        gain, fam = best_gain(T, P)
        # The axis layer is a ROTATED copy of the cap directions where one has been found
        # (scripts/axis_rotate.py, 2026-08-27), and the earlier lattice-shell layer where it
        # has not.  The rotated ones are two to three orders of magnitude larger, because the
        # shell rule could only take vectors of the lattice itself, while a rotation is free
        # of it and only has to clear the 30-degree caps.  No count is quoted here: the one
        # that was (92 954) had never been a shipped value, and the shipped value has since
        # changed twice.  scripts/verify_poles.py prints them, and checks them.
        # THREE kinds of layer exist and the largest wins.  A rotated copy of the cap
        # directions cannot exceed |Z|, and at k = 19, 20 and 21 the caps are the minimal
        # vectors of Lambda_k while tau(k) is larger; scripts/axis_record.py takes the record
        # configuration of R^k as the source instead, which is what those three dimensions
        # use.  scripts/verify_record.py checks them, scripts/verify_poles.py the others.
        # A LAYER BELONGS TO A CAP SET, NOT TO A DIMENSION.  Every kind below stores the
        # poles as indices into W, or as a rotation of W, so a layer built against one W is
        # not merely stale against another -- it indexes the wrong vectors.  When the cap
        # set changed at k = 12 and 13 (Lambda_k to the Kappa section K_k) the old layers
        # stayed on disk and stayed correct for the set they were built from, and the
        # candidate list is therefore keyed on SRC[k] so that the two can never be mixed.
        if SRC[k] == 'kappa':
            _files = [('rotated(kappa)', 'poles_kap_%d.npz' % k),
                      ('record(kappa)', 'poles_kaprec_%d.npz' % k)]
        else:
            _files = [('record', 'poles_rec_%d.npz' % k),
                      # k = 11: a rational rotation of a configuration that lives in
                      # Z[sqrt2].  scripts/verify11.py checks it, in that ring.
                      ('sqrt2', 'poles_sqrt2_%d.npz' % k),
                      ('rotated', 'poles_rot_%d.npz' % k)]
        # ... and the file is checked to have been built from that W, not merely named
        # for it: a layer is indices into a cap set, and reading one against another set is
        # a well-typed way to get an answer to a question nobody asked.  capsha.py.
        _Wk, _ck = _capset(k, SRC[k])
        cand = {}
        for _kind, _fn in _files:
            _p = os.path.join(DATA, _fn)
            if os.path.exists(_p):
                _z = np.load(_p)
                assert _Wk is None or capsha.matches(_z, _Wk, _ck), (
                    'k=%d: %s was built from a cap set other than %s' % (k, _fn, SRC[k]))
                cand[_kind] = int(len(_z['keep']))
        pf = os.path.join(DATA, 'poles_%d.npy' % k)
        if os.path.exists(pf) and SRC[k] != 'kappa':
            cand['shell'] = 2 * len(np.load(pf))
        kind, npol = max(cand.items(), key=lambda t: t[1]) if cand else ('none', 0)
        assert not cand or npol == max(cand.values()), 'a smaller layer was chosen'
        POLESRC[k] = kind
        gain += npol
        how = "%5d T + %3d P (%s)%s" % (T, P, SRC[k],
                                        " +%d poles" % npol if npol else "")
        # The triple scheme is CHOSEN here rather than compared, so check that it really is
        # the better of the two.  The alternative is antipodal pairs on the full Lambda_k,
        # which needs lambda = tau_lat(k)/2 classes and admits tau(k) poles.  At k = 23 that
        # would need 46575 classes against the 32000 available, so it is not even feasible.
        lam_alt = ANTIP_FULL[k] // 2
        if lam_alt <= NCLS:
            gain_alt = best_gain(0, lam_alt)[0] + TAU[k]
            assert gain >= gain_alt, (
                "k=%d: antipodal pairs would give %d, better than the triples' %d"
                % (k, gain_alt, gain))
    new = max(N + gain, run)
    run = new
    old = N + TAU[k]
    rows.append((72 + k, old, new, gain))
    TABLE.append({'dim': 72 + k, 'k': k, 'W': SRC[k], 'scheme': how.strip(), 'gain': int(gain),
                  'value': int(new), 'previous': int(old), 'family': fam,
                  # T and P are what scripts/make_readme_table.py needs to rebuild
                  # the family and crossover tables of README.md from the shipped
                  # .npy files; without them it would have to re-derive the parts.
                  'T': int(T), 'P': int(P), 'poles': int(npol),
                  'poles_from': POLESRC.get(k, 'none')})
    print("   %2d %2d  %-30s %6d  +%-11d %-16d  %d"
          % (72 + k, k, how, need, new - old, new, old))

print()
prev = D80
for d, old, new, g in rows:
    assert new >= prev, (d, new, prev)
    prev = new
print("   monotone above dimension 80 (%d): True" % D80)
own = [d for d, old, new, g in rows if N + g > D80]
print("   dimensions with a construction of their own: %d-%d" % (min(own), max(own)))
print("   dimension 81 now clears dimension 80 by %d on its own." % (N + rows[0][3] - D80))
print()
print("   Previously this package claimed +60112 ... +20493000 from Caro-Wei's 221 and")
print("   antipodal pairs, and dimensions 81-87 were carried by dimension 80.")
print("   Dimension 95 is now a factor %.2f on that." % (rows[-1][3] / 20493000))

# The README table is GENERATED from this, so it cannot drift from the driver.
import json
json.dump(TABLE, open(os.path.join(HERE, 'rows81-95.json'), 'w'), indent=1)
print()
print("   wrote rows81-95.json (%d rows); scripts/make_readme_table.py turns it into the "
      "markdown" % len(TABLE))
print("   table in README.md, so the two cannot disagree.")
