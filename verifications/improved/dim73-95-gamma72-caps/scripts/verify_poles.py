#!/usr/bin/env python3
"""Exact verification of the axis layers, k = 9..23.

A pole is a point (0,a) with a a unit vector of R^k.  At cap level t = 2/3 the three
conditions are

  pole . equator  = 0                                          always
  pole . cap      = sqrt(1-t) <z,a> = <z,a>/sqrt3  <= 1/2   <=>  4<z,a>^2 <= 3 |z|^2 |a|^2
  pole . pole     = <a,a'> <= 1/2                           <=>  2<a,a'> <= |a||a'|

TWO KINDS OF LAYER ARE VERIFIED HERE, and each dimension uses whichever it ships.

  ROTATED (data/poles_rot_<k>.npz, since 2026-08-27).  The poles are a rotated copy of the
  cap directions themselves: a = Q w with Q = N/D rational and Q^T C Q = C, so the pole-pole
  condition is inherited from W being a 60-degree code and only the cap condition has to be
  checked.  Because D can have fifteen digits, 4X^2 would overflow; 3 m^2 D^2 is never a
  perfect square, so for the INTEGER X = (Nw)^T(Cz) the test "4X^2 <= 3 m^2 D^2" is exactly
  "2|X| <= isqrt(3 m^2 D^2)" -- one big-integer square root, then int64 throughout.
  These layers keep between 99% and 100% of the direction set; the caps are small
  and the smallest layer, at k = 18, is 99.68% of it.  audit.py check 5p reads that
  floor back out of this sentence and compares it with the shipped layers, because a
  measured number in a docstring is what goes stale with nothing noticing.
  The rotation must also PRESERVE THE SPAN of the directions:
  several of the sets are written in more coordinates than they span, and an isometry of
  those coordinates is not an isometry of the axis space.  That is checked here on a basis
  of the span, which settles it because Q is linear.

  RECORD (data/poles_rec_<k>.npz, since 2026-08-29).  At k = 19, 20 and 21 the poles are a
  rotated copy of the RECORD configuration of R^k rather than of the cap directions, which a
  copy of the caps cannot reach because tau(k) is larger than |Z| there.  Those three k are
  SKIPPED below and checked by scripts/verify_record.py instead; their rotated copies are
  still shipped and still correct, but they are no longer what the dimension uses, so
  checking their counts against rows81-95.json here would report a disagreement that is not
  one.  Dimension 83's layer is checked by scripts/verify11.py, in Z[sqrt2].

  KAPPA (data/poles_kap_<k>.npz, since 2026-08-30).  At k = 12 and 13 the cap directions
  are no longer Lambda_k but the Kappa sections K_12 and K_13 (scripts/kappa_caps.py), which
  are larger -- 756 against 648 and 918 against 906 -- and split perfectly into zero-sum
  triples.  A layer is indices into W together with a rotation OF W, so a layer built on
  Lambda_k says nothing about K_k: the old poles_rot_12.npz and poles_rec_13.npz were not
  merely stale against the new cap set, they indexed vectors that are no longer there, and
  they were deleted rather than left to be checked against the wrong configuration.  These
  two layers are rotated copies and are checked by exactly the ROTATED rule below.

  SHELL (data/poles_<k>.npy, the earlier rule).  The poles are taken of norm M = 2m in the
  same lattice, where |a-w|^2 = M + m - 2<a,w> being a lattice norm forces |<a,w>| <= M/2 = m
  against a threshold of sqrt(3Mm)/2 = 1.2247 m.  Nothing is assumed: every inequality is
  checked.  These are what the rotated layers replaced, and they are still verified so the
  comparison in the report is against something checked rather than remembered.

Both +a and -a are poles, so distinct directions need the two-sided form.  k = 11 has no
layer at all: dimension 83 uses the 604-point record realised over Z[sqrt2], whose Gram is
irrational, and neither rule is set up for it."""
import json
import os
import sys
from fractions import Fraction as F
from math import gcd, isqrt

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import capsha                                                      # noqa: E402

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
DATA = os.path.join(HERE, 'data')
print(__doc__)

ROWS = {r['k']: r for r in json.load(open(os.path.join(HERE, 'rows81-95.json')))}
# tau(k) is READ, not retyped: this file carried tau(20) = 17400 and tau(21) = 27720 until
# 2026-08-27, which are Lambda_20 and Lambda_21, both 2048 below the best known.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'common'))
from published import TAU                                # noqa: E402


def span_basis(W):
    """indices of rows of W spanning its row space, by fraction-free elimination

    The gcd step is not cosmetic: without it the entries double at every round and the
    24-coordinate direction sets produce numbers thousands of digits long."""
    rows_, piv, idx = [], [], []
    A = W.shape[1]
    for i in range(len(W)):
        r = [int(x) for x in W[i]]
        for p, b in zip(piv, rows_):
            if r[p]:
                r = [b[p] * x - r[p] * y for x, y in zip(r, b)]
                g = 0
                for x in r:
                    g = gcd(g, x)
                if g > 1:
                    r = [x // g for x in r]
        p = next((j for j in range(A) if r[j]), None)
        if p is None:
            continue
        piv.append(p)
        rows_.append(r)
        idx.append(i)
        if len(idx) == A:
            break
    return idx


def _solve_rows(B, xs):
    """Solve y B = x over Q for each x, and return None for any x not in the row space.

    One elimination on the augmented system, in Fractions.  A cofactor expansion here would
    be exact and unrunnable: these direction sets live in 24 coordinates."""
    k, A = len(B), len(B[0])
    M = [[F(int(B[r][col])) for r in range(k)] for col in range(A)]      # A x k
    rhs = [[F(int(x[col])) for x in xs] for col in range(A)]             # A x len(xs)
    piv = []
    r = 0
    for col in range(k):
        p = next((i2 for i2 in range(r, A) if M[i2][col]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        rhs[r], rhs[p] = rhs[p], rhs[r]
        pv = M[r][col]
        M[r] = [v / pv for v in M[r]]
        rhs[r] = [v / pv for v in rhs[r]]
        for i2 in range(A):
            if i2 != r and M[i2][col]:
                f = M[i2][col]
                M[i2] = [a - f * b for a, b in zip(M[i2], M[r])]
                rhs[i2] = [a - f * b for a, b in zip(rhs[i2], rhs[r])]
        piv.append(col)
        r += 1
    if r != k:
        return None                       # B was not independent -- caller's bug
    for i2 in range(r, A):                # the rows below carry the consistency conditions
        for t in range(len(xs)):
            if rhs[i2][t] != 0:
                return None
    return True


def poles_stay_in_the_axis_space(W, N):
    """Does the rotation N map the span of W into itself?

    The claim for dimension 72+k puts the cap directions AND the poles in one R^k.  Several
    direction sets are written in more coordinates than they span -- Lambda_21's is 27720
    vectors in Lambda_24's 24 coordinates spanning 21 -- and a rotation of those coordinates
    need not preserve the span.  That is what went wrong on 2026-08-27 at k = 13, 15, 17, 18,
    21, 22 and 23: every inner product was right and the configuration needed R^96 rather
    than R^93.  Q is linear, so checking a basis of the span settles it."""
    idx = span_basis(W)
    B = [[int(v) for v in W[i]] for i in idx]
    xs = [[int(v) for v in (W[i].astype(object) @ N.T)] for i in idx]     # D * Q w_i
    return _solve_rows(B, xs) is True, len(idx)


def directions(k):
    """(W, c, m): the cap directions dimension 72+k uses, their metric, and their norm.

    The choice must be the one final81-95.py records, not one re-derived here: the two rules
    already disagree at k = 11, and an axis layer checked against the wrong W would verify
    nothing."""
    src = ROWS[k]['W']
    if src.startswith('record'):
        W = np.load(os.path.join(DATA, 'rec_%d_W.npy' % k)).astype(np.int64)
        c = np.load(os.path.join(DATA, 'rec_%d_c.npy' % k)).astype(np.int64)
    elif src == 'kappa':
        W = np.load(os.path.join(DATA, 'kap%d_W.npy' % k)).astype(np.int64)
        c = np.load(os.path.join(DATA, 'kap%d_c.npy' % k)).astype(np.int64)
    else:
        f = 'lambda9_W.npy' if k == 9 else 'lam%d_W.npy' % k
        W = np.load(os.path.join(DATA, f)).astype(np.int64)
        c = np.ones(W.shape[1], np.int64)
    keep = (np.abs(W).sum(0) > 0) | (c != 1)
    W, c = W[:, keep], c[keep]
    m = int((W[0] * c) @ W[0])
    assert (np.einsum('ij,ij->i', W * c, W) == m).all(), 'k=%d: W is not of one norm' % k
    return W, c, m


ok = True
tot = 0
print("   k   |W|   poles  kind     cap condition                    pole-pole      tau(k)")
SUPERSEDED = []
for k in range(9, 24):
    # The layer file follows the CAP SET, not the dimension: see KAPPA above.
    if ROWS[k]['W'] == 'kappa':
        rot = os.path.join(DATA, 'poles_kap_%d.npz' % k)
        shell = rec = ''
    else:
        rot = os.path.join(DATA, 'poles_rot_%d.npz' % k)
        shell = os.path.join(DATA, 'poles_%d.npy' % k)
        rec = os.path.join(DATA, 'poles_rec_%d.npz' % k)
    if not os.path.exists(rot) and not (shell and os.path.exists(shell)):
        continue
    if rec and os.path.exists(rec):
        # A THIRD kind of layer exists at this k, built from the record configuration of R^k
        # rather than from the cap directions, and it is the one the dimension uses.  The
        # rotated copy below is still correct, but it is no longer what rows81-95.json
        # records, so checking its count against the row would report a disagreement that is
        # not one.  scripts/verify_record.py checks the layer that IS used.
        SUPERSEDED.append(k)
        continue
    W, c, m = directions(k)
    CZ = (W * c).T

    # W itself must be a 60-degree code, or the rotated copy inherits nothing.  Blocked --
    # the whole Gram is 18.5 GiB at k = 22 -- and computed in FLOAT, which is exact here and
    # only here: every entry of <w, w'> is a sum of at most 24 products of integers bounded
    # by 4 and 3, so it is an integer well under 2^53 and float64 represents each partial sum
    # exactly.  The bound is asserted rather than assumed, and the result is checked to be
    # integral.  numpy has no BLAS for int64, so this is the difference between seconds and
    # half an hour at k = 23; the CAP condition below cannot use it, because there the
    # products run to 10^18, and stays in int64.
    _bnd = int(np.abs(W).max()) ** 2 * int(np.abs(c).max()) * W.shape[1]
    assert _bnd < 2 ** 50, 'k=%d: the direction Gram is too large for exact float64' % k
    Wf, CZf = W.astype(np.float64), CZ.astype(np.float64)
    _BS = max(1, int(8e6 // max(1, len(W))))
    _mx = -float(m)
    for _s in range(0, len(W), _BS):
        _B = Wf[_s:_s + _BS] @ CZf
        assert (_B == np.rint(_B)).all(), 'the float Gram is not integral'
        for _t in range(len(_B)):
            _B[_t, _s + _t] = -m
        _mx = max(_mx, float(_B.max()))
    _mx = int(_mx)
    assert 2 * _mx <= m, 'k=%d: the cap directions are not a 60-degree code' % k

    if os.path.exists(rot):
        z = np.load(rot)
        # A layer is indices into W and a rotation OF W.  Where the file records which W it
        # was built from, that has to be the W being used here -- see scripts/capsha.py.
        assert capsha.matches(z, W, c),             'k=%d: %s was built from a different cap set' % (k, os.path.basename(rot))
        N, D, kp = z['N'].astype(np.int64), int(z['D'][0]), z['keep']
        # Q = N/D is an isometry of the form: N^T C N = D^2 C, in Python integers because
        # D^2 has thirty digits.
        No = N.astype(object)
        C = np.zeros((len(c), len(c)), dtype=object)
        for t in range(len(c)):
            C[t, t] = int(c[t])
        assert ((No.T @ C @ No) == (D * D) * C).all(), 'k=%d: N^T C N != D^2 C' % k
        lim = isqrt(3 * m * m * D * D)
        assert lim * lim < 3 * m * m * D * D < (lim + 1) ** 2, \
            'k=%d: 3 m^2 D^2 is a perfect square, so the integer test is not equivalent' % k
        assert 2 * m * D < 2 ** 62, 'k=%d: 2mD would overflow int64' % k
        # X = (N w_i)^T C w_j runs to 10^18, so this one cannot be done in float directly.
        # It CAN be done in two float passes, which matters: numpy has no BLAS for int64, and
        # the straight integer product is 2.1e11 operations at k = 23 -- half an hour against
        # a minute.  Split each entry of U = W[kp] N^T at bit 31, U = 2^31 A + R, and note
        #
        #   |A_ij| <= max|U| / 2^31,   |R_ij| < 2^31,   |C w_j| <= cmax * wmax  per coordinate,
        #
        # so every entry of A@CZ and of R@CZ is an integer bounded by 2^31 * cmax * wmax * k,
        # which the assert below pins under 2^50: float64 holds those exactly, and each
        # partial sum on the way is an integer of the same size, so the products are exact.
        # They are then recombined in int64, with |A@CZ| asserted small enough that the shift
        # cannot overflow.  Both asserts are checks, not assumptions -- if either failed the
        # script would stop rather than return a wrong answer.
        A = W[kp] @ N.T                       # D * (rotated pole directions)
        SH = 31
        Ahi, Alo = A >> SH, A - ((A >> SH) << SH)
        _cw = int(np.abs(CZ).max())
        _bd = (int(np.abs(Ahi).max()) + 2 ** SH) * _cw * W.shape[1]
        assert _bd < 2 ** 50, 'k=%d: the limb product is too large for exact float64' % k
        Hf, Lf, CZf2 = (Ahi.astype(np.float64), Alo.astype(np.float64),
                        CZ.astype(np.float64))
        worst = 0
        _B = max(1, int(8e6 // max(1, len(W))))
        for s in range(0, len(kp), _B):
            XH = Hf[s:s + _B] @ CZf2
            XL = Lf[s:s + _B] @ CZf2
            assert (XH == np.rint(XH)).all() and (XL == np.rint(XL)).all(), 'not integral'
            XHi = XH.astype(np.int64)
            assert int(np.abs(XHi).max()) < 2 ** 31, 'k=%d: the shift would overflow' % k
            X = (XHi << SH) + XL.astype(np.int64)
            worst = max(worst, int(np.abs(X).max()))
        cap_ok = 2 * worst <= lim
        # Pole-pole needs no scan of its own.  <Qw, Qw'> = <w, w'> because Q is an isometry
        # of the form, and the kept poles are indexed by a SUBSET of W, whose whole Gram was
        # just bounded above by m/2.  A second pass over the kept poles would re-derive that
        # from the same numbers; what is checked instead is the step that is not obvious --
        # that Q really is an isometry, which is the N^T C N = D^2 C above.
        _pm = _mx
        pp_ok = 2 * _pm <= m
        # AND the poles must be in the axis space.  Nothing above asks this, and at
        # k = 13, 15, 17, 18, 21, 22 and 23 the first rotations built for these layers
        # failed it: those direction sets are written in more coordinates than they span --
        # Lambda_21's is 27720 vectors in Lambda_24's 24 coordinates spanning 21 -- and a
        # rotation of the AMBIENT coordinates is not an isometry of R^k.  Every inner
        # product was correct and the configuration needed R^96 instead of R^93.
        _span_ok, _nb = poles_stay_in_the_axis_space(W, N)
        if not _span_ok:
            ok = False
            print('   k=%d: the rotation does not preserve the span of the directions, so '
                  'the poles are NOT in R^%d' % (k, k))
        npol, kind = len(kp), 'rotated'
        left = '2*%d <= %d' % (worst, lim)
        right = '2*%d <= %d' % (_pm, m)
    else:
        P = np.load(shell).astype(np.int64)
        Pg, Wg = P * c, W * c
        M = int((Pg[0] * P[0]).sum())
        assert (np.einsum('ij,ij->i', Pg, P) == M).all()
        a1 = 0
        for s in range(0, len(P), 512):
            B = Pg[s:s + 512] @ W.T
            a1 = max(a1, int((4 * B * B).max()))
        cap_ok = a1 <= 3 * M * m
        IP = Pg @ P.T
        np.fill_diagonal(IP, 0)
        a2 = int((4 * IP * IP).max())
        pp_ok = a2 <= M * M
        npol, kind = 2 * len(P), 'shell'
        left = '%d <= %d' % (a1, 3 * M * m)
        right = '%d <= %d' % (a2, M * M)

    good = cap_ok and pp_ok
    ok = ok and good
    tot += npol
    print("  %2d %6d %6d  %-8s %-32s %-14s %6d%s"
          % (k, len(W), npol, kind, left, right, TAU[k],
             '  AT THE CEILING' if npol == TAU[k] else ''))
    if not good:
        print("      *** k=%d FAILS a pole condition ***" % k)
    if ROWS[k]['poles'] != npol:
        ok = False
        print("      *** k=%d: rows81-95.json says %d poles, this layer has %d ***"
              % (k, ROWS[k]['poles'], npol))


if SUPERSEDED:
    print('  k = %s: the axis layer comes from the RECORD configuration of R^k, not from a'
          % ', '.join(str(k) for k in SUPERSEDED))
    print('  rotated copy of the caps; scripts/verify_record.py checks those, exactly.')


# ---- dimensions 75-80: the root systems, rebuilt here rather than imported ---------------
# The roots are their coordinate vectors in the simple-root basis and the Gram is the Cartan
# matrix, so E_6 and E_7 need no embedding in R^8.  The rotation is a G-isometry:
# N^T G N = D^2 G.  Built from the Cartan matrix by reflection closure, sorted, so the index
# set in the shipped file means the same thing here as it did there -- if the two orders ever
# disagreed the check below would fail rather than pass on the wrong object.
CARTAN = {3: [(0, 1), (1, 2)], 4: [(0, 1), (1, 2), (1, 3)],
          5: [(0, 1), (1, 2), (2, 3), (2, 4)],
          6: [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5)],
          7: [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6)],
          8: [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)]}
RTAU = {3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240}


def root_system(k):
    G = 2 * np.eye(k, dtype=np.int64)
    for i, j in CARTAN[k]:
        G[i, j] = G[j, i] = -1
    seen = set()
    for r in np.eye(k, dtype=np.int64):
        seen.add(tuple(int(x) for x in r))
        seen.add(tuple(-int(x) for x in r))
    frontier = list(seen)
    while frontier:
        new = []
        for r in frontier:
            rv = np.array(r, dtype=np.int64)
            for i in range(k):
                e = np.zeros(k, dtype=np.int64)
                e[i] = 1
                t = tuple(int(x) for x in rv - (rv @ G[i]) * e)
                if t not in seen:
                    seen.add(t)
                    new.append(t)
        frontier = new
    R = np.array(sorted(seen), dtype=np.int64)
    assert len(R) == RTAU[k], (k, len(R))
    assert (np.einsum('ij,jk,ik->i', R, G, R) == 2).all(), k
    return R, G


_any = False
for k in range(3, 9):
    rf = os.path.join(DATA, 'poles_rot_k%d.npz' % k)
    if not os.path.exists(rf):
        continue
    if not _any:
        print()
        print("   k   dim  roots  poles  cap condition                 pole-pole   tau(k)")
        _any = True
    R, G = root_system(k)
    off = R @ G @ R.T
    np.fill_diagonal(off, -2)
    _roff = int(off.max())
    assert 2 * _roff <= 2, 'k=%d: the roots are not a 60-degree code' % k
    z = np.load(rf)
    N, D, kp = z['N'].astype(np.int64), int(z['D'][0]), z['keep']
    No, Go = N.astype(object), G.astype(object)
    assert ((No.T @ Go @ No) == (D * D) * Go).all(), 'k=%d: N^T G N != D^2 G' % k
    # A root system has full rank in its own root coordinates, so this cannot fail here --
    # which is the point of asking: the same question at k = 13, 15, 17, 18, 21, 22 and 23,
    # where the coordinates outnumber the span, is what caught seven bad layers.
    assert poles_stay_in_the_axis_space(R, No)[0], 'k=%d: the poles are not in R^%d' % (k, k)
    m = 2
    lim = isqrt(3 * m * m * D * D)
    assert lim * lim < 3 * m * m * D * D < (lim + 1) ** 2, k
    X = (R[kp] @ N.T) @ (R @ G).T
    w = int(np.abs(X).max())
    cap_ok = 2 * w <= lim
    pp_ok = 2 * _roff <= m          # inherited: the kept poles index a subset of the roots
    ok = ok and cap_ok and pp_ok
    tot += len(kp)
    print("  %2d  %3d  %5d %6d  2*%-12d <= %-12d 2*%-3d <= %-3d %6d%s"
          % (k, 72 + k, len(R), len(kp), w, lim, _roff, m, RTAU[k],
             '  AT THE CEILING' if len(kp) == RTAU[k] else ''))
    if not (cap_ok and pp_ok):
        print("      *** k=%d FAILS a pole condition ***" % k)

print()
print("ALL AXIS LAYERS VERIFIED EXACTLY -- %d poles in total across the dimensions above"
      % tot if ok else "*** SOME AXIS LAYER FAILED ***")
sys.exit(0 if ok else 1)
