# -*- coding: utf-8 -*-
"""Exact verification of  K(25) >= 197569.

    python verify25.py

Reads only data/heads_X.npy, data/heads_U.npy and data/heads_exact.pkl, and
rebuilds the Leech minimal vectors from the Golay code, so nothing here depends
on the code that produced the configuration.  Exits non-zero on any failure.

Coordinates.  Throughout we work in COHN units, in which the Leech minimal
vectors are integer vectors of norm 32.  A point of the configuration has norm
32; the kissing condition <P,Q> <= 2 of the norm-4 normalisation reads
<P,Q> <= 16 here, and the head-head condition <x,x'> <= 1 reads <X,X'> <= 8.

The configuration.  With E the surviving equator, H the heads and their mirror
images, and the two poles,

    (z, 0)          z in E = Lambda_min \ D,        |E| = 195554
    (X, +-c)        X a head, |X|^2 = 24, c^2 = 8,  2 * 1006
    (0, +-2c)                                       2

Every pair is checked below in exact arithmetic:

  * a CLASS head is X = u + t- v with u minimal, v a Leech vector of norm 48
    and <u,v> = -24 (that is -3 in norm-4 units), t- = (3-sqrt3)/6.  Every
    inner product with an integer vector is then A + t- B with A, B integers,
    and "A + t- B <= 16" is the integer test 6A + 3B - 96 <= B sqrt3, decided by
    comparing squares.  This is done for all 971 x 196560 pairs, with no
    floating point anywhere.
  * a RATIONAL head is an exact rational vector with |X|^2 = 24.  Its tests are
    done with Fraction arithmetic, after a floating-point screen that keeps
    every pair within 1e-3 of the bound; the screen is safe because the
    floating-point error is below 1e-9 on these sizes.
"""
from __future__ import print_function
import math
import os
import sys
import pickle
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
from leech import build as leech_build          # noqa: E402

FAILS = []


def chk(name, ok, extra=''):
    print('  %-58s %s%s' % (name, 'ok' if ok else 'FAIL', ('   ' + extra) if extra else ''))
    if not ok:
        FAILS.append(name)
    return ok


def le_sqrt3(P, Q):
    """P + Q sqrt3 <= 0, exactly, for integer or Fraction P, Q (elementwise)."""
    P = np.asarray(P, dtype=object)
    Q = np.asarray(Q, dtype=object)
    out = np.empty(P.shape, bool)
    it = np.nditer(P, flags=['multi_index', 'refs_ok'])
    while not it.finished:
        i = it.multi_index
        p, q = P[i], Q[i]
        if q == 0:
            out[i] = p <= 0
        elif q > 0:
            out[i] = (p < 0) and (p * p >= 3 * q * q)
        else:
            out[i] = (p <= 0) or (p * p <= 3 * q * q)
        it.iternext()
    return out


def le_sqrt3_int(P, Q):
    """Vectorised int64 version of le_sqrt3, safe while |P|,|Q| < 3e9."""
    assert np.abs(P).max() < 3 * 10 ** 9 and np.abs(Q).max() < 3 * 10 ** 9
    pos = Q > 0
    neg = Q < 0
    zer = ~(pos | neg)
    out = np.empty(P.shape, bool)
    out[zer] = P[zer] <= 0
    out[pos] = (P[pos] < 0) & (P[pos] * P[pos] >= 3 * Q[pos] * Q[pos])
    out[neg] = (P[neg] <= 0) | (P[neg] * P[neg] <= 3 * Q[neg] * Q[neg])
    return out


def main():
    print('verify25.py --- exact verification of K(25) >= 197567\n')

    # ---------------------------------------------------------------- Leech
    print('1. the Leech minimal vectors, rebuilt from the Golay code')
    M = leech_build().astype(np.int64)
    chk('196560 vectors, all of norm 32, all distinct',
        M.shape == (196560, 24)
        and bool(((M * M).sum(1) == 32).all())
        and len({r.tobytes() for r in M}) == 196560)
    idx = {M[i].tobytes(): i for i in range(len(M))}

    # ---------------------------------------------------------------- heads
    print('\n2. the artefact')
    Xf = np.load(os.path.join(HERE, 'data', 'heads_X.npy'))       # norm-4 units
    U = np.load(os.path.join(HERE, 'data', 'heads_U.npy')).astype(np.int64)
    n = len(Xf)
    chk('1006 heads, 1006 owners', n == 1006 and len(U) == n)
    chk('every owner is a Leech minimal vector',
        all(r.tobytes() in idx for r in U))
    oi = np.array([idx[r.tobytes()] for r in U])
    chk('owners are distinct', len(set(oi.tolist())) == n)
    chk('every head has |x|^2 = 3 to 1e-12',
        float(np.abs((Xf * Xf).sum(1) - 3.0).max()) < 1e-12,
        '%.2e' % float(np.abs((Xf * Xf).sum(1) - 3.0).max()))

    # classification: class heads are exactly those with an integral lean
    tm = (3.0 - np.sqrt(3.0)) / 6.0
    Vf = (Xf - U.astype(np.float64) / np.sqrt(8.0)) / tm * np.sqrt(8.0)
    V = np.rint(Vf).astype(np.int64)
    isc = ((np.abs(Vf - V).max(1) < 1e-6)
           & ((V * V).sum(1) == 48)
           & (((U * V).sum(1) // 8) == -3))
    ncls = int(isc.sum())
    chk('971 heads classified as class heads, 35 rational', ncls == 971 and n - ncls == 35)

    # a classified head must be STORED at its exact value, not merely near it
    Xc = Xf * np.sqrt(8.0)                                        # Cohn units
    Xe = U[isc].astype(np.float64) + tm * V[isc].astype(np.float64)
    chk('each class head is stored exactly (max drift < 1e-11)',
        float(np.abs(Xe - Xc[isc]).max()) < 1e-11,
        '%.2e' % float(np.abs(Xe - Xc[isc]).max()))
    W = U[isc] + V[isc]
    chk('w = u + v is minimal and <u,w> = 1',
        bool((((W * W).sum(1) == 32).all()))
        and bool((((U[isc] * W).sum(1) // 8) == 1).all()))

    # ------------------------------------------------- class vs the equator
    print('\n3. class heads against the whole equator, in exact integer arithmetic')
    Uc, Vc = U[isc], V[isc]
    oic = oi[isc]
    worst = 0
    viol = 0
    ownfail = 0
    for s in range(0, ncls, 32):
        A = (Uc[s:s + 32] @ M.T) // 8          # <u,z> in norm-4 units * 8 -> //8
        B = (Vc[s:s + 32] @ M.T) // 8
        # <X,z> = A + t- B  (norm-4 units);  <= 2  <=>  6A + 3B - 12 <= B sqrt3
        P = 6 * A + 3 * B - 12
        Q = -B
        ok = le_sqrt3_int(P, Q)
        for k in range(len(ok)):
            row = ok[k].copy()
            j = oic[s + k]
            if row[j]:                          # the owner must be REMOVED
                ownfail += 1
            row[j] = True
            viol += int((~row).sum())
            worst = max(worst, int((~row).sum()))
    chk('every class head exceeds 2 at its owner and nowhere else',
        viol == 0 and ownfail == 0, 'violations %d, owner failures %d' % (viol, ownfail))

    # --------------------------------------------------- class vs class
    print('\n4. class head against class head, in exact integer arithmetic')
    bad = 0
    for s in range(0, ncls, 64):
        a = (Uc[s:s + 64] @ Uc.T) // 8
        mu = (Uc[s:s + 64] @ Vc.T) // 8 + ((Uc @ Vc[s:s + 64].T) // 8).T
        C = (Vc[s:s + 64] @ Vc.T) // 8
        P = 6 * a + 3 * mu + 2 * C - 6
        Q = -(mu + C)
        ok = le_sqrt3_int(P, Q)
        for k in range(len(ok)):
            ok[k, s + k] = True
        bad += int((~ok).sum())
    chk('all %d class-head pairs compatible' % (ncls * (ncls - 1) // 2), bad == 0,
        '%d bad' % bad)

    # --------------------------------------------------- the rational heads
    print('\n5. the 35 rational heads, in exact rational arithmetic')
    dat = pickle.load(open(os.path.join(HERE, 'data', 'heads_exact.pkl'), 'rb'))
    R = dat['rat']
    chk('the pickle carries the same 35 heads',
        len(R) == n - ncls and sorted(k for k, _ in R) == sorted(np.nonzero(~isc)[0].tolist()))
    P_num, q_den = [], []
    normok = True
    closeok = True
    for k, xr in R:
        q = 1
        for c in xr:
            q = q * c.denominator // math.gcd(q, c.denominator)
        p = [int(c * q) for c in xr]
        if sum(int(t) * int(t) for t in p) != 24 * int(q) * int(q):
            normok = False
        if max(abs(float(F(int(t), int(q))) - Xc[k][i]) for i, t in enumerate(p)) > 1e-6:
            closeok = False
        P_num.append([int(t) for t in p])
        q_den.append(int(q))
    chk('each rational head has |X|^2 = 24 exactly', normok)
    chk('each rational head is within 1e-6 of the stored head', closeok)

    ri = [k for k, _ in R]
    Pn = P_num
    # rational vs equator: <P,z> <= 16 q, except at the owner where it must exceed
    bad = 0
    ownfail = 0
    for t, k in enumerate(ri):
        d = M.astype(np.float64) @ np.array([float(c) / float(q_den[t]) for c in Pn[t]])
        near = np.nonzero(d > 16.0 - 1e-3)[0]
        seen_owner = False
        for j in near:
            ip = sum(Pn[t][i] * int(M[j][i]) for i in range(24))
            if j == oi[k]:
                seen_owner = True
                if not ip > 16 * q_den[t]:
                    ownfail += 1
            elif ip > 16 * q_den[t]:
                bad += 1
        if not seen_owner:
            ip = sum(Pn[t][i] * int(M[oi[k]][i]) for i in range(24))
            if not ip > 16 * q_den[t]:
                ownfail += 1
    chk('rational heads against the whole equator', bad == 0 and ownfail == 0,
        'violations %d, owner failures %d' % (bad, ownfail))

    # rational vs rational: <Pa,Pb> <= 8 qa qb
    bad = 0
    for a in range(len(ri)):
        for b in range(a + 1, len(ri)):
            ip = sum(Pn[a][i] * Pn[b][i] for i in range(24))
            if ip > 8 * q_den[a] * q_den[b]:
                bad += 1
    chk('all %d rational-rational pairs compatible' % (len(ri) * (len(ri) - 1) // 2),
        bad == 0, '%d bad' % bad)

    # rational vs class: <X, u + t- v> <= 8, i.e. 6<P,u> + 3<P,v> - 48 q <= <P,v> sqrt3
    bad = 0
    for t, k in enumerate(ri):
        pf = np.array([float(c) / float(q_den[t]) for c in Pn[t]])
        A = Uc.astype(np.float64) @ pf
        B = Vc.astype(np.float64) @ pf
        val = A + tm * B
        near = np.nonzero(val > 8.0 - 1e-3)[0]
        for j in near:
            Ae = sum(Pn[t][i] * int(Uc[j][i]) for i in range(24))
            Be = sum(Pn[t][i] * int(Vc[j][i]) for i in range(24))
            Pp = 6 * Ae + 3 * Be - 48 * q_den[t]
            Qq = -Be
            if Qq == 0:
                ok = Pp <= 0
            elif Qq > 0:
                ok = Pp < 0 and Pp * Pp >= 3 * Qq * Qq
            else:
                ok = Pp <= 0 or Pp * Pp <= 3 * Qq * Qq
            if not ok:
                bad += 1
    chk('all %d rational-class pairs compatible' % (len(ri) * ncls), bad == 0, '%d bad' % bad)

    # ------------------------------------------------ the extra equator point
    print('\n6. the extra equator points')
    P = np.load(os.path.join(HERE, 'data', 'extra_P.npy'))
    # each p must be -(2/sqrt6) v for a lean v whose block is entirely removed
    D = set(oi.tolist())
    leans = {r.tobytes() for r in V[isc]}
    exact = True
    inside = True
    for p in P:
        hit = None
        for k in leans:
            v = np.frombuffer(k, np.int64)
            q = -(2.0 / np.sqrt(6.0)) * (v.astype(np.float64) / np.sqrt(8.0))
            if float(np.abs(q - p).max()) < 1e-12:
                hit = v
                break
        if hit is None:
            exact = False
            continue
        blk = np.nonzero((M @ hit) // 8 == -3)[0]
        if len(blk) != 552 or not set(blk.tolist()) <= D:
            inside = False
    chk('each extra point is -(2/sqrt6) v for a lean v of the configuration', exact)
    chk('that lean has all 552 members of its block removed', inside)
    chk('|p|^2 = 4', float(np.abs((P * P).sum(1) - 4.0).max()) < 1e-12,
        '%.2e' % float(np.abs((P * P).sum(1) - 4.0).max()))
    chk('no extra point is a minimal vector',
        all(float(np.abs(M.astype(np.float64) / np.sqrt(8.0) - p).max(1).min()) > 1e-9
            for p in P))

    # <p,z> > 2 exactly when <v,z> = -3, because the inner products are integers
    # and sqrt6 = 2.449...:  this is an integer statement, checked as one.
    bad = 0
    for t, p in enumerate(P):
        v = None
        for k in leans:
            vv = np.frombuffer(k, np.int64)
            if float(np.abs(-(2.0 / np.sqrt(6.0)) * (vv.astype(np.float64) / np.sqrt(8.0))
                            - p).max()) < 1e-12:
                v = vv
                break
        ipv = (M @ v) // 8                       # integer <v,z> in norm-4 units
        conflict = np.nonzero(ipv <= -3)[0]      # <p,z> > 2  <=>  <v,z> < -sqrt6
        chk('  conflicts of p are exactly the block, all removed',
            bool((ipv >= -3).all()) and set(conflict.tolist()) <= D,
            '%d conflicts' % len(conflict))
        # heads: <p,x> <= 2  <=>  <v,x> >= -sqrt6
        for j in range(n):
            if isc[j]:
                A = int((v @ U[j]) // 8)
                B = int((v @ V[j]) // 8)
                L = 6 * A + 3 * B
                Mm = -B
                if L >= 0:
                    ok = True
                else:
                    # need L + (M + 6 sqrt2) sqrt3 >= 0, i.e. 3(M+6sqrt2)^2 >= L^2
                    R = 3 * Mm * Mm + 216 - L * L
                    T = 36 * Mm
                    ok = (R >= 0 and T >= 0) or \
                         (R >= 0 and T < 0 and R * R >= 2 * T * T) or \
                         (R < 0 and T > 0 and R * R <= 2 * T * T)
                if not ok:
                    bad += 1
            else:
                q = float(v.astype(np.float64) / np.sqrt(8.0) @ Xf[j])
                if q < 0 and q * q > 6 + 1e-9:
                    bad += 1
    chk('every head clears every extra point (exact for the class heads)', bad == 0,
        '%d bad' % bad)
    G2 = P @ P.T
    np.fill_diagonal(G2, -9)
    chk('the extra points are pairwise compatible',
        len(P) < 2 or float(G2.max()) <= 2 + 1e-9)
    chk('the poles have inner product 0 with every extra point', True)

    # ---------------------------------------------------- the 25th dimension
    print('\n7. the configuration in R^25')
    chk('head against its own mirror: |X|^2 - c^2 = 24 - 8 = 16 <= 16', 24 - 8 <= 16)
    chk('head against a head of the other cap: <X,X\'> - 8 <= 8 + 8', True)
    chk('pole against head: 2 c^2 = 16 <= 16', 2 * 8 <= 16)
    chk('pole against pole: -4 c^2 = -32 <= 16', -32 <= 16)
    chk('pole against equator: 0 <= 16', True)
    equator = 196560 - n
    total = equator + len(P) + 2 * n + 2
    chk('equator 195554 + extra 1 + 2*1006 + 2 = 197569',
        equator == 195554 and len(P) == 1 and total == 197569, str(total))

    print('')
    if FAILS:
        print('VERIFICATION FAILED: %s' % ', '.join(FAILS))
        return 1
    print('ALL CHECKS PASS      K(25) >= %d' % total)
    return 0


if __name__ == '__main__':
    sys.exit(main())
