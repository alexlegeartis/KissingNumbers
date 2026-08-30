#!/usr/bin/env python3
"""Exact verification of the axis layers built from the RECORD configuration of R^k.

    python verify_record.py

scripts/axis_rotate.py makes the poles a rotated copy of the CAP DIRECTIONS, so its layers
can never exceed |Z|.  At k = 19, 20 and 21 the caps are the minimal vectors of Lambda_k and
tau(k) is larger, so those layers were short by construction.  scripts/axis_record.py builds
the record configuration instead.  This file checks the result and assumes none of it.

WHAT A LAYER HAS TO SATISFY.  A pole is a point (0,a) of the sphere of R^(72+k).  At cap
level t = 2/3 it needs

    4 <a,z>^2 <= 3 |a|^2 |z|^2      against every cap direction z   (at least 30 degrees)
    <a,a'>    <= 1/2 |a| |a'|       against every other pole        (at least 60 degrees)

and the second is ONE-SIDED: an obtuse pair is legal and -a is a legal pole beside a.  The
record configuration contains complementary sign vectors, which sit at 180 degrees, so a
two-sided test would reject the very thing being verified.

WHAT IS REBUILT HERE.  Everything except the code words and the rotation:

  * the octads, read out of data/lam<k>_W.npy at k = 19 and 20, where the cap file is already
    written in k coordinates, and built from the Golay code at k = 21, where it is not;
  * the odd-sign base, 4 C(n,2) pair vectors and 128 sign patterns of odd weight per octad;
  * the extra vectors, the +-1 vectors of the shipped code words.

The shipped file carries only the code words, the frame M, and the rotation N/D -- three
things a search produced -- and every inequality below is checked in exact integer
arithmetic.

THE FRAME.  At k = 21 the cap directions are 27 720 vectors in 24 coordinates spanning 21, so
a configuration built in R^21 has to be carried in by a rational similarity M with
M^T M = mu I.  Then a pole is M Q t, and

    <w, M Q t> = <M^T w, Q t>,

so the caps are pulled back by M^T and the test is the one above in R^k.  The pole-pole
condition needs no scan: <M Q t, M Q t'> = mu <t, t'> because M^T M = mu I and N^T N = D^2 I,
both checked, so the 60-degree condition on the poles IS the 60-degree condition on the
source, which is checked directly.

WHY THE FRAME EXISTS ONLY AT SOME k.  The discriminant of the axis space -- the determinant
of any Gram of it, modulo rational squares -- is the product of the norms of any orthogonal
basis.  A frame of common norm mu therefore forces disc = mu^k.  At k = 18 and 22 the
discriminant is 3 and the rank is even, so mu^k is a square and NO frame exists: those axis
spaces cannot hold a rational copy of anything built in the standard R^k.  At k = 21 the
discriminant is 1, so mu must be a square -- which is why an orthogonal frame among the cap
directions, of norm 32, always stops one short, and why the frame shipped here is built from
norm-4 vectors instead.  The discriminant is recomputed below.
"""
import hashlib
import json
import os
import sys
from fractions import Fraction as F
from math import isqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', '..', 'common'))
from axis_record import (census, dual_code, halved, caps, octads_shortened,   # noqa: E402
                         odd_base, pair_vectors, vectors_of)
import axis_cohn                                                              # noqa: E402
from published import TAU                                                     # noqa: E402

PKG = os.path.join(HERE, '..')
DATA = os.path.join(PKG, 'data')
ROWS = {r['k']: r for r in json.load(open(os.path.join(PKG, 'rows81-95.json')))}


def squarefree(n):
    n = abs(int(n))
    out, d = 1, 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e % 2:
            out *= d
        d += 1
    return out * n


def discriminant(W):
    """the determinant of a Gram of the span, modulo rational squares"""
    rows, piv, basis = [], [], []
    for r in W:
        v = [F(int(x)) for x in r]
        for c, b in zip(piv, rows):
            if v[c]:
                f = v[c]
                v = [a - f * y for a, y in zip(v, b)]
        c = next((i for i, x in enumerate(v) if x), None)
        if c is None:
            continue
        basis.append([int(x) for x in r])
        rows.append([x / v[c] for x in v])
        piv.append(c)
    B = np.array(basis, dtype=object)
    G = (B @ B.T).tolist()
    n = len(G)
    A = [[F(int(G[i][j])) for j in range(n)] for i in range(n)]
    det = F(1)
    for col in range(n):
        p = next((r for r in range(col, n) if A[r][col]), None)
        if p is None:
            return 0, n
        if p != col:
            A[col], A[p] = A[p], A[col]
            det = -det
        det *= A[col][col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for r in range(col + 1, n):
            if A[r][col]:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    assert det.denominator == 1
    return squarefree(int(det)), n


def in_span(W, M):
    """every column of M lies in the row space of W, by exact rational elimination"""
    rows, piv = [], []
    for r in W:
        v = [F(int(x)) for x in r]
        for c, b in zip(piv, rows):
            if v[c]:
                f = v[c]
                v = [a - f * y for a, y in zip(v, b)]
        c = next((i for i, x in enumerate(v) if x), None)
        if c is None:
            continue
        rows.append([x / v[c] for x in v])
        piv.append(c)
        if len(piv) == W.shape[1]:
            break
    for col in range(M.shape[1]):
        v = [F(int(M[i][col])) for i in range(M.shape[0])]
        for c, b in zip(piv, rows):
            if v[c]:
                f = v[c]
                v = [a - f * y for a, y in zip(v, b)]
        if any(x for x in v):
            return False, col
    return True, len(piv)


def published_source(k, cs):
    """the record configuration of R^k as data/pole_src_<k>.npy, checked against the data set

    The k whose layers come from scripts/axis_cohn.py do not rebuild their source from the
    Golay code: it is a PUBLISHED configuration, given as exact integer coordinates by
    H. Cohn's dimensions1-24.txt.  The package ships it so that this file needs nothing
    outside itself, and re-reads the data set to compare whenever it can be found -- the
    shipped copy is then checked against its source rather than against itself.
    """
    X = np.load(os.path.join(DATA, 'pole_src_%d.npy' % k)).astype(np.int64)
    c = np.load(os.path.join(DATA, 'pole_src_%d_c.npy' % k)).astype(np.int64)
    assert (c == cs).all(), 'k=%d: pole_src_%d_c.npy disagrees with the layer file' % (k, k)
    assert X.shape[1] == k, 'k=%d: the shipped source has %d coordinates' % (k, X.shape[1])
    ds = axis_cohn.data_set()
    if ds:
        Y, cy, drop = axis_cohn.record(k, ds)
        assert Y is not None, "k=%d: the data set's block has no integral point" % k
        assert (cy == c).all(), 'k=%d: the data set gives a different metric' % k
        assert Y.shape == X.shape and (Y == X).all(), \
            'k=%d: data/pole_src_%d.npy is NOT the configuration in the data set' % (k, k)
        # A block is not always wholly integral: dimension 13 writes 48 of its 1154 points
        # with a 2*sqrt(3) in them and those are not usable here.  The count that was
        # dropped is part of what is checked, so that a shipped source silently missing
        # points would not pass as "matches the data set".
        assert len(X) + drop == TAU[k], \
            ('k=%d: the data set gives %d integral and %d other points, tau(k) = %d'
             % (k, len(X), drop, TAU[k]))
        print('   k=%d: the source matches dimensions1-24.txt exactly (%d points%s)'
              % (k, len(X),
                 '' if not drop else ', the other %d of tau(%d) = %d being irrational'
                 % (drop, k, TAU[k])))
    else:
        print('   k=%d: dimensions1-24.txt not found; the shipped source is checked but not'
              ' compared with the data set' % k)
    return X, 0, 0


def rebuild(k, code, ncols):
    """the record configuration of R^k, from the cap file or the Golay code"""
    if ncols == k:
        V = halved(k)
        oct_, parity, pairs = census(V)
        assert parity == 0, 'the shipped octads are not of even sign parity'
    else:
        oct_ = octads_shortened(k)
        pairs = pair_vectors(k)
    Z = odd_base(oct_, pairs, k)
    E = vectors_of(code, k)
    # the extras must be legal: their words lie in the dual of the octad code, and they are
    # pairwise far enough apart.  Both are re-derived, not taken from the file.
    D, rank, ker = dual_code(oct_, k)
    Dset = set(D)
    dmin = (k + 3) // 4
    for c in code:
        assert int(c) in Dset, 'a code word is not in the dual of the octad code'
    cs = sorted(int(c) for c in code)
    for i in range(len(cs)):
        for j in range(i + 1, len(cs)):
            assert bin(cs[i] ^ cs[j]).count('1') >= dmin, \
                'two code words are closer than distance %d' % dmin
    return np.concatenate([Z, E]), len(pairs), len(Z)


def is_code(T, cs=None):
    """the 60-degree condition on the source, one-sided, exact, blocked

    Every entry of <t,t'> is a sum of at most 24 products of integers bounded by 2 -- or by
    the source's metric weights, where it has one -- so it is an integer far below 2^53 and
    float64 represents it and every partial sum exactly.  The bound is asserted and the
    result checked to be integral, as elsewhere in this package.
    """
    c = np.ones(T.shape[1], dtype=np.int64) if cs is None else cs.astype(np.int64)
    _bnd = int(np.abs(T).max()) ** 2 * T.shape[1] * int(c.max())
    assert _bnd < 2 ** 50, 'the source Gram is too large for exact float64'
    norms = ((T * c) * T).sum(1)
    Tf = (T * c).astype(np.float64)
    Uf = T.astype(np.float64)
    B = max(1, int(8e6 // max(1, len(T))))
    worst, bad = -10 ** 9, 0
    for s in range(0, len(T), B):
        X = Tf[s:s + B] @ Uf.T
        assert (X == np.rint(X)).all(), 'the float Gram is not integral'
        Xi = X.astype(np.int64)
        for t in range(Xi.shape[0]):
            Xi[t, s + t] = -10 ** 9                 # a pole against itself is not a pair
        worst = max(worst, int(Xi.max()))
        na = norms[s:s + B][:, None]
        bad += int(((Xi > 0) & (4 * Xi * Xi > na * norms[None, :])).sum())
    return worst, bad


def check(k):
    """verify data/poles_rec_<k>.npz; (npol, ok)"""
    p = os.path.join(DATA, 'poles_rec_%d.npz' % k)
    if not os.path.exists(p):
        return None
    z = np.load(p)
    N, D, kp = z['N'].astype(np.int64), int(z['D'][0]), z['keep']
    code = [int(c) for c in z['code']]
    M = z['M'].astype(np.int64)
    # THE SOURCE'S OWN METRIC.  The configurations rebuilt from the Golay code live in the
    # standard R^k and carry no 'cs'.  The published ones do: Cohn's data set writes the
    # record configuration of R^15, R^17 and R^18 with inner-product coefficients
    # diag(1^14,2), diag(1^16,2) and diag(1^16,2,6), and at k = 18 that is the whole reason
    # a frame exists at all -- see the note on the discriminant above.
    cs = z['cs'].astype(np.int64) if 'cs' in z.files else None
    Cs = np.eye(M.shape[1], dtype=np.int64) if cs is None else np.diag(cs)
    W = caps(k)
    m = int(W[0] @ W[0])
    assert (np.einsum('ij,ij->i', W, W) == m).all(), 'the cap directions are not of one norm'
    ok = True

    # the frame is a similarity onto the axis space
    G = M.T @ M
    assert int(G[0, 0]) % int(Cs[0, 0]) == 0, 'M^T M is not a multiple of diag(c)'
    mu = int(G[0, 0]) // int(Cs[0, 0])
    assert (G == mu * Cs).all(), 'M^T M is not mu diag(c)'
    Zk = W @ M
    # Every cap lies in the span, so it is its own projection through the frame:
    # w = M diag(1/(mu c)) M^T w, hence sum_i (M^T w)_i^2 / c_i = mu <w,w>.  Cleared of
    # denominators that is the identity below, and for c = 1 it is the plain norm check it
    # replaces.  A frame whose image were SMALLER than the span would fail it.
    _L = 1
    for _v in np.diag(Cs):
        _L = _L * int(_v) // np.gcd(_L, int(_v))
    _w = np.array([_L // int(v) for v in np.diag(Cs)], dtype=np.int64)
    assert ((Zk.astype(object) ** 2) @ _w == mu * m * _L).all(), \
        'the pull-back changed a cap norm: the frame does not span the axis space'
    # AND the poles must be IN the axis space.  They lie in the image of M, so that image has
    # to be the span of the cap directions -- checked here rather than inferred, because a
    # rotation that carried the poles out of R^k is a mistake this package has already made,
    # with every inner product in the layer correct.
    _ok, _w = in_span(W, M)
    assert _ok, 'column %d of the frame is not in the span of the cap directions' % _w
    sf, rank = discriminant(W)
    assert rank == M.shape[1], 'the frame has %d columns for a span of %d' % (M.shape[1], rank)
    _dc = 1
    for _v in np.diag(Cs):
        _dc *= int(_v)
    assert squarefree(mu ** rank * _dc) == sf, \
        ('a frame of norm %d for a metric of discriminant %d is impossible in a space of '
         'discriminant %d' % (mu, squarefree(_dc), sf))

    # the rotation is an isometry of the source's own metric
    No, Co = N.astype(object), Cs.astype(object)
    assert ((No.T @ Co @ No) == (D * D) * Co).all(), 'N^T diag(c) N != D^2 diag(c)'

    # the source
    if cs is None:
        T, npair, nbase = rebuild(k, code, W.shape[1])
    else:
        T, npair, nbase = published_source(k, cs)
    # The source can be SHORT of tau(k) and still be legal: at k = 13 forty-eight of the
    # 1154 points of the record configuration are written with a 2*sqrt(3) in them and are
    # not rational, so 1106 are used.  It can never be LONGER -- that would mean a 60-degree
    # code beating the published record -- and a source shorter than it should be can only
    # cost poles, never invent them, which is why the exact split is checked against the
    # data set in published_source() and only bounded here.
    assert len(T) <= TAU[k], 'the source has %d vectors, tau(%d) = %d' % (len(T), k, TAU[k])
    # THE INDEX SET INDEXES INTO THIS.  The source is rebuilt here rather than shipped, which
    # makes the ORDER of the rebuilt rows part of the contract: a different order would make
    # `keep` select different poles, and the check below would then be checking something
    # else.  scripts/seal_record.py stamps the file with the digest of the source the search
    # actually used; compare it.
    assert 'sha' in z.files, 'the layer is unsealed: run scripts/seal_record.py'
    _sha = np.frombuffer(hashlib.sha256(np.ascontiguousarray(T, dtype=np.int64).tobytes())
                         .digest(), dtype=np.uint8)
    assert (_sha == z['sha']).all(),         'k=%d: the rebuilt source is not the one the index set was made for' % k
    worst, bad = is_code(T, cs)
    if bad:
        ok = False
        print('   k=%d: the source is NOT a 60-degree code, %d pairs break it' % (k, bad))

    # the cap condition on the kept poles, exactly.  X = (N t)^T (M^T w) is an integer and
    # the test is 4 X^2 <= 3 (mu m) |t|^2 D^2; 3 (mu m) |t|^2 D^2 is not a perfect square, so
    # it is 2|X| <= isqrt(...) with one big-integer root per norm and int64 after.
    A = T[kp] @ N.T
    norms = ((T * np.diag(Cs)) * T).sum(1)[kp]
    lim = {}
    for na in sorted(set(int(v) for v in norms)):
        L = 3 * mu * m * na * D * D
        r = isqrt(L)
        assert r * r < L < (r + 1) ** 2, 'k=%d: 3 mu m |t|^2 D^2 is a perfect square' % k
        lim[na] = r
    SH = 31
    Ahi, Alo = A >> SH, A - ((A >> SH) << SH)
    _bd = (int(np.abs(Ahi).max()) + 2 ** SH) * int(np.abs(Zk).max()) * len(N)
    assert _bd < 2 ** 50, 'k=%d: the limb product is too large for exact float64' % k
    Hf, Lf, Zf = Ahi.astype(np.float64), Alo.astype(np.float64), Zk.astype(np.float64)
    lo = np.array([lim[int(v)] for v in norms], dtype=np.int64)
    B = max(1, int(8e6 // max(1, len(Zk))))
    capbad, mx = 0, 0
    for s in range(0, len(kp), B):
        XH = Hf[s:s + B] @ Zf.T
        XL = Lf[s:s + B] @ Zf.T
        assert (XH == np.rint(XH)).all() and (XL == np.rint(XL)).all(), 'not integral'
        XHi = XH.astype(np.int64)
        assert int(np.abs(XHi).max()) < 2 ** 31, 'k=%d: the shift would overflow' % k
        X = (XHi << SH) + XL.astype(np.int64)
        w = np.abs(X).max(1)
        mx = max(mx, int(w.max()))
        # the comparison DOUBLES w, and w is already a recombined two-limb product; assert
        # the headroom rather than trust it, since an overflow here would silently turn a
        # pole inside a cap into a pole that passes.
        assert 2 * int(w.max()) < 2 ** 63 - 1, 'the doubled inner product would overflow'
        capbad += int((2 * w > lo[s:s + B]).sum())
    if capbad:
        ok = False
        print('   k=%d: %d kept poles are inside a cap' % (k, capbad))
    npol = len(kp)
    print('  %2d  %3d %7d %7d %7d %7d  %3d  %s'
          % (k, 72 + k, len(W), len(T), npol, TAU[k] - npol, mu, D))
    print('           cap condition, every kept pole against every cap:  2*%d <= %d%s'
          % (mx, int(lo.min()), '   AT THE CEILING' if npol == TAU[k] else ''))
    if ROWS[k]['poles'] != npol:
        ok = False
        print('      *** k=%d: rows81-95.json says %d poles, this layer has %d ***'
              % (k, ROWS[k]['poles'], npol))
    return npol, ok


def main():
    print(__doc__)
    print('   k  dim     |Z|  source   poles   short   mu  denominator')
    ok, tot, seen = True, 0, 0
    for k in range(9, 24):
        got = check(k)
        if got is None:
            continue
        seen += 1
        tot += got[0]
        ok = ok and got[1]
    if not seen:
        print('   no data/poles_rec_<k>.npz found')
        return 1
    print()
    print('%s -- %d poles in %d dimension(s)'
          % ('ALL RECORD AXIS LAYERS VERIFIED EXACTLY' if ok else '*** A LAYER FAILED ***',
             tot, seen))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
