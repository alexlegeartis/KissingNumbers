#!/usr/bin/env python3
"""Re-derive run3's OWN certificates by a second, independent arithmetic path.

WHICH NUMBERS.  These are the certificates run3's default path produces: 627 822 180 at
k = 3, which is the programme's optimum and the claim, and 361 275 383 at k = 4, which is 97
BELOW the optimum and is no longer the claim.  Dimension 68 is claimed at 361 275 480 from
common/exact_vertex.py, which solves the dual vertex instead of rounding a float; verify68.py
checks that.  This script is about the ARITHMETIC of run3's path, and both of its targets are
valid lower bounds, so it is left checking exactly what it always checked.

    python scripts/recheck_cert.py          # every LP-certified claim, about a minute
    python scripts/recheck_cert.py 68       # one of them

WHY.  The certificate in ../../../common/kpoint_lp.py earns its speed by rewriting the problem:
it scales every constraint to integers, divides each row by its gcd, and carries the dual as one
integer vector over a single common denominator.  Each of those steps is somewhere to be wrong,
and a wrong one would inflate a claim silently -- the numbers involved are far too large to
sanity-check by eye.

So this script does none of them.  It takes the integer dual that the certificate settled on,
maps it back to the ORIGINAL rows, and evaluates weak duality elementwise in Fraction
arithmetic, one pass:

    slack_j  =  c_j - sum_i y_i A_ij           must be >= 0 wherever x_j is unbounded
    bound    =  y.b - sum over bounded j of u_j * max(0, -slack_j)
    claim    <= ceil(bound)

The map back is the part worth stating: _integralise returns row i scaled by L_i, so the dual
a_i/D of the integral system corresponds to y_i = (a_i/D) * L_i on the original rows.  Getting
that factor wrong is the difference between a certificate and nonsense, which is exactly the
kind of thing this pass exists to catch -- and did, on its first run.

This is a check on the ARITHMETIC.  That the linear programme is the right programme is a
separate matter, settled by scripts/leech_truth.py, which counts the Leech cross-sections
outright and finds the same programme exact at k = 2, 3 and 4."""
import os
import sys
import time
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', '..', 'common'))
import kpoint_lp   # noqa: E402

M1 = {0: 2603658750, 1: 1512243200, 2: 280928256, 3: 13959168, 4: 127800}

# Written out above rather than imported, because this file exists to be a second
# path.  But theta.py derives the same five numbers from Venkov's theorem, and if
# that derivation is ever corrected this script would go on certifying against the
# old distribution without a word.  So: state the agreement.
from theta import onepoint   # noqa: E402
_venkov = onepoint(72, 8, 6218175600, 11)[0]
assert M1 == {c: int(_venkov[c]) for c in _venkov}, (
    'the one-point distribution here no longer matches theta.py: %s vs %s'
    % (M1, {c: int(_venkov[c]) for c in _venkov}))

CLAIMS = {
    68: ([[8, 0, -4, 4], [0, 8, -4, 4], [-4, -4, 8, -4], [4, 4, -4, 8]], 361275383),
    69: ([[8, -4, -4], [-4, 8, 0], [-4, 0, 8]], 627822180),
    # Dimension 70 is the third LP-certified claim and goes through the same
    # rewriting; 71 is an exact count of the cross-section and needs no certificate.
    70: ([[8, 4], [4, 8]], 1249778250),
}
DENOMS = (10 ** 6, 10 ** 9, 10 ** 12, 10 ** 15, 10 ** 18)
CAP = {}
_orig = kpoint_lp._cert_integer


def _replicate(A, b, L, ub, m, idx, tags, Z, sense, y, D):
    """the integer certificate for one D, returning (value, integer dual) or None"""
    nr = len(A)
    SG = int(sense)
    zj = idx[Z]
    free = [k for k in range(m) if ub[k] is None]
    ai = next(i for i, t in enumerate(tags) if t[0] == 'mom' and all(q == 0 for q in t[1]))
    e0 = A[ai][0][1]
    a = [int(round(float(y[i]) / float(L[i]) * D)) for i in range(nr)]
    S = [0] * m
    for i in range(nr):
        if a[i]:
            for k, v in A[i]:
                S[k] += a[i] * v
    T = [(D * SG if k == zj else 0) - S[k] for k in range(m)]
    if min(T[k] for k in free) < 0:
        sh = min(T[k] for k in free) // e0
        a[ai] += sh
        for k in range(m):
            T[k] -= sh * e0
    if min(T[k] for k in free) < 0:
        return None
    num = sum(a[i] * b[i] for i in range(nr))
    for k in range(m):
        if ub[k] is not None and T[k] < 0:
            num -= int(ub[k]) * (-T[k])
    return (-((-num) // D) if sense > 0 else (-num) // D), a


def spy(rows, rhs, ub, m, idx, tags, Z, sense, y):
    """let the real certificate run, and record the dual it won with"""
    A, b, L = kpoint_lp._integralise(rows, rhs)
    best = bestD = None
    for D in DENOMS:
        got = _replicate(A, b, L, ub, m, idx, tags, Z, sense, y, D)
        if got is None:
            continue
        if best is None or (got[0] > best if sense > 0 else got[0] < best):
            best, bestD = got[0], (D, got[1])
    out = _orig(rows, rhs, ub, m, idx, tags, Z, sense, y)
    CAP.setdefault('runs', []).append(
        dict(rows=rows, rhs=rhs, ub=ub, m=m, idx=idx, Z=Z, sense=sense,
             cert=out, D=bestD, L=L, replicated=best))
    return out


def independent(run):
    """weak duality, elementwise in Fractions, against the original rows"""
    rows, rhs, ub, m, idx = run['rows'], run['rhs'], run['ub'], run['m'], run['idx']
    Z, sense, L = run['Z'], run['sense'], run['L']
    D, a = run['D']
    y = [F(int(a[i]), D) * F(L[i]) for i in range(len(a))]
    nz = [i for i in range(len(rows)) if y[i]]
    SG, zj = F(sense), idx[Z]
    worst = None
    pen = F(0)
    for j in range(m):
        acc = F(0)
        for i in nz:
            v = rows[i][j]
            if v:
                acc += y[i] * v
        sl = (SG if j == zj else F(0)) - acc
        if ub[j] is None:
            if worst is None or sl < worst:
                worst = sl
        elif sl < 0:
            pen += F(int(ub[j])) * (-sl)
    bnd = sum((y[i] * rhs[i] for i in nz), F(0)) - pen
    v = (-((-bnd.numerator) // bnd.denominator) if sense > 0
         else (-bnd.numerator) // bnd.denominator)
    return worst, v, len(nz)


if __name__ == '__main__':
    print(__doc__)
    want = [int(x) for x in sys.argv[1:]] or sorted(CLAIMS)
    kpoint_lp._cert_integer = spy
    fail = []
    for dim in want:
        G, claim = CLAIMS[dim]
        CAP.clear()
        t0 = time.time()
        r = kpoint_lp.run3(72, 6218175600, 8, [1, 2, 3, 4], M1, G, tdes=11)
        cert = r[2] if r and r[0] == 'cert' else None
        print()
        print("=" * 86)
        print("dimension %d: k = %d, claim %d" % (dim, len(G), claim))
        print("=" * 86)
        print("  certificate from the engine: %s   (%.0fs)" % (cert, time.time() - t0))
        runs = [z for z in CAP['runs'] if z['cert'] is not None]
        print("  dual candidates certified: %d, values %s"
              % (len(runs), [z['cert'] for z in runs]))
        best = max(runs, key=lambda z: z['cert'])
        assert best['replicated'] == best['cert'], (best['replicated'], best['cert'])
        D, a = best['D']
        worst, v, nz = independent(best)
        print("  winning dual: common denominator %d, %d of %d entries nonzero"
              % (D, nz, len(a)))
        print("  minimum slack over the free cells: %s" % worst)
        print("  independent bound from that dual: %d" % v)
        ok = (worst >= 0) and v == best['cert'] == cert == claim
        print("  %-4s dual feasible, bound reproduced, and equal to the claim"
              % ('PASS' if ok else 'FAIL'))
        if not ok:
            fail.append(dim)
    print()
    print("=" * 86)
    if fail:
        print("*** DISAGREEMENT in dimensions %s -- do not use those values ***"
              % ', '.join(map(str, fail)))
    else:
        print("ALL CHECKS PASSED -- %d certificate%s re-derived independently"
              % (len(want), "" if len(want) == 1 else "s"))
    print("=" * 86)
    raise SystemExit(1 if fail else 0)
