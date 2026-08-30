#!/usr/bin/env python3
"""The EXACT optimum of the k-point programme, by solving for the dual vertex itself.

    from exact_vertex import exact_optimum
    exact_optimum(72, 6218175600, 8, [1,2,3,4], M1, G)     # -> (bound, report)

WHY THIS EXISTS.  run3 proposes a dual with a floating-point solve, rounds it, and repairs it
into exact feasibility by lowering one coordinate.  That is sound -- an imprecise dual weakens
the bound and cannot inflate it -- but it LOSES value, and how much it loses depends on the
conditioning, which depends on how the Gram matrix happens to be presented.  At k = 3 the
sixteen bases of one scaled A_3 came out as thirteen different certificates that way; at k = 4
the loss against the true optimum was 97.

TWO STEPS FIX IT.

(1) REDUCE.  Most rows are identifications, not equations: n[c] = n[-c] from the antipodal
    symmetry, n[c] = n[c'] or n[c] = 0 or n[c] = 1 from the translations.  Substituting them
    out is exact -- it leaves the feasible set alone -- and it takes k = 4 from 3473 rows and
    1681 columns down to 617 and 520, dense and well scaled.  HiGHS then returns the same
    optimum from every method and every presolve setting, where before it returned points
    violating the constraints by several whole vectors.

(2) TAKE THE VERTEX, NOT A ROUNDING OF IT.  At the optimum, complementary slackness makes the
    dual tight exactly on the active columns, so y solves a linear system with rational data.
    Read the active set off the numeric dual, solve that system in Fractions, and then VERIFY
    A^T y <= c over every column exactly.  The verification is what makes it a certificate:
    a wrong active set fails the check and nothing is returned.

The result is floor(b.y), the largest integer this programme can certify.  Weak duality is
unchanged (Lemma 9 of the paper); only the witness is better.

CHECKED against the two claims run3 already certifies exactly, dimensions 70 and 69, which
this reproduces to the unit, and against the Leech cross-sections that leech_truth.py counts
outright.
"""
from fractions import Fraction as F

import numpy as np
from scipy.optimize import linprog

from kpoint_lp import build_system


def reduce_system(rows, rhs, ub, m):
    """Substitute out every row that is an identification.  Returns
    (rows, rhs, ub, groups, fixed_value_of_zero_class, representative_of, ok)."""
    parent = list(range(m))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    fixed, used, ok = {}, [False] * len(rows), True
    for r, (row, rv) in enumerate(zip(rows, rhs)):
        nz = [(j, x) for j, x in enumerate(row) if x]
        if len(nz) == 2 and rv == 0 and nz[0][1] == -nz[1][1]:
            ra, rb = find(nz[0][0]), find(nz[1][0])
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
            used[r] = True
        elif len(nz) == 1:
            fixed[nz[0][0]] = F(rv) / nz[0][1]
            used[r] = True

    cls = {}
    for j, v in fixed.items():
        rj = find(j)
        if rj in cls and cls[rj] != v:
            ok = False
        cls[rj] = v
    reps = sorted({find(j) for j in range(m)} - set(cls))
    pos = {r: t for t, r in enumerate(reps)}
    groups = [[] for _ in reps]
    for j in range(m):
        if find(j) in pos:
            groups[pos[find(j)]].append(j)

    nrows, nrhs = [], []
    for r, (row, rv) in enumerate(zip(rows, rhs)):
        if used[r]:
            continue
        acc = [F(0)] * len(reps)
        const = F(0)
        for j, x in enumerate(row):
            if not x:
                continue
            rj = find(j)
            if rj in cls:
                const += x * cls[rj]
            else:
                acc[pos[rj]] += x
        b = F(rv) - const
        if not any(acc):
            if b != 0:
                ok = False
            continue
        nrows.append(acc)
        nrhs.append(b)
    nub = []
    for g in groups:
        caps = [ub[j] for j in g if ub[j] is not None]
        nub.append(min(caps) if caps else None)
    return nrows, nrhs, nub, groups, cls, find, pos, ok


def _exact_dual(rows, rhs, ub, m, zt, N):
    """(bound, note) from the exact dual vertex, or (None, why)."""
    nr = len(rows)
    sc = [max(abs(x) for x in row if x) for row in rows]
    Af = np.array([[float(x / s) for x in row] for row, s in zip(rows, sc)])
    bf = np.array([float(F(rv) / s) / N for rv, s in zip(rhs, sc)])
    c = np.zeros(m)
    c[zt] = 1.0
    bnds = [(0.0, None if u is None else float(u) / N) for u in ub]
    best = None
    for meth in ('highs-ipm', 'highs', 'highs-ds'):
        r = linprog(c=c, A_eq=Af, b_eq=bf, bounds=bnds, method=meth,
                    options={'presolve': False})
        if r.success and (best is None or r.fun > best.fun):
            best = r
    if best is None:
        return None, 'no numeric solve succeeded'
    marg = np.asarray(best.eqlin.marginals) / np.array(sc)
    sl = c - marg @ np.array([[float(x) for x in row] for row in rows])
    scale = max(1.0, float(np.abs(sl).max()))
    B = [j for j in range(m) if abs(sl[j]) <= 1e-9 * scale]
    if not B:
        return None, 'no active columns'

    # solve  sum_i y_i A_ij = c_j  for j in B, exactly
    aug = [[F(rows[i][j]) for i in range(nr)] + [F(1) if j == zt else F(0)] for j in B]
    piv, r0 = [], 0
    for col in range(nr):
        pr = next((t for t in range(r0, len(aug)) if aug[t][col]), None)
        if pr is None:
            continue
        aug[r0], aug[pr] = aug[pr], aug[r0]
        pv = aug[r0][col]
        aug[r0] = [x / pv for x in aug[r0]]
        for t in range(len(aug)):
            if t != r0 and aug[t][col]:
                f = aug[t][col]
                aug[t] = [x - f * z for x, z in zip(aug[t], aug[r0])]
        piv.append(col)
        r0 += 1
        if r0 == len(aug):
            break
    for t in range(r0, len(aug)):
        if aug[t][nr] != 0:
            return None, 'the active system is inconsistent'
    y = [F(0)] * nr
    for t, col in enumerate(piv):
        y[col] = aug[t][nr]

    slack = []
    for j in range(m):
        s = F(1) if j == zt else F(0)
        for i in range(nr):
            if y[i] and rows[i][j]:
                s -= y[i] * rows[i][j]
        slack.append(s)
    free = [j for j in range(m) if ub[j] is None]
    bad = sum(1 for j in free if slack[j] < 0)
    if bad:
        return None, 'dual infeasible on %d columns' % bad
    pen = sum(F(ub[j]) * max(F(0), -slack[j]) for j in range(m) if ub[j] is not None)
    val = sum(y[i] * F(rhs[i]) for i in range(nr) if y[i]) - pen
    return -((-val.numerator) // val.denominator), ('exact vertex, %d active columns, '
                                                    '%d dual coordinates' % (len(B), len(piv)))


def exact_optimum(n, N, mu, IPS, M1, G, AMAX=2, tdes=11):
    """(bound, report). The bound is floor of the programme's exact optimum, or None."""
    S = build_system(n, N, mu, IPS, M1, G, AMAX, tdes)
    m, idx = S['m'], S['idx']
    zero = tuple([0] * S['K'])
    if zero not in idx:
        return None, 'the all-zero cell is not in the system'
    rows, rhs, ub, groups, cls, find, pos, ok = reduce_system(S['rows'], S['rhs'], S['ub'], m)
    if not ok:
        return None, 'the identifications are contradictory'
    zj = find(idx[zero])
    if zj in cls:
        return int(cls[zj]), 'the zero cell is fixed by substitution'
    v, note = _exact_dual(rows, rhs, ub, len(groups), pos[zj], N)
    return v, '%d cells -> %d variables, %d rows; %s' % (m, len(groups), len(rows), note)

if __name__ == '__main__':
    import sys
    import time
    print(__doc__)
    print("=" * 78)
    print("VALIDATION: the exact vertex against answers that are already known")
    print("=" * 78)
    M1L = {0: 93150, 1: 47104, 2: 4600}
    M1G = {0: 2603658750, 1: 1512243200, 2: 280928256, 3: 13959168, 4: 127800}
    CASES = [
        # the Leech cross-sections leech_truth.py COUNTS outright -- not certificates
        ((24, 196560, 4, [1, 2], M1L, [[4, 0], [0, 4]]), 43164,
         'Leech k=2 orthogonal  -> dim 22, counted'),
        ((24, 196560, 4, [1, 2], M1L, [[4, 1], [1, 4]]), 44550,
         'Leech k=2 cos 1/4     -> dim 22, counted'),
        ((24, 196560, 4, [1, 2], M1L, [[4, 2], [2, 4]]), 49896,
         'Leech k=2 60 degrees  -> dim 22, counted'),
        # P_48, a second lattice, where minimum = maximum so the programme PINS the
        # cross-section rather than bounding it
        ((48, 52416000, 6, [1, 2, 3], {0: 23766960, 1: 12608784, 2: 1678887, 3: 36848},
          [[6, 3], [3, 6]]), 12309600, 'P_48  k=2 60 degrees  -> dim 46, pinned'),
        # and the two Gamma_72 claims run3 already certifies exactly
        ((72, 6218175600, 8, [1, 2, 3, 4], M1G, [[8, -4], [-4, 8]]), 1249778250,
         'Gamma_72 k=2          -> dim 70, claimed'),
        ((72, 6218175600, 8, [1, 2, 3, 4], M1G, [[8, -4, -4], [-4, 8, 0], [-4, 0, 8]]),
         627822180, 'Gamma_72 k=3          -> dim 69, claimed'),
    ]
    ok = True
    for args, want, label in CASES:
        t0 = time.time()
        v, rep = exact_optimum(*args)
        good = (v == want)
        ok = ok and good
        print("   %-38s %-12s expected %-12s %s   (%.0fs)"
              % (label, v, want, "MATCH" if good else "*** MISMATCH ***", time.time() - t0))
        print("      %s" % rep)
    print()
    print("A counted value is the strongest check there is: 49896 is the kissing number of")
    print("Lambda_22, obtained by enumeration, and the programme's optimum equals it.")
    print()
    print("ALL CHECKS PASSED" if ok else "*** FAILED ***")
    sys.exit(0 if ok else 1)
