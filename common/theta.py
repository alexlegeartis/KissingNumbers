#!/usr/bin/env python3
"""Extremal even unimodular lattices: theta series and the one-point distribution.

    from theta import extremal_theta, onepoint

Two exact computations, both in rational arithmetic and both with a built-in surplus
check, used by the dimension 46/47 and 70/71 packages.

extremal_theta(n)
    The theta series of an *extremal* even unimodular lattice of dimension n is the unique
    modular form of weight n/2 for SL_2(Z) with constant term 1 whose coefficients of
    q^1 .. q^(n/24) all vanish.  It is computed here from E_4 and Delta by exact linear
    algebra over Q.  For n = 48 it gives minimum 6 and 52416000 minimal vectors; for
    n = 72, minimum 8 and 6218175600.  The lattice does not have to be named: those two
    numbers are properties of the *class*, which is why the dimension 46/47 statements hold
    for all four of P_48p, P_48q, P_48m, P_48n and the dimension 70/71 ones for Gamma_72.

onepoint(n, mu, N, t)
    VENKOV's theorem: for an extremal even unimodular lattice of dimension n = 0 (mod 24)
    the minimal vectors, rescaled to the unit sphere, form a spherical 11-design (7-design
    for E_8's dimension 8).  Hence for a fixed minimal vector u and every p <= t/2

        sum_v <u/|u|, v/|v|>^(2p)  =  N * (2p-1)!! / (n(n+2)...(n+2p-2))     exactly.

    The possible values of c = <u,v> are the integers with |c| <= mu/2 -- because
    |u-v|^2 = 2mu - 2c must be a lattice norm, hence >= mu -- together with c = +-mu for
    v = +-u.  That is mu/2 + 1 unknowns against t/2 + 1 equations, so the system is
    OVER-DETERMINED; the surplus equations are a genuine check, and they hold.  The
    function returns (solution, surplus_ok, number_of_surplus_equations).

    n[0] counts the minimal vectors orthogonal to u.  They span an (n-1)-dimensional
    space, all have norm mu, and any two are at angle >= 60 degrees, so tau(n-1) >= n[0].

    | lattice   | n  | mu | n[0]       | is                    |
    |-----------|----|----|------------|-----------------------|
    | E_8       |  8 |  2 | 126        | E_7                   |
    | Leech     | 24 |  4 | 93150      | Lambda_23             |
    | P_48      | 48 |  6 | 23766960   | dimension 47          |
    | Gamma_72  | 72 |  8 | 2603658750 | dimension 71          |

    python theta.py    checks all four against their known values.
"""
from fractions import Fraction as F

P = 40      # number of q-coefficients carried


def _mul(a, b):
    r = [0] * P
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if i + j < P and y:
                    r[i + j] += x * y
    return r


def _sig3(m):
    return sum(d ** 3 for d in range(1, m + 1) if m % d == 0)


E4 = [1] + [240 * _sig3(m) for m in range(1, P)]
DELTA = [0] * P
DELTA[1] = 1
for _nn in range(1, P):
    _f = [0] * P
    _f[0] = 1
    if _nn < P:
        _f[_nn] = -1
    for _ in range(24):
        DELTA = _mul(DELTA, _f)


def _power(a, e):
    r = [0] * P
    r[0] = 1
    for _ in range(e):
        r = _mul(r, a)
    return r


def extremal_theta(n):
    """q-expansion of the extremal even unimodular theta series of dimension n."""
    w = n // 2
    basis = []
    b = 0
    while 12 * b <= w:
        if (w - 12 * b) % 4 == 0:
            basis.append(((w - 12 * b) // 4, b))
        b += 1
    S = [_power(E4, a) if bb == 0 else _mul(_power(E4, a), _power(DELTA, bb))
         for a, bb in basis]
    k = n // 24
    m = len(S)
    assert m == k + 1, (n, m, k)
    A = [[F(S[j][i]) for j in range(m)] for i in range(k + 1)]
    rhs = [F(1)] + [F(0)] * k
    for i in range(m):                                   # exact Gauss-Jordan over Q
        p = next(r for r in range(i, m) if A[r][i] != 0)
        A[i], A[p] = A[p], A[i]
        rhs[i], rhs[p] = rhs[p], rhs[i]
        pv = A[i][i]
        A[i] = [x / pv for x in A[i]]
        rhs[i] = rhs[i] / pv
        for r in range(m):
            if r != i and A[r][i] != 0:
                f2 = A[r][i]
                A[r] = [x - f2 * y for x, y in zip(A[r], A[i])]
                rhs[r] -= f2 * rhs[i]
    return [sum(rhs[j] * S[j][i] for j in range(m)) for i in range(P)]


def minimum_and_kissing(n):
    """(mu, N) read off the extremal theta series: the first non-zero coefficient."""
    th = extremal_theta(n)
    for i in range(1, P):
        if th[i] != 0:
            assert th[i].denominator == 1, th[i]
            return 2 * i, int(th[i])
    raise ValueError(n)


def onepoint(n, mu, N, t):
    """Venkov one-point distribution.  t is the design order (11, or 7 for E_8).

    Returns (sol, surplus_ok, n_surplus) where sol[c] = #{v minimal : <u,v> = c} for
    c = 0..mu/2 and each non-zero c is attained by that many v of each sign."""
    PM = t // 2
    hm = mu // 2
    unk = list(range(0, hm + 1))
    rows, rhs = [], []
    for p in range(0, PM + 1):
        num = F(1)
        for j in range(1, 2 * p, 2):
            num *= j
        den = F(1)
        for j in range(p):
            den *= (n + 2 * j)
        # +-u contribute 2 * (mu/mu)^{2p} = 2 to every moment; move them to the right side
        rows.append([F(2 if c > 0 else 1) * F(c, mu) ** (2 * p) for c in unk])
        rhs.append(F(N) * F(num) / den - 2)
    K = len(unk)
    A = [rows[i][:] for i in range(K)]
    B = [rhs[i] for i in range(K)]
    for i in range(K):
        p = next(r for r in range(i, K) if A[r][i] != 0)
        A[i], A[p] = A[p], A[i]
        B[i], B[p] = B[p], B[i]
        pv = A[i][i]
        A[i] = [x / pv for x in A[i]]
        B[i] = B[i] / pv
        for r in range(K):
            if r != i and A[r][i] != 0:
                f = A[r][i]
                A[r] = [x - f * y for x, y in zip(A[r], A[i])]
                B[r] -= f * B[i]
    sol = dict(zip(unk, B))
    surplus = list(range(K, len(rows)))
    ok = all(sum(rows[i][j] * sol[unk[j]] for j in range(K)) == rhs[i] for i in surplus)
    return sol, ok, len(surplus)


if __name__ == '__main__':
    print(__doc__)
    ok = True
    print("=" * 78)
    print("extremal theta series")
    print("=" * 78)
    for n, mu_e, N_e in ((8, 2, 240), (24, 4, 196560), (48, 6, 52416000),
                         (72, 8, 6218175600)):
        mu, N = minimum_and_kissing(n)
        good = (mu, N) == (mu_e, N_e)
        ok = ok and good
        print("   n = %-3d  minimum %d, kissing number %-12d  expected (%d, %d)  %s"
              % (n, mu, N, mu_e, N_e, mu_e, N_e, ) if False else
              "   n = %-3d  minimum %d, kissing number %-12d  expected (%d, %d)  %s"
              % (n, mu, N, mu_e, N_e, "MATCH" if good else "*** MISMATCH ***"))

    print()
    print("=" * 78)
    print("one-point distributions, and the cross-sections they force")
    print("=" * 78)
    for n, mu, N, t, truth, name in ((8, 2, 240, 7, 126, "E_8      -> E_7"),
                                     (24, 4, 196560, 11, 93150, "Leech    -> Lambda_23"),
                                     (48, 6, 52416000, 11, 23766960, "P_48     -> dim 47"),
                                     (72, 8, 6218175600, 11, 2603658750, "Gamma_72 -> dim 71")):
        sol, sok, ns = onepoint(n, mu, N, t)
        integral = all(v.denominator == 1 and v >= 0 for v in sol.values())
        total = sum((2 if c else 1) * sol[c] for c in sol) + 2
        good = sol[0] == truth and sok and integral and total == N
        ok = ok and good
        print("   %-22s n[0] = %-12s known %-12s %s"
              % (name, sol[0], truth, "MATCH" if sol[0] == truth else "*** MISMATCH ***"))
        print("      %s  |  integral and non-negative: %s  |  total = N: %s  |  "
              "%d surplus equation(s) hold: %s"
              % (" ".join("n[%d]=%s" % (c, sol[c]) for c in sorted(sol)),
                 integral, total == N, ns, sok))

    print()
    print("ALL CHECKS PASSED" if ok else "*** FAILURES ***")
    raise SystemExit(0 if ok else 1)
