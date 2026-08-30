"""Bachoc-Vallentin three-point SDP bound for the class problem:
   antipodal spherical codes in S^23 whose distinct-point inner products lie in
   {-1, -1/4, 0, 1/4} (i.e. sets of LINES in R^24 with pairwise |cos| in {0,1/4}).
   Feasibility SDP in the triple distribution, bisected on M = 2 * (#lines)."""
import numpy as np, itertools, sys
from fractions import Fraction as Fr
import cvxpy as cp

n = 24
SIG = [-1.0, -0.25, 0.0, 0.25, 1.0]

def geg(d, K):
    """list of numpy poly evaluators: returns function t->P_k^d(t) values, k=0..K"""
    def ev(t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        P = [np.ones_like(t), t.copy()]
        for k in range(1, K):
            P.append(((2*k+d-2)*t*P[k] - k*P[k-1])/(k+d-2))
        return np.array(P[:K+1])
    return ev

def Yk(k, m, u, v, t):
    g2 = (1-u*u)*(1-v*v)
    if g2 <= 1e-14:
        if k > 0: return np.zeros((m+1, m+1))
        pk = 1.0; pref = 1.0
    else:
        g = np.sqrt(g2); arg = (t-u*v)/g
        arg = min(max(arg, -1.0), 1.0)
        pk = float(geg(n-1, k)(arg)[k][0]); pref = g**k
    I = np.arange(m+1)
    return np.outer(u**I, v**I)*pref*pk

def realizable(u, v, t):
    return 1 + 2*u*v*t - u*u - v*v - t*t >= -1e-12

def build(M, dmax=6, K=30, extra=True, verbose=False):
    trips = [tr for tr in itertools.product(SIG, repeat=3) if realizable(*tr)]
    # S_3 orbits on the three positions
    orb = {}
    for tr in trips:
        key = tuple(sorted(tr))
        orb.setdefault(key, []).append(tr)
    keys = sorted(orb)
    F = {kk: cp.Variable(nonneg=True) for kk in keys}          # value on the whole orbit
    Fv = {tr: F[tuple(sorted(tr))] for tr in trips}
    a  = {s: cp.Variable(nonneg=True) for s in SIG}
    con = [a[1.0] == 1, a[-1.0] == 1, a[0.25] == a[-0.25],
           cp.sum([a[s] for s in SIG]) == M]
    # degeneracy identities
    for tr in trips:
        u,v,t = tr
        if abs(u) == 1:
            con.append(Fv[tr] == (a[v] if abs(t-u*v) < 1e-12 else 0))
        elif abs(v) == 1:
            con.append(Fv[tr] == (a[u] if abs(t-u*v) < 1e-12 else 0))
        elif abs(t) == 1:
            con.append(Fv[tr] == (a[u] if abs(v-u*t) < 1e-12 else 0))
    # marginals
    for u in SIG:
        con.append(cp.sum([Fv[tr] for tr in trips if tr[0] == u]) == M*a[u])
    # two-point LP
    ev = geg(n, K)
    for k in range(1, K+1):
        con.append(cp.sum([a[s]*float(ev(s)[k][0]) for s in SIG]) >= 0)
    # three-point PSD
    for k in range(1, dmax+1):
        Mk = sum(Fv[tr]*Yk(k, dmax, *tr) for tr in trips)
        con.append(Mk >> 0)
    if extra:
        for u in SIG:
            con.append(cp.sum([Fv[tr] for tr in trips if tr[0] == u and tr[1] == u]) >= a[u]**2)
    return cp.Problem(cp.Minimize(0), con)

if __name__ == "__main__":
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    lo, hi = 100.0, 2000.0                    # M = 2*lines
    def feas(M):
        p = build(M, dmax=dmax)
        try:
            p.solve(solver=cp.SCS, eps=1e-8, max_iters=200000, verbose=False)
        except Exception as e:
            print("solver error", e); return True
        return p.status not in ("infeasible", "infeasible_inaccurate")
    print("dmax =", dmax)
    for M in [400, 500, 600, 700, 800, 850, 851, 900, 1000]:
        print("M = %5d  feasible: %s" % (M, feas(float(M))), flush=True)
