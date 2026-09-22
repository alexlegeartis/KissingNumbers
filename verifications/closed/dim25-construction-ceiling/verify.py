# -*- coding: utf-8 -*-
"""CLOSED: the CEILING of the dimension-25 lens-head CONSTRUCTION.

    python verify.py [SECONDS]        SECONDS = the CP-SAT budget, default 900

What this file proves is a statement about a construction, not about K(25).  The
construction is: the Leech equator in R^24 x {0}, minus the points a head removes; two
caps of heads at height +-1 over norm-6 leans; the two poles; and the non-lattice equator
point.  No removal is ever shared -- a theorem, re-checked below -- so the total is

    196560 + H + 3,      H = (block family) + (supplement),

and H is the only quantity in the whole construction that can still move.  This package
bounds H:

  * the block family is capped at 971 over the ENTIRE norm-6 shell (a full-shell search,
    performed elsewhere; what is verified here is that the shipped 971-head core really
    is a legal configuration with 971 distinct owners);
  * the supplement is a maximum independent set in ONE explicit graph -- the 1188 heads
    that survive against that core, a COMPLETED enumeration over all 77 350 candidate
    owners -- and CP-SAT bounds the independence number of that graph.

Ceiling = 196560 + 971 + (that bound) + 3.  With the bound 89 recorded by the pool
computation this is 197623, and the shipped 197579 sits 44 below it.

Every count printed below is DERIVED from the artefacts.  Nothing asserts a magic number:
this repository has already shipped a verifier that asserted `n == 1006` and thereby
rejected a better configuration, so the assertions here are RELATIONS --
len(owners) == len(heads), equator == 196560 - H, total == equator + 2H + 2 + |E| -- and
the interesting integers are printed, not demanded.

Coordinates.  The Leech minimal vectors are integer vectors of norm 32 (Cohn units); a
head is stored in norm-4 units, where |x|^2 = 3 and an equator point is z = w / sqrt8 with
|z|^2 = 4.  The kissing conditions are <x,z> <= 2 (head against equator) and <x,x'> <= 1
(head against head).

Exactness.  Every core head is a CLASS head x = (u + t- v) / sqrt8 with u minimal, v of
norm 48 and <u,v> = -24, and t- = (3 - sqrt3)/6.  Each of its inner products with an
integer vector is therefore A + t- B with A, B integers, and the test reduces to
"P + Q sqrt3 <= 0" for integers P, Q, decided by comparing squares.  The core is checked
that way, with no floating point deciding anything.  The 1188 pool heads are NOT class
heads (1 of 1188 is); they came out of a numerical search and are stored as floats, so
their checks are floating point, and the worst margin of each is printed so that a reader
can see how close to the bound the data really sits.
"""
from __future__ import print_function

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from leech import load as leech_load                                    # noqa: E402

DATA = os.path.join(HERE, 'data')
SHIPPED = os.path.join(HERE, '..', '..', 'improved', 'dim25-lens-heads', 'data')

SQ8 = np.sqrt(8.0)
TM = (3.0 - np.sqrt(3.0)) / 6.0
TOL = 1e-9

# The screen the pool enumeration ran: a minimal vector z is a CANDIDATE OWNER when the
# core reaches no further than this towards it.  It is an input here, not a derivation --
# it is reproduced only because it identifies WHICH 971-head core the 1188-head pool was
# enumerated over (a second, unrelated 971-head core in the working files screens to a
# different number).
SCREEN = 1.9352350988087748

FAILS = []


def chk(name, ok, extra=''):
    print('  %-62s %s%s' % (name, 'ok' if ok else 'FAIL', ('   ' + extra) if extra else ''))
    if not ok:
        FAILS.append(name)
    return ok


def le_sqrt3_int(P, Q):
    """P + Q sqrt3 <= 0, exactly, elementwise, for int64 arrays with |P|,|Q| < 3e9.

    Same three-case decision as verifications/improved/dim25-lens-heads/verify.py; it is
    twelve lines and lives in both because neither package should have to import the
    other to be checkable on its own.
    """
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
    budget = float(sys.argv[1]) if len(sys.argv) > 1 else 900.0
    t_start = time.time()
    print('verify.py --- the ceiling of the dimension-25 lens-head construction\n')

    # ------------------------------------------------------------------ 1. the equator
    print('1. the Leech minimal vectors, rebuilt from the Golay code by common/leech.py')
    M = leech_load().astype(np.int64)
    NMIN = len(M)
    chk('%d vectors, all of norm 32, all distinct' % NMIN,
        M.shape[1] == 24
        and bool(((M * M).sum(1) == 32).all())
        and len({r.tobytes() for r in M}) == NMIN)
    idx = {M[i].tobytes(): i for i in range(NMIN)}
    Zf = M.astype(np.float64) / SQ8                       # norm-4 units, |z|^2 = 4

    # ------------------------------------------------------------------ 2. the core
    print('\n2. the block-family core')
    X = np.load(os.path.join(DATA, 'rec971_X.npy'))       # norm-4 units, |x|^2 = 3
    own = np.load(os.path.join(DATA, 'rec971_own.npy')).astype(np.int64)
    n = len(X)
    chk('%d core heads, %d owner indices' % (n, len(own)), len(own) == n)
    chk('every owner index names a minimal vector',
        bool(own.min() >= 0) and bool(own.max() < NMIN))
    chk('the %d owners are distinct' % n, len(set(own.tolist())) == n)
    dn = float(np.abs((X * X).sum(1) - 3.0).max())
    chk('every core head has |x|^2 = 3', dn < 1e-12, '%.2e' % dn)

    # the lean: x = (u + t- v)/sqrt8 with v integral of norm 48 and <u,v> = -24
    Vf = (X * SQ8 - M[own].astype(np.float64)) / TM
    V = np.rint(Vf).astype(np.int64)
    drift = float(np.abs(Vf - V).max())
    chk('every core head is a CLASS head (u + t- v)/sqrt8, v integral', drift < 1e-6,
        '%.2e' % drift)
    chk('every lean has |v|^2 = 48 (norm 6) and <u,v> = -24 (that is -3)',
        bool(((V * V).sum(1) == 48).all())
        and bool((((M[own] * V).sum(1)) == -24).all()))
    W = M[own] + V
    chk('u + v is again a minimal vector, at <u,u+v> = 8 (that is 1)',
        bool(((W * W).sum(1) == 32).all())
        and bool((((M[own] * W).sum(1)) == 8).all()))

    leans = {}
    for i, r in enumerate(V):
        leans.setdefault(r.tobytes(), []).append(i)
    L = np.array([np.frombuffer(k, np.int64) for k in leans])
    sizes = sorted((len(v) for v in leans.values()), reverse=True)
    print('   %d distinct leans carry the %d heads, in blocks %s'
          % (len(L), n, ' + '.join(str(s) for s in sizes)))
    GL = (L @ L.T) // 8
    print('   the lean Gram, in norm-4 units (diagonal 6 = norm 6):')
    for row in GL:
        print('      ' + ' '.join('%3d' % v for v in row))
    chk('the blocks account for every head and none twice', sum(sizes) == n)
    chk('no block exceeds 552 = A(23,1/5), the proved size of a whole block',
        max(sizes) <= 552, 'largest block %d' % max(sizes))

    # -------------------------------------------- (a) the core is pairwise compatible
    print('\n3. (a) the core is a legal head set: <x,x> <= 1 on every pair, exactly')
    bad = 0
    tight = 0
    for s in range(0, n, 64):
        a = (M[own][s:s + 64] @ M[own].T) // 8
        mu = ((M[own][s:s + 64] @ V.T) // 8) + ((M[own] @ V[s:s + 64].T) // 8).T
        C = (V[s:s + 64] @ V.T) // 8
        P = 6 * a + 3 * mu + 2 * C - 6
        Q = -(mu + C)
        okm = le_sqrt3_int(P, Q)
        eq = (P == 0) & (Q == 0)
        for k in range(len(okm)):
            okm[k, s + k] = True
            eq[k, s + k] = False
        bad += int((~okm).sum())
        tight += int(eq.sum())
    chk('all %d core pairs meet at <x,x> <= 1' % (n * (n - 1) // 2), bad == 0,
        '%d violations' % bad)
    print('   %d of them meet at exactly 1 (the within-block value)' % (tight // 2))

    # -------------------------------------------- (b) removals are singletons, distinct
    print('\n4. (b) every core head removes exactly its own owner, exactly')
    nrem = np.zeros(n, np.int64)
    wrong_owner = 0
    for s in range(0, n, 32):
        A = (M[own][s:s + 32] @ M.T) // 8
        B = (V[s:s + 32] @ M.T) // 8
        # <x,z> = A + t- B  in norm-4 units;  <= 2  <=>  (6A + 3B - 12) + (-B) sqrt3 <= 0
        okm = le_sqrt3_int(6 * A + 3 * B - 12, -B)
        rem = ~okm
        nrem[s:s + 32] = rem.sum(1)
        for k in range(len(rem)):
            j = int(own[s + k])
            if not (rem[k].sum() == 1 and rem[k][j]):
                wrong_owner += 1
    chk('each head removes exactly one equator point, and it is its owner',
        wrong_owner == 0 and bool((nrem == 1).all()),
        'removal counts seen: %s' % sorted(set(nrem.tolist())))
    removed = set(own.tolist())
    chk('no two heads share a removal, so the deletions total %d' % n,
        len(removed) == n and int(nrem.sum()) == n)
    equator_core = NMIN - n
    chk('equator under the core alone = %d - %d = %d' % (NMIN, n, equator_core),
        equator_core == NMIN - n)

    # WHY no removal is ever shared, measured rather than asserted: one spherical
    # 7-design gives <x,p> >= (18 + sqrt30)/7 = 3.353889... for a head x that removes p,
    # while two COMPATIBLE heads sharing a removal p would both have <x,p> <= 10/3 =
    # 3.333333...  The two intervals do not meet, so a shared removal is impossible; the
    # data must therefore sit above 10/3, and it does.
    ip_own = (X * Zf[own]).sum(1)
    lo = float(ip_own.min())
    chk('every head reaches its owner at more than 10/3 = 3.333333, so no two heads',
        lo > 10.0 / 3.0, 'min %.12f' % lo)
    print('   can ever share a removal (the 7-design floor for a removal is')
    print('   (18 + sqrt30)/7 = %.12f, and the data sits at %.12f)'
          % ((18.0 + np.sqrt(30.0)) / 7.0, lo))

    # ------------------------------------------------------------------ (c) the pool
    print('\n5. (c) the supplement pool: the complete set of heads admissible against it')
    Pl = np.load(os.path.join(DATA, 'p137_poolX.npy'))
    Uo = np.load(os.path.join(DATA, 'p137_poolU.npy')).astype(np.int64)
    m = len(Pl)
    chk('%d pool heads, %d owner indices' % (m, len(Uo)), len(Uo) == m)
    chk('the %d pool owners are distinct' % m, len(set(Uo.tolist())) == m)
    chk('no pool owner is a core owner',
        len(set(Uo.tolist()) & removed) == 0,
        '%d shared' % len(set(Uo.tolist()) & removed))
    dnp = float(np.abs((Pl * Pl).sum(1) - 3.0).max())
    chk('every pool head has |x|^2 = 3', dnp < 1e-9, '%.2e' % dnp)

    worst_other = -9.0
    worst_owner = 9.0
    counts = set()
    for s in range(0, m, 64):
        G = Pl[s:s + 64] @ Zf.T
        counts |= set((G > 2.0 + TOL).sum(1).tolist())
        for k in range(len(G)):
            j = int(Uo[s + k])
            worst_owner = min(worst_owner, float(G[k, j]))
            g = G[k].copy()
            g[j] = -9.0
            worst_other = max(worst_other, float(g.max()))
    chk('each pool head too removes exactly one equator point, its own owner',
        counts == {1}, 'removal counts seen: %s' % sorted(counts))
    chk('a pool head too reaches its owner above 10/3, so it shares no removal',
        worst_owner > 10.0 / 3.0, 'min %.12f' % worst_owner)
    print('   the largest value a pool head reaches anywhere but its own owner is')
    print('   %.12f, over the bound 2 by %.2e' % (worst_other, worst_other - 2.0))

    Vp = np.rint((Pl * SQ8 - M[Uo].astype(np.float64)) / TM).astype(np.int64)
    cls = ((np.abs((Pl * SQ8 - M[Uo].astype(np.float64)) / TM - Vp).max(1) < 1e-6)
           & ((Vp * Vp).sum(1) == 48)
           & ((M[Uo] * Vp).sum(1) == -24))
    print('   %d of the %d pool heads has exact CLASS form; the rest came out of a numerical'
          % (int(cls.sum()), m))
    print('   search and are stored as plain floats, so their checks above are floating point')

    core_max = -9.0
    for s in range(0, m, 128):
        core_max = max(core_max, float((Pl[s:s + 128] @ X.T).max()))
    chk('every pool head is compatible with the WHOLE core, <x,x> <= 1',
        core_max <= 1.0 + TOL, 'max %.15f' % core_max)

    scr = np.full(NMIN, -9.0)
    for s in range(0, n, 64):
        scr = np.maximum(scr, (X[s:s + 64] @ Zf.T).max(0))
    ncand = int((scr <= SCREEN).sum())
    inscr = int((scr[Uo] <= SCREEN).sum())
    print('   the candidate-owner screen of THIS core holds %d minimal vectors,' % ncand)
    print('   and %d of the %d pool owners lie in it' % (inscr, m))
    chk('every pool owner is a candidate owner of this core', inscr == m)

    # ------------------------------------------------ (d) the graph and its CP-SAT bound
    print('\n6. (d) the pool conflict graph, and its independence number from above')
    Gm = Pl @ Pl.T
    np.fill_diagonal(Gm, -9.0)
    Cf = Gm > 1.0 + TOL
    same = Uo[:, None] == Uo[None, :]
    np.fill_diagonal(same, False)
    Cf |= same
    chk('the conflict relation is symmetric and loop-free',
        bool((Cf == Cf.T).all()) and not bool(np.diag(Cf).any()))
    dens = float(Cf.mean())
    print('   %d vertices, %d edges, density %.4f' % (m, int(Cf.sum()) // 2, dens))

    try:
        from ortools.sat.python import cp_model
    except ImportError:
        chk('OR-Tools is installed, so the independence bound can be obtained', False,
            'ortools is MISSING -- no bound was proved here, and none is claimed')
        print('\nVERIFICATION FAILED: %s' % ', '.join(FAILS))
        return 1

    mod = cp_model.CpModel()
    xv = [mod.NewBoolVar('x%d' % i) for i in range(m)]
    cov = np.zeros_like(Cf)
    nc = 0
    for v in np.argsort(-Cf.sum(1)):
        nb = np.nonzero(Cf[v] & ~cov[v])[0]
        if len(nb) == 0:
            continue
        cl = [int(v)]
        for u in nb:
            if all(Cf[u][c] for c in cl):
                cl.append(int(u))
        if len(cl) >= 2:
            mod.AddAtMostOne([xv[c] for c in cl])
            nc += 1
            for a in cl:
                for b in cl:
                    cov[a][b] = True
    rem = np.nonzero(np.triu(Cf & ~cov, 1))
    for a, b in zip(*rem):
        mod.AddAtMostOne([xv[int(a)], xv[int(b)]])
    print('   model: %d clique constraints covering the graph, %d residual pairs'
          % (nc, len(rem[0])))
    mod.Maximize(sum(xv))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = budget
    solver.parameters.num_search_workers = 8
    t0 = time.time()
    st = solver.Solve(mod)
    secs = time.time() - t0
    name = solver.StatusName(st)
    print('   CP-SAT %s in %.0f s of a %.0f s budget' % (name, secs, budget))
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        chk('CP-SAT returned a usable bound on the independence number', False,
            'status %s -- NO bound was obtained, and none is claimed' % name)
        print('\nVERIFICATION FAILED: %s' % ', '.join(FAILS))
        return 1
    realised = int(round(solver.ObjectiveValue()))
    raw = float(solver.BestObjectiveBound())
    B = int(np.floor(raw + 1e-6))
    print('   independent set realised %d, upper bound %.4f -> alpha <= %d'
          % (realised, raw, B))
    chk('the bound is an upper bound on the set the same solve realised', B >= realised)
    if name == 'OPTIMAL':
        print('   the solve CLOSED: alpha of this graph is exactly %d' % B)
    else:
        print('   the solve did NOT close; %d is the best bound it proved, which is what'
              % B)
        print('   the ceiling below is built from -- it is an upper bound, not alpha')
    print('   (the bound this package records is 89; this run proved %d)' % B)

    # --------------------------------------------------------------- (e) the arithmetic
    print('\n7. (e) the ceiling, and the shipped configuration against it')
    Xs = np.load(os.path.join(SHIPPED, 'heads_X.npy'))
    Us = np.load(os.path.join(SHIPPED, 'heads_U.npy')).astype(np.int64)
    Es = np.load(os.path.join(SHIPPED, 'extra_P.npy'))
    Hs = len(Xs)
    chk('the shipped package carries %d heads and %d owners' % (Hs, len(Us)),
        len(Us) == Hs)
    si = set(idx[r.tobytes()] for r in Us)
    chk('the shipped owners are distinct and all minimal', len(si) == Hs)
    chk('the shipped configuration contains this whole core', removed <= si)
    supp = len(si & set(Uo.tolist()))
    chk('the rest of it, %d heads, comes from this very pool' % supp,
        len(removed) + supp == Hs, '%d + %d = %d' % (len(removed), supp, Hs))
    chk('so the shipped supplement is an independent set in the graph above, of size %d'
        % supp, supp <= B, '%d <= %d' % (supp, B))
    sel = np.array([i for i in range(m) if int(Uo[i]) in si], dtype=np.int64)
    nsc = int(cls[sel].sum()) if len(sel) else 0
    print('   two accountings of the same %d heads: %d core + %d pool here, and'
          % (Hs, n, supp))
    print('   %d class + %d other in the shipped package, because %d of the %d supplementary'
          % (n + nsc, Hs - n - nsc, nsc, supp))
    print('   heads turns out to have exact class form over a further lean')

    nE = len(Es)
    eq_ship = NMIN - Hs
    tot_ship = eq_ship + 2 * Hs + 2 + nE
    chk('shipped: equator %d + 2 x %d heads + 2 poles + %d extra = %d'
        % (eq_ship, Hs, nE, tot_ship),
        eq_ship == NMIN - Hs and tot_ship == eq_ship + 2 * Hs + 2 + nE, str(tot_ship))

    Hmax = n + B
    eq_max = NMIN - Hmax
    tot_max = eq_max + 2 * Hmax + 2 + nE
    chk('ceiling: H <= %d + %d = %d, so the total is %d + %d + %d = %d'
        % (n, B, Hmax, NMIN, Hmax, 2 + nE, tot_max),
        Hmax == n + B and tot_max == NMIN + Hmax + 2 + nE, str(tot_max))
    chk('the shipped configuration does not exceed the ceiling', tot_ship <= tot_max,
        '%d <= %d' % (tot_ship, tot_max))
    print('   the shipped %d is %d below the ceiling %d' % (tot_ship, tot_max - tot_ship, tot_max))
    print('   what is left in this construction: %d heads, worth %d points'
          % (Hmax - Hs, tot_max - tot_ship))

    print('\n   ran in %.0f s' % (time.time() - t_start))
    if FAILS:
        print('\nVERIFICATION FAILED: %s' % ', '.join(FAILS))
        return 1
    print('\nALL CHECKS PASS')
    print('   the construction stops at %d + %d + %d + %d = %d'
          % (NMIN, n, B, 2 + nE, tot_max))
    print('   the shipped configuration is %d, within %d of it'
          % (tot_ship, tot_max - tot_ship))
    return 0


if __name__ == '__main__':
    sys.exit(main())
