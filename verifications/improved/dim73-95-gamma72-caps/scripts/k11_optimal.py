#!/usr/bin/env python3
"""Dimension 83's partition is OPTIMAL, not merely the best found.  CP-SAT, in exact integers.

    python k11_optimal.py            needs OR-Tools; not part of run_all.py

The k = 11 cap set is the 604-point record configuration of R^11, realised over Z[sqrt2] as
(A + B sqrt2)/3 with an integer metric c, so 9<x,y> = P + Q sqrt2 with P, Q integers and the
scaled norm N = 36.  The gain is proportional to 2T + P over a partition into zero-sum triples
and pairs at cosine <= -1/2 (KNOWLEDGE.md section 50), and the shipped certificate gives
84 triples + 176 pairs = 344 units.

KNOWLEDGE.md section 106 recorded 3000 restarts of the packing all returning 84 triples and
read that invariance as a ceiling.  It is one -- but an invariance across a sweep is a
measurement of the sweep, and the same greedy was nine, thirteen, twenty-five and sixty
triples short at k = 17, 19, 21 and 23 (section 107).  So the question is asked exactly here
instead.  The instance is small enough to settle outright: 604 points, 1896 zero-sum triples,
20 006 admissible pairs.

  * maximise 2T + P over triples and pairs, each point used at most once   ->  344, OPTIMAL
  * maximise T alone                                                       ->   84, OPTIMAL

so 84 + 176 is the exact optimum and dimension 83's value cannot be improved by repacking this
configuration.  Both solves take seconds.

Everything below is integer arithmetic.  A pair part needs cos <= -1/2, i.e. <x,y> <= -N/2,
i.e. P + N/2 <= -Q sqrt2, decided by sign analysis and P^2 against 2Q^2 -- the same comparator
scripts/verify11.py uses, and the one whose sign was reversed there until 2026-08-30.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')


def le(L, R):
    """L <= R sqrt2, for integers L and R"""
    if R >= 0 and L <= 0:
        return True
    if R < 0 and L > 0:
        return False
    return (L * L <= 2 * R * R) if L > 0 else (L * L >= 2 * R * R)


def main():
    try:
        from ortools.sat.python import cp_model
    except ImportError:
        print('OR-Tools is not installed; this script is not part of run_all.py.')
        return 0
    print(__doc__)
    Z = np.load(os.path.join(DATA, 'rec11_exact.npz'))
    A, B, c = (Z['A'].astype(np.int64), Z['B'].astype(np.int64), Z['c'].astype(np.int64))
    n = len(A)
    P = (A * c) @ A.T + 2 * ((B * c) @ B.T)
    Q = (A * c) @ B.T + (B * c) @ A.T
    N = int(P[0, 0])
    half = -N // 2
    assert (np.diag(P) == N).all() and (np.diag(Q) == 0).all(), 'norms are not equal'
    key = {(A[i].tobytes(), B[i].tobytes()): i for i in range(n)}

    TRI = set()
    for i in range(n):
        for q in np.nonzero((P[i] == half) & (Q[i] == 0))[0]:
            t = key.get(((-(A[i] + A[q])).tobytes(), (-(B[i] + B[q])).tobytes()))
            if t is not None:
                TRI.add(tuple(sorted((i, int(q), t))))
    TRI = sorted(TRI)
    PAIR = [(i, j) for i in range(n) for j in range(i + 1, n)
            if le(int(P[i, j]) - half, -int(Q[i, j]))]
    print('   %d points, %d zero-sum triples, %d admissible pairs' % (n, len(TRI), len(PAIR)))
    # the predicate must be able to REJECT: witnesses at cos = +1/2 and cos = 0 from this Gram
    for want in (-half, 0):
        w = np.argwhere((P == want) & (Q == 0) & ~np.eye(n, dtype=bool))
        assert len(w), 'no cos = %d/%d witness in this Gram' % (want, N)
        i, j = int(w[0][0]), int(w[0][1])
        assert not le(int(P[i, j]) - half, -int(Q[i, j])), 'the pair test accepts cos >= 0'
    print('   the pair test rejects a cos = +1/2 and a cos = 0 witness: True')

    shipped = np.load(os.path.join(DATA, 'rec11_parts.npz'))
    T0, P0 = len(shipped['triples']), len(shipped['pairs'])
    print('   shipped: %d triples + %d pairs = %d units' % (T0, P0, 2 * T0 + P0))
    print()
    for name, obj in (('2T + P', 'units'), ('T', 'triples')):
        m = cp_model.CpModel()
        xt = [m.NewBoolVar('t%d' % i) for i in range(len(TRI))]
        xp = [m.NewBoolVar('p%d' % i) for i in range(len(PAIR))]
        inc = [[] for _ in range(n)]
        for i, t in enumerate(TRI):
            for v in t:
                inc[v].append(xt[i])
        for i, p in enumerate(PAIR):
            for v in p:
                inc[v].append(xp[i])
        for e in inc:
            m.AddAtMostOne(e)
        m.Maximize(2 * sum(xt) + sum(xp) if obj == 'units' else sum(xt))
        s = cp_model.CpSolver()
        s.parameters.max_time_in_seconds = 900.0
        s.parameters.num_search_workers = 8
        t0 = time.time()
        r = s.Solve(m)
        ok = r == cp_model.OPTIMAL
        val, bound = int(s.ObjectiveValue()), int(s.BestObjectiveBound())
        print('   maximise %-7s : %-8s %d (bound %d) in %.0f s'
              % (name, s.StatusName(r), val, bound, time.time() - t0))
        assert ok, 'not solved to optimality'
        assert val == (2 * T0 + P0 if obj == 'units' else T0), (
            'the shipped partition is not optimal: %d against %d' % (T0 if obj else 0, val))
    print()
    print('   84 triples is the exact maximum and 84 + 176 the exact optimum on 2T + P,')
    print('   so dimension 83 cannot be improved by repacking this configuration.')
    print()
    print('ALL CHECKS PASSED')
    return 0


if __name__ == '__main__':
    sys.exit(main())
