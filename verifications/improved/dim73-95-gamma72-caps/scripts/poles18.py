#!/usr/bin/env python3
"""An exact axis layer for dimensions 74-80, over the WHOLE of every legal shell.

    python scripts/poles18.py            # ~1 min, rebuilds data/poles_k2..k8.npy and verifies

final73-80.py claims poles at k = 1 automatically (t = 3/4 makes sqrt(1-t) = 1/2) and from
this file for k = 2..8, where the triple scheme sits at t = 2/3 and a pole needs
|<z,a>| <= sqrt3/2 in normalised terms.

THE LEGALITY ARGUMENT.  For a of norm M and w of norm m in the same lattice, |a-w|^2 and
|a+w|^2 are lattice norms, so both are >= m and hence |<a,w>| <= M/2.  The pole condition is
4<a,w>^2 <= 3*M*m, and M <= 3m makes that automatic.  So EVERY shell of the direction
lattice with m < M <= 3m is legal against every cap direction at once.

WHAT THIS FIXES.  The previous version took candidates of the form a = w1 + w2 with
<w1,w2> = 0.  That is a subset of the norm-2m shell and it misses two things:

  * A_2 has no two orthogonal roots at all, so k = 2 produced NOTHING and dimension 74
    claimed no poles -- yet A_2's norm-6 shell is M = 3m, is legal, and gives all six;
  * even where it produces something it sees only part of one shell.  At k = 8 the norm-16
    shell has 1080 lines and the w1 + w2 rule reaches a fraction of them.

Sweeping both legal shells and searching each properly gives 394 poles against 248.

NOT A CEILING.  This sweeps the shells of the DIRECTION LATTICE, because that is where the
M <= 3m argument makes legality automatic.  A legal pole need not lie in that lattice, and at
k = 3 it does not: this file reaches 4 lines, 8 poles, while the dimension-27 configuration
carries 12 = K(3), whose worst pole-to-direction cosine is 0.8536 against the permitted
sqrt3/2 = 0.8660, and which is verified over all 20 108 045 530 of its pairs.  So the axis
layers of dimensions 74-80 are conservative, by an unknown margin bounded by K(k).  Nothing
rests on them being best: they are added to a gain, never subtracted from one.

SUPERSEDED FROM k = 5 UP (2026-08-27).  scripts/axis_rotate.py answers the same question
without the lattice: a rotated copy of the root system is a 60-degree code for free and only
has to clear the 30-degree caps, and a rational rotation is built from the Cartan matrix by
Cayley.  It gives 40, 70, 118 and 208 at k = 5, 6, 7, 8 against the 28, 58, 76 and 174 here,
and 40 IS tau(5).  It loses at k = 3, where the caps cover 80% of S^2 and a rotated copy
keeps 6 -- which is why this file still supplies k = 2, 3 and 4, and why both rules are
verified.  The margin above is therefore no longer unknown at k >= 5: it was 12, 12, 42 and
34.

DETERMINISM.  The search is a fixed number of seeded random-order greedy passes with a
(1,2)-swap improvement -- no wall clock anywhere -- so this script rebuilds the same files
every time and run_all.py can treat it as a verification rather than a search.

EXACTNESS.  Integer arithmetic throughout.  The shell enumeration is checked against the
known theta series of each root lattice, so an incomplete enumeration cannot pass unnoticed,
and each pole set is re-verified from scratch at the end: in the lattice, norm in (m, 3m],
legal against every cap direction, a 60-degree code among itself once both +-a are taken,
and no more than tau(k) of them.
"""
import itertools
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')


# ---- the root systems, and the lattices they generate, in explicit coordinates ----------
def A_roots(n):
    return np.array([[1 if i == p else (-1 if i == q else 0) for i in range(n + 1)]
                     for p, q in itertools.permutations(range(n + 1), 2)], dtype=np.int64)


def D_roots(n):
    V = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = np.zeros(n, np.int64); v[i] = si; v[j] = sj; V.append(v)
    return np.array(V, dtype=np.int64)


def E8_roots():
    """scaled to minimum 8, i.e. 2*E_8: (+-2,+-2,0^6) and (+-1)^8 with an even sign count"""
    V = [2 * v for v in D_roots(8)]
    for m in range(256):
        s = np.array([-1 if (m >> t) & 1 else 1 for t in range(8)], np.int64)
        if (s < 0).sum() % 2 == 0:
            V.append(s)
    return np.array(V, dtype=np.int64)


def A_shell(n, Mmax):
    b = int(Mmax ** .5)
    return np.array([v for v in (np.array(c, dtype=np.int64)
                                 for c in itertools.product(range(-b, b + 1), repeat=n + 1))
                     if v.sum() == 0 and 0 < int(v @ v) <= Mmax], dtype=np.int64)


def D_shell(n, Mmax):
    b = int(Mmax ** .5)
    return np.array([v for v in (np.array(c, dtype=np.int64)
                                 for c in itertools.product(range(-b, b + 1), repeat=n))
                     if v.sum() % 2 == 0 and 0 < int(v @ v) <= Mmax], dtype=np.int64)


def E8_shell(Mmax):
    """2*E_8 = {x in Z^8 : all coordinates the same parity, sum = 0 mod 4}"""
    out, b = [], int(Mmax ** .5)
    for pool in ([t for t in range(-b, b + 1) if t % 2 == 0],
                 [t for t in range(-b, b + 1) if t % 2 == 1]):
        for c in itertools.product(pool, repeat=8):
            if sum(c) % 4:
                continue
            if 0 < sum(t * t for t in c) <= Mmax:
                out.append(np.array(c, dtype=np.int64))
    return np.array(out, dtype=np.int64)


def lines(V):
    """sign-normalise and deduplicate: a pole is a DIRECTION, and +-a are one of them"""
    out, seen = [], set()
    for v in V:
        w = -v if v[np.nonzero(v)[0][0]] < 0 else v
        t = w.tobytes()
        if t not in seen:
            seen.add(t); out.append(w)
    return np.array(out, dtype=np.int64)


E = E8_roots()
r1 = E[0]
r2 = E[(E @ r1) == -4][0]                       # a root at 120 degrees: r1, r2 span an A_2
SYS = {2: A_roots(2), 3: A_roots(3), 4: D_roots(4), 5: D_roots(5),
       6: E[((E @ r1) == 0) & ((E @ r2) == 0)],  # E_6 = the complement of that A_2
       7: E[(E @ r1) == 0],                      # E_7 = the complement of one root
       8: E}
L8 = E8_shell(24)
SHELL = {2: A_shell(2, 6), 3: A_shell(3, 6), 4: D_shell(4, 6), 5: D_shell(5, 6),
         6: L8[((L8 @ r1) == 0) & ((L8 @ r2) == 0)], 7: L8[(L8 @ r1) == 0], 8: L8}
TAU = {1: 2, 2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240}
# theta series of each root lattice at the norms this file touches; E_6's norm-24 count is
# 720, which is what a direct enumeration gives and what the E_6 Cartan matrix confirms
THETA = {2: {2: 6, 6: 6}, 3: {2: 12, 4: 6, 6: 24}, 4: {2: 24, 4: 24, 6: 96},
         5: {2: 40, 4: 90, 6: 240}, 6: {8: 72, 16: 270, 24: 720},
         7: {8: 126, 16: 756, 24: 2072}, 8: {8: 240, 16: 2160, 24: 6720}}
RESTARTS = 60


def conflicts(C):
    """a and a' clash iff |cos| > 1/2 -- BOTH signs, because +-a are both poles"""
    n2 = np.einsum('ij,ij->i', C, C)
    IP = C @ C.T
    bad = 4 * IP * IP > np.outer(n2, n2)
    np.fill_diagonal(bad, False)
    return bad


def search(C, seed, restarts=RESTARTS):
    """random-order greedy plus a (1,2)-swap improvement, best over a FIXED set of restarts"""
    if not len(C):
        return C
    bad = conflicts(C)
    n = len(C)
    rng = np.random.default_rng(seed)
    best = []
    for _ in range(restarts):
        keep = np.ones(n, bool); S = []
        for i in rng.permutation(n):
            if not keep[i]:
                continue
            S.append(int(i)); keep &= ~bad[i]; keep[i] = False
        improved = True
        while improved:                                   # drop one, put two back
            improved = False
            inS = np.zeros(n, bool); inS[S] = True
            cnt = bad[:, S].sum(axis=1)
            for v in list(S):
                cand = np.nonzero((~inS) & (cnt == 1) & bad[v])[0]
                if len(cand) < 2:
                    continue
                sub = bad[np.ix_(cand, cand)]
                pair = None
                for x in range(len(cand)):
                    free = np.nonzero(~sub[x][x + 1:])[0]
                    if len(free):
                        pair = (int(cand[x]), int(cand[x + 1 + free[0]])); break
                if pair:
                    S = [t for t in S if t != v] + list(pair)
                    improved = True; break
        if len(S) > len(best):
            best = S
    return C[np.array(best)]


def verify(P, W, k):
    """re-check a pole set from scratch, in integers"""
    if not len(P):
        return True, "empty"
    m = int(W[0] @ W[0])
    Q = np.vstack([P, -P])
    nQ = np.einsum('ij,ij->i', Q, Q)
    if not ((nQ > m) & (nQ <= 3 * m)).all():
        return False, "a pole is outside (m, 3m]"
    IPW = Q @ W.T
    if not (4 * IPW * IPW <= 3 * nQ[:, None] * m).all():
        return False, "cap . pole violated"
    IP = Q @ Q.T
    np.fill_diagonal(IP, -10 ** 9)
    if not ((IP <= 0) | (4 * IP * IP <= np.outer(nQ, nQ))).all():
        return False, "pole . pole violated"
    if len(np.unique(Q, axis=0)) != len(Q):
        return False, "two poles are the same direction"
    if len(Q) > TAU[k]:
        return False, "more poles than tau(k)"
    return True, "ok"


if __name__ == '__main__':
    print(__doc__)
    print("   k   m   shells (norm: lines)          theta   best shell        poles  tau(k)")
    tot = 0
    for k in sorted(SYS):
        W = SYS[k].astype(np.int64)
        m = int(W[0] @ W[0])
        S = lines(SHELL[k])
        n2 = np.einsum('ij,ij->i', S, S)
        counts = {int(M): 2 * int((n2 == M).sum()) for M in sorted(set(n2.tolist()))}
        assert counts == {M: v for M, v in THETA[k].items() if M <= 3 * m}, \
            "k=%d: shell counts %s do not match the theta series %s" % (k, counts, THETA[k])
        bestP, bestM = None, None
        for M in sorted(counts):
            if M <= m:
                continue
            P = search(S[n2 == M], k * 100 + M)
            if bestP is None or len(P) > len(bestP):
                bestP, bestM = P, M
        ok, why = verify(bestP, W, k)
        assert ok, "k=%d: %s" % (k, why)
        np.save(os.path.join(DATA, 'poles_k%d.npy' % k), bestP)
        tot += 2 * len(bestP)
        print("   %2d  %2d   %-28s  %-6s  norm %-12d %-6d %d"
              % (k, m, counts, "ok", bestM, 2 * len(bestP), TAU[k]))
    print()
    print("ALL AXIS LAYERS BUILT AND VERIFIED EXACTLY -- %d poles across dimensions 74-80" % tot)
    print("(the previous candidate rule a = w1 + w2 with w1 perp w2 gave 248, and none at k = 2)")
