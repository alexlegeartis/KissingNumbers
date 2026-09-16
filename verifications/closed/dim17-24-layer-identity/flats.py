"""THE FLAT LAYER, dimensions 19 to 23 -- what it is, how big it can be, and where it dies.

    python flats.py

The second family in Cohn-Li (arXiv:2411.04916) is, in R^n with 17 <= n <= 23,

    BASE   4*C(n,2) vectors  (+-2,+-2,0^{n-2}),
           128*|C| vectors  (+-1)^8 on the support of a block, ODD number of minus signs,
           C = the weight-8 words of the Golay code SHORTENED at the 24-n deleted positions;
    FLATS  sqrt(8/n) * (+-1)^n.

Everything has norm^2 = 8, so the 60-degree condition is "pairwise inner product <= 4".
Three facts, each derived here rather than assumed:

  * a flat with minus-signs on c meets a block O in  sqrt(8/n)(8-2k),  k = #disagreements,
    and the block carries EVERY odd pattern, so k = 0 happens unless <c,O> = 0.  The blocks
    span the shortened Golay code; the Golay is self-dual and shortening-then-dualising is
    puncturing; so c ranges over the (24-n)-PUNCTURED Golay code, 4096 words.
  * two flats at Hamming distance d meet in (8/n)(n-2d) <= 4  iff  d >= n/4.
  * so the flat layer is a maximum independent set in Cay(F_2^12, S), S = the punctured
    Golay's nonzero words of weight < ceil(n/4).

This script computes S, settles alpha EXACTLY where the structure settles it, and then
BUILDS the dimension-20 and dimension-21 configurations and checks every pair.

Two things worth being precise about, because they are easy to state wrongly:

  * BASE = tau(Lambda_n) only for n = 19, 20, 21.  In 17, 18, 22 and 23 this family's base is
    strictly smaller than the laminated lattice (4384/3744, 6500/6372, 43164, 65780 against
    5346, 7398, 49896, 93150), so in those dimensions the record is NOT base + flats.  For 17
    and 18 it is the OTHER Cohn-Li family, R^16 (+) R^k over the odd Barnes-Wall lattice; for
    22 and 23 it is the lattice itself.
  * for n <= 18 the size of C depends on WHICH coordinates are deleted -- the Golay octads are
    a 5-design, so deleting up to 5 points gives a count independent of the choice (210, 130,
    78 for n = 21, 20, 19) and deleting 6 or 7 does not (45 or 46, 25 or 30).
"""
import os, sys, itertools, collections

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', '..', 'common'))
from golay import golay24

ok_all = []


def ok(cond, msg):
    ok_all.append(bool(cond))
    print(("  OK   " if cond else "  FAIL ") + msg)
    return bool(cond)


G = golay24().astype(np.uint8)
ok(G.shape == (4096, 24), "the Golay code: 4096 words")
ok(dict(zip(*[x.tolist() for x in np.unique(G.sum(1), return_counts=True)]))
   == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}, "weight distribution 1/759/2576/759/1")


def f2_rank(M):
    M = M.copy() % 2
    r = 0
    for c in range(M.shape[1]):
        p = next((i for i in range(r, len(M)) if M[i, c]), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        for i in range(len(M)):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
    return r


TAU_LAMBDA = {17: 5346, 18: 7398, 19: 10668, 20: 17400, 21: 27720,
              22: 49896, 23: 93150}
RECORD = {19: 11948, 20: 19448, 21: 29768, 22: 49896, 23: 93150}

print("\nthe base of this family, and how it compares with the laminated lattice:")
print("   n   |C|   4C(n,2)   base    tau(Lambda_n)   base is")
BASE = {}
for n in range(17, 24):
    P = 24 - n
    blocks = G[(G.sum(1) == 8) & (G[:, :P].sum(1) == 0)][:, P:]
    base = 4 * (n * (n - 1) // 2) + 128 * len(blocks)
    BASE[n] = (base, blocks)
    print("   %2d  %4d   %6d  %6d   %10d    %s"
          % (n, len(blocks), 4 * (n * (n - 1) // 2), base, TAU_LAMBDA[n],
             "EQUAL" if base == TAU_LAMBDA[n] else "short by %d" % (TAU_LAMBDA[n] - base)))
ok(all(BASE[n][0] == TAU_LAMBDA[n] for n in (19, 20, 21)),
   "base = tau(Lambda_n) exactly for n = 19, 20, 21")
ok(all(BASE[n][0] < TAU_LAMBDA[n] for n in (17, 18, 22, 23)),
   "and is strictly smaller for n = 17, 18, 22, 23")

# For n <= 18 the size of C depends on WHICH coordinates are deleted, because the Golay
# octads are only a 5-design.  Checked over EVERY deletion set, not sampled.
OCT = G[G.sum(1) == 8].astype(np.int16)
for n, P, want in ((21, 3, {210}), (20, 4, {130}), (19, 5, {78}), (18, 6, {45, 46}), (17, 7, {25, 30})):
    seen = set()
    subs = list(itertools.combinations(range(24), P))
    for s in range(0, len(subs), 20000):
        blk = subs[s:s + 20000]
        M = np.zeros((len(blk), 24), dtype=np.int16)
        for r, D in enumerate(blk):
            M[r, list(D)] = 1
        ov = M @ OCT.T
        seen |= set(int(x) for x in (ov == 0).sum(1))
    ok(seen == want, "   n = %2d: over all %d deletion sets |C| takes exactly the values %s"
       % (n, len(subs), sorted(seen)))
print("   so for n >= 19 the count is independent of the choice (5-design), and for n <= 18 it")
print("   is not -- which is why a shortened-Golay base must say which positions it deleted.")

# ---------------------------------------------------------------- the ambient and S
print("\nthe flat layer as an independent set:")
STRUCT = {}
for n in (19, 20, 21, 22, 23):
    P = 24 - n
    amb = np.unique(G[:, P:], axis=0)
    blocks = BASE[n][1]
    need = -(-n // 4)
    S = np.array([w for w in amb if 0 < int(w.sum()) < need], dtype=np.uint8)
    rk_blocks = f2_rank(blocks.astype(np.int8))
    STRUCT[n] = (amb, S, need)
    print("   n = %2d: ambient = the %d-punctured Golay, %d words, d_min = %d; blocks span a "
          "[%d,%d] code" % (n, P, len(amb), min(int(w.sum()) for w in amb if w.sum()), n, rk_blocks))
    print("           flats need d >= %d; |S| = %d, weights %s"
          % (need, len(S), dict(zip(*[x.tolist() for x in np.unique(S.sum(1), return_counts=True)]))
             if len(S) else {}))
    ok(len(amb) == 4096 and rk_blocks == 12 - P,
       "   n = %2d: ambient has 4096 words and the blocks span the shortened code (dim %d)"
       % (n, 12 - P))


def cayley(amb, S):
    """Components, degree and bipartiteness of Cay(ambient, S).  Exact, by BFS."""
    idx = {tuple(w): i for i, w in enumerate(amb)}
    nb = [[] for _ in amb]
    for i, w in enumerate(amb):
        for s in S:
            j = idx.get(tuple(np.bitwise_xor(w, s)))
            if j is not None:
                nb[i].append(j)
    colour = [-1] * len(amb)
    comps, bip = [], True
    for start in range(len(amb)):
        if colour[start] >= 0:
            continue
        colour[start] = 0
        stack, comp = [start], [start]
        while stack:
            u = stack.pop()
            for v in nb[u]:
                if colour[v] < 0:
                    colour[v] = 1 - colour[u]
                    stack.append(v)
                    comp.append(v)
                elif colour[v] == colour[u]:
                    bip = False
        comps.append(comp)
    deg = set(len(x) for x in nb)
    return comps, deg, bip, colour


print("\ndimension 20 and dimension 21: alpha is forced, and equals Cohn-Li's value")
FLAT = {}
for n in (20, 21):
    amb, S, need = STRUCT[n]
    comps, deg, bip, colour = cayley(amb, S)
    sizes = collections.Counter(len(c) for c in comps)
    print("   n = %2d: %d components, sizes %s; degrees %s; bipartite %s"
          % (n, len(comps), dict(sizes), sorted(deg), bip))
    ok(bip and deg == {len(S)}, "   n = %2d: the graph is %d-regular and BIPARTITE, so every "
                                "component is regular bipartite and has a perfect matching; by "
                                "Koenig alpha = |V|/2 = 2048" % (n, len(S)))
    # both sides are independent; take the one containing 0 and check the halves are equal
    side = [i for i in range(len(amb)) if colour[i] == 0]
    half_ok = all(sum(1 for i in c if colour[i] == 0) == len(c) // 2 for c in comps)
    ok(len(side) == 2048 and half_ok,
       "   n = %2d: an explicit maximum flat layer, 2048 words, one side of each component"
       % n)
    F = amb[side]
    D = (F[:, None, :] != F[None, :, :]).sum(2)
    np.fill_diagonal(D, n)
    ok(int(D.min()) >= need, "   n = %2d: every pair of those flats is at distance >= %d"
       % (n, need))
    FLAT[n] = F

print("\n   n = 19: |S| = %d has words of weight 3 AND 4, so the graph is NOT bipartite; its"
      % len(STRUCT[19][1]))
print("           four components of 1024 have alpha = 320 each, 4 x 320 = 1280, which is")
print("           certified in verifications/closed/dim17-23-cohn-li-mechanism/ (Ho's 11948).")

print("\n   n = 22, 23: S is EMPTY -- the punctured Golay already has d_min >= ceil(n/4) -- so")
print("           the whole ambient is a flat layer, 4096 words, and yet:")
for n in (22, 23):
    tot = BASE[n][0] + 4096
    ok(tot < TAU_LAMBDA[n], "   n = %2d: %d + 4096 = %d < %d = tau(Lambda_%d), so this family "
       "loses to the lattice outright" % (n, BASE[n][0], tot, TAU_LAMBDA[n], n))

# ---------------------------------------------------------------- build 20 and 21, exactly
print("\nthe configurations themselves, every pair checked in integer arithmetic:")
for n in (20, 21):
    base_n, blocks = BASE[n]
    V = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (2, -2):
            for sj in (2, -2):
                v = [0] * n
                v[i] = si
                v[j] = sj
                V.append(v)
    for b in blocks:
        sup = np.nonzero(b)[0]
        for m in range(256):
            if bin(m).count('1') % 2 == 1:                     # ODD number of minus signs
                v = [0] * n
                for t, p in enumerate(sup):
                    v[p] = -1 if (m >> t) & 1 else 1
                V.append(v)
    U = np.array(V, dtype=np.int32)
    ok(len(U) == base_n and (int((U * U).sum(1).min()), int((U * U).sum(1).max())) == (8, 8),
       "   n = %2d: base built, %d vectors, all of norm^2 = 8" % (n, base_n))

    worst = -99
    for s in range(0, len(U), 2000):                            # chunked: the Gram is big
        Ch = U[s:s + 2000] @ U.T
        for t in range(len(Ch)):
            Ch[t, s + t] = -99
        worst = max(worst, int(Ch.max()))
    ok(worst <= 4, "   n = %2d: base is a 60-degree code (max inner product %d <= 4)"
       % (n, worst))
    ok(len(set(map(tuple, U.tolist()))) == base_n, "   n = %2d: base vectors all distinct" % n)

    Sg = np.where(FLAT[n] == 1, -1, 1).astype(np.int32)         # flat sign vectors, +-1
    # base . flat = sqrt(8/n) * (U . s);  <= 4  iff  (U.s) <= sqrt(2n)  (negatives are free)
    M = U @ Sg.T
    lim = int(np.floor((2 * n) ** 0.5))
    ok(int(M.max()) <= lim, "   n = %2d: base-against-flat, max U.s = %d <= floor(sqrt(2n)) = %d"
       % (n, int(M.max()), lim))
    print("   n = %2d: total = %d + %d = %d   (published record %d)"
          % (n, base_n, len(FLAT[n]), base_n + len(FLAT[n]), RECORD[n]))
    ok(base_n + len(FLAT[n]) == RECORD[n],
       "   n = %2d: the construction reaches the published record exactly" % n)

print()
print("So in 20 and 21 the flat layer is EXACTLY 2048 and Cohn-Li's numbers are the ceiling of")
print("this family; in 22 and 23 the family loses to the laminated lattice by 2636 and 23274.")
print()
if all(ok_all):
    print("ALL CHECKS PASSED -- %d of %d" % (sum(ok_all), len(ok_all)))
else:
    print("FAILURES: %d of %d checks failed" % (len(ok_all) - sum(ok_all), len(ok_all)))
    sys.exit(1)
