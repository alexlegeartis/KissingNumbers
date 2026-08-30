#!/usr/bin/env python3
"""Dimensions 73-80: the cap construction over Gamma_72 with EXPLICIT classes and
ZERO-SUM TRIPLES.

Two independent improvements over the first version of this package:

  (1) Explicit coordinates.  Gamma_72's Gram matrix and six automorphisms are published in
      the Catalogue of Lattices (Nebe, entry "Gamma72").  Verified here: 72x72, symmetric,
      diagonal all 8, determinant exactly 1 (Bareiss, integer arithmetic) -- even unimodular.
      Searching its minimal lines gives a class of 2106 lines where Caro-Wei could only
      guarantee 221, and 120 pairwise disjoint classes (2096-2106 each, 252 358 lines in
      total) come essentially free as random automorphism images of one class.

  (2) Zero-sum triples.  Gamma_72's class condition is |<u,u'>| <= 2 of 8, i.e. |cos| <= 1/4,
      which is exactly the threshold that lets the cap level drop from t = 3/4 to t = 2/3.
      At t = 2/3 the "same line, two cap directions" constraint becomes <z,z'> <= -1/2, so a
      class line may be used at a zero-sum TRIPLE of directions rather than an antipodal
      PAIR -- worth 4 points per line instead of 2.  Three is the maximum possible, since n
      points pairwise at cosine <= -1/2 satisfy n <= 1 + 1/(1/2) = 3.

      Gain = 2 * sum_i |C_i| * (|Z_i| - 1),  so a triple part is worth double a pair part.

      The parts are a partition of a kissing configuration W of R^k into zero-sum triples
      (plus leftovers); triples.py computes and verifies them for A_1, A_2, A_3, D_4, D_5,
      E_6, E_7, E_8, reaching the optimum ceil(|W|/3) in every case that allows it.

CALIBRATION.  The same scheme over the Leech, with its 248-line class, reproduces Cohn's
table in dimensions 26, 28, 29, 30, 31 exactly and beats it in 25 (+2) and 27 (+496) -- the
two records this project already holds.  See triples.py and capcheck2.py.

    python final.py
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N = 6218175600
TAU = {1: 2, 2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240}
SYS = {1: "A_1", 2: "A_2", 3: "A_3", 4: "D_4", 5: "D_5", 6: "E_6", 7: "E_7", 8: "E_8"}
# parts computed and verified by triples.py: (number of triples, number of pairs)
PARTS = {1: (0, 1), 2: (2, 0), 3: (4, 0), 4: (8, 0), 5: (12, 2), 6: (24, 0),
         7: (42, 0), 8: (80, 0)}

G = np.load(os.path.join(HERE, 'data', 'gamma72_gram.npy')).astype(np.int64)
Z = np.load(os.path.join(HERE, 'data', 'disjoint_classes.npz'))
keys = sorted(Z.files)
classes = [Z[k].astype(np.int64) for k in keys]

print(__doc__)
print("=" * 78)
print("VERIFICATION of the explicit data")
print("=" * 78)
assert G.shape == (72, 72) and (G == G.T).all() and (np.diag(G) == 8).all()


def det_bareiss(M):
    M = [[int(x) for x in r] for r in M]
    n = len(M)
    prev = 1
    sign = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


d = det_bareiss(G)
print("   Gram 72x72, symmetric, diagonal all 8, det = %d  ->  even unimodular: %s"
      % (d, d == 1))
seen = set()
worst = 0
for C in classes:
    nrm = np.einsum('ij,jk,ik->i', C, G, C)
    assert (nrm == 8).all()
    IP = C @ G @ C.T
    np.fill_diagonal(IP, 0)
    worst = max(worst, int(np.abs(IP).max()))
    for r in C:
        t = tuple(int(x) for x in r)
        assert t not in seen and tuple(-x for x in t) not in seen
        seen.add(t)
sizes = sorted((len(c) for c in classes), reverse=True)
print("   %d disjoint classes, sizes %d..%d, %d lines total, all of norm 8"
      % (len(classes), sizes[-1], sizes[0], sum(sizes)))
print("   largest |inner product| inside any class: %d  (a class needs <= 2): %s"
      % (worst, worst <= 2))
print("   all %d lines distinct and no antipodal collision: True" % len(seen))

print()
print("=" * 78)
print("RESULTS")
print("=" * 78)
print("  dim  k  W     scheme                     classes  gain        this work         previous")
out = []
# At t = 2/3 a pole (0,w) meets a cap at sqrt(1/3)*<z,w>, so a pole must sit at least 30
# degrees from EVERY cap direction.  Such directions do exist for k >= 2 (the root systems'
# min-max cosine is below sqrt(3)/2), but a full tau(k)-point pole configuration is only
# verified here for k = 2 (the 30-degree-rotated hexagon) and k = 4.  We therefore report
# the triple scheme WITHOUT poles -- conservative, and irrelevant at this scale since
# tau(k) <= 240 against gains in the hundreds of thousands.  For k = 1 no pole direction
# exists at all at t = 2/3, so there we use the pair scheme at t = 3/4, where b = 1/2 makes
# poles automatic.
LAM = {1: 1, 2: 3, 3: 6, 4: 12, 5: 20, 6: 36, 7: 63, 8: 120}
for k in range(1, 9):
    tri, pair = PARTS[k]
    need = tri + pair
    assert need <= len(sizes), (k, need)
    use = sizes[:need]
    g_tri = 4 * sum(use[:tri]) + 2 * sum(use[tri:tri + pair])
    # The axis layer, in whichever form the package ships for this k.  Two rules produce
    # one, and each dimension takes the larger:
    #   shell  (data/poles_k<k>.npy) -- every norm-2m direction of the same root lattice is
    #          automatically at >= 30 degrees from every cap direction (scripts/poles18.py);
    #   rotated (data/poles_rot_k<k>.npz, 2026-08-27) -- a rational G-isometry of the root
    #          system itself, which inherits the 60-degree condition and only has to clear
    #          the 30-degree caps (scripts/axis_rotate.py).
    # The rotation wins from k = 5 up, and at k = 3 as well once the search is a DESCENT
    # rather than sampling: sampling kept 6 of the 12 there, against the shell rule's 8, and
    # descent reaches all 12, which is tau(3).  It ties at k = 4, where the shell rule
    # already attains tau(4) = 24.  Only the files that beat the shell rule exist at all, so
    # the branch below is the whole of the choice.  Both count +-a.
    _rf = os.path.join(HERE, 'data', 'poles_rot_k%d.npz' % k)
    _pf = os.path.join(HERE, 'data', 'poles_k%d.npy' % k)
    if os.path.exists(_rf):
        _np = int(len(np.load(_rf)['keep']))
    else:
        _np = 2 * len(np.load(_pf)) if os.path.exists(_pf) else 0
    g_tri += _np          # t = 2/3, no poles
    g_pair = 2 * sum(sizes[:LAM[k]]) + TAU[k]                          # t = 3/4, with poles
    if g_tri >= g_pair:
        gain, how = g_tri, "%2d tri + %d pair + %d poles, t=2/3" % (tri, pair, _np)
    else:
        gain, how = g_pair, "%2d pairs + %d poles, t=3/4" % (LAM[k], TAU[k])
    new = N + gain
    old = N + TAU[k]
    out.append((72 + k, new, old, gain, how, SYS[k], _np))
    print("   %2d  %d  %-4s %-26s %4d   +%-9d %-16d  %d"
          % (72 + k, k, SYS[k], how, need, gain, new, old))
print()
print("   For comparison, the Caro-Wei + antipodal-pairs version of this package gave")
print("   gains of 442, 1326, 2652, 5304, 8840, 15912, 27846, 53040.")
print("   Dimension 80: %d vs 53040, a factor %.2f."
      % (out[-1][3], out[-1][3] / 53040))
print()
prev = 0
for d_, new, old, g, _h, _w, _p in out:
    assert new > prev
    prev = new
print("   monotone in the dimension: True")

# Dimension 80 is the floor that final81-95.py asserts monotonicity against.  It used to be a
# constant there and had drifted 98 low -- the axis layer that verify_poles.py certifies was
# added here but not copied across -- so it is handed over rather than repeated.
import json
json.dump({'dim80': int(out[-1][1]),
           'rows': [[int(d_), int(new), int(old), int(g)]
                    for d_, new, old, g, _h, _w, _p in out],
           # scripts/make_readme_table.py builds the 73-80 rows of README.md from these,
           # the way it already builds 81-95 from rows81-95.json.  The hand-written block
           # had drifted a pole layer behind the driver.
           'table': [{'dim': int(d_), 'k': int(d_) - 72, 'W': _w, 'scheme': _h,
                      'gain': int(g), 'value': int(new), 'previous': int(old),
                      'poles': int(_p)}
                     for d_, new, old, g, _h, _w, _p in out]},
          open(os.path.join(HERE, 'rows73-80.json'), 'w'), indent=1)
print("   wrote rows73-80.json (dimension 80 = %d), read by final81-95.py" % out[-1][1])
