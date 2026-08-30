"""alpha(G_0) exactly -- the number that decides dimension 19.

maxbip(G_0) <= 2*alpha(G_0), and Ho's record is maxbip(G_0) = 320.  So alpha(G_0) = 160 closes
dimension 19 at tau(19) = 11948, while alpha(G_0) >= 161 is the only door left open.
alpha(G_1) = 80 was proven optimal on the 256-vertex quotient, giving alpha(G_0) >= 160; the
lift cannot do better because |A|, |B| <= alpha(G_1) caps that family at 2*160 = 320 exactly.

Three strengthenings make the 512-vertex model tractable where the plain one stalls:

  * vertex-transitivity: fix vertex 0 in the set;
  * maximality: every maximum independent set is maximal, so each vertex is in the set or has
    a neighbour in it;
  * girth 5: G_0 is triangle-free with no 4-cycles, so every 5-cycle carries at most two
    chosen vertices -- a large family of genuinely tight cuts that edge constraints miss.

    python alpha512.py [seconds] [max5cycles]
"""
import sys, os, collections, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'common'))
import numpy as np
from capcheck import golay
from ortools.sat.python import cp_model

TL = float(sys.argv[1]) if len(sys.argv) > 1 else 900.0
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
GOL = np.array(sorted(golay()), dtype=np.uint8)[:, 5:]
W = GOL.sum(1)
codes = [int("".join(map(str, r)), 2) for r in GOL]
Sraw = sorted(codes[i] for i in range(len(codes)) if 0 < W[i] <= 4)
even = [c for c in Sraw if bin(c).count('1') % 2 == 0]
sp, bas = {0}, []
for c in even:
    if c not in sp:
        bas.append(c); sp |= {x ^ c for x in sp}
D = len(bas); NV = 1 << D
idx = {}
for m in range(NV):
    v = 0
    for b in range(D):
        if m >> b & 1: v ^= bas[b]
    idx[v] = m
S0 = sorted(idx[c] for c in even)
S0set = set(S0)
print(__doc__)
print("G_0: %d vertices, %d-regular" % (NV, len(S0)))

m = cp_model.CpModel()
x = [m.NewBoolVar("x%d" % i) for i in range(NV)]
for v in range(NV):
    for s in S0:
        u = v ^ s
        if u > v:
            m.Add(x[v] + x[u] <= 1)
    if os.environ.get('NOMAX') != '1':
        m.Add(x[v] + sum(x[v ^ s] for s in S0) >= 1)      # maximality
if os.environ.get('NOSYM') != '1':
    m.Add(x[0] == 1)

# 5-cycles: v, v^a, v^a^b, v^a^b^c, v^a^b^c^d with a^b^c^d^e = 0, all in S0
n5 = 0
seen = set()
for a, b, c, d in itertools.combinations(S0, 4):
    for perm in ((a, b, c, d), (a, b, d, c), (a, c, b, d), (a, c, d, b), (a, d, b, c),
                 (a, d, c, b)):
        e = perm[0] ^ perm[1] ^ perm[2] ^ perm[3]
        if e in S0set:
            for v in range(NV):
                cyc = tuple(sorted((v, v ^ perm[0], v ^ perm[0] ^ perm[1],
                                    v ^ perm[0] ^ perm[1] ^ perm[2], v ^ e)))
                if len(set(cyc)) == 5 and cyc not in seen:
                    seen.add(cyc)
                    m.Add(sum(x[t] for t in cyc) <= 2)
                    n5 += 1
                    if n5 >= CAP:
                        break
            if n5 >= CAP:
                break
    if n5 >= CAP:
        break
print("5-cycle cuts added: %d" % n5, flush=True)

m.Maximize(sum(x))
s = cp_model.CpSolver()
s.parameters.max_time_in_seconds = TL
s.parameters.num_search_workers = 8
s.parameters.log_search_progress = False
st = s.Solve(m)
val, ub = int(s.ObjectiveValue()), int(s.BestObjectiveBound())
print("\nalpha(G_0): best %d, proven upper bound %d, status %s" % (val, ub, s.StatusName(st)))
print("=> maxbip(G_0) <= 2*%d = %d   (Ho's record is 320)" % (ub, 2 * ub))
if 2 * ub <= 320:
    print("*** DIMENSION 19 IS CLOSED: tau(19) = 11948 is the ceiling of this construction ***")
elif val > 160:
    print("*** alpha(G_0) > 160 -- the door to a new dimension-19 record is open ***")
