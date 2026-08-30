"""alpha with every fibre capped at 22 -- the only regime where 193 could live.

The heavy model (61440 vertex 4-clique cuts) spends its budget building rather than searching.
The fibre cuts dominate them anyway: for any 4 MUTUALLY ADJACENT fibres the exact independence
number is 32, and Q has 120 such 4-cliques, so 120 constraints replace 61440.  With
vertex-transitivity (x_0 = 1) and maximality that is the whole model.

Progress is printed as it comes so the bracket on alpha is visible while it runs.

    tau(17) = 5346 + 2*alpha ; alpha = 192 closes dimension 17, alpha >= 193 is a record.

    python cap22_det.py [seconds]
"""
import sys, os, itertools, collections
import numpy as np
from ortools.sat.python import cp_model

TL = float(sys.argv[1]) if len(sys.argv) > 1 else 2400.0
print(__doc__)
H = os.path.dirname(os.path.abspath(__file__))
SUP = [tuple(int(v) for v in r) for r in np.load(os.path.join(H, '..', 'data', 'd17_supports.npy'))]
rm1 = {0}
for f in range(1, 16):
    for c in (0, 1):
        v = sum(1 << x for x in range(16) if bin(x & f).count('1') % 2 == c)
        rm1 |= {y ^ v for y in rm1}
RM2 = [x for x in range(1 << 16) if all(bin(x & y).count('1') % 2 == 0 for y in rm1)]
m0 = sum(1 << i for i in SUP[0])
Kl = sorted(b for b in RM2 if bin(b & m0).count('1') % 2 == 0)
idx = {b: i for i, b in enumerate(Kl)}
N = len(Kl)
W4 = [b for b in Kl if bin(b).count('1') == 4]
W4s = set(W4)
fibid, reps = {}, []
for b in Kl:
    if b not in fibid:
        for r in rm1:
            fibid[b ^ r] = len(reps)
        reps.append(b)
Q = len(reps)
qadj = [set(d for d in range(Q) if d != c and
            any((reps[c] ^ reps[d] ^ r) in W4s for r in rm1)) for c in range(Q)]
q4 = [cl for cl in itertools.combinations(range(Q), 4)
      if all(cl[j] in qadj[cl[i]] for i in range(4) for j in range(i + 1, 4))]
print("graph %d vertices, %d-regular ; quotient %d fibres with %d 4-cliques"
      % (N, len(W4), Q, len(q4)), flush=True)

m = cp_model.CpModel()
x = [m.NewBoolVar("x%d" % i) for i in range(N)]
adj = [[idx[Kl[i] ^ w] for w in W4] for i in range(N)]
for i in range(N):
    for j in adj[i]:
        if j > i:
            m.Add(x[i] + x[j] <= 1)
    m.Add(x[i] + sum(x[j] for j in adj[i]) >= 1)
byfib = collections.defaultdict(list)
for i, b in enumerate(Kl):
    byfib[fibid[b]].append(i)
a = [m.NewIntVar(0, 32, "a%d" % c) for c in range(Q)]
for c in range(Q):
    m.Add(a[c] == sum(x[i] for i in byfib[c]))
for cl in q4:
    m.Add(sum(a[c] for c in cl) <= 32)
m.Add(x[0] == 1)
# PROVEN (section 74e/74f): a_c >= 29 forces every neighbour of c to zero, so the whole set is
# at most 29 + alpha(G[R(c)]) = 189.  Hence any set larger than 192 has every a_c <= 28.
# Tightened: the exact ILP count of block-translates fitting in the complement of A_c gives
#   a_c + P(32-a_c) + 160 <= 192 for every a_c >= 23, so a 193-set has every a_c <= 22.
for c in range(Q):
    m.Add(a[c] <= 22)
m.Maximize(sum(x))


class CB(cp_model.CpSolverSolutionCallback):
    def __init__(self):
        super().__init__()
        self.n = 0

    def on_solution_callback(self):
        self.n += 1
        print("   solution %d: |S| = %d   bound %d   (%.0fs)"
              % (self.n, int(self.ObjectiveValue()), int(self.BestObjectiveBound()),
                 self.WallTime()), flush=True)


s = cp_model.CpSolver()
s.parameters.max_time_in_seconds = TL
s.parameters.num_search_workers = 1   # single worker => DETERMINISTIC, so the
s.parameters.random_seed = 1          # dual bound below is reproducible exactly
st = s.Solve(m, CB())
best, ub = int(s.ObjectiveValue()), int(s.BestObjectiveBound())
print("\nalpha = %d, proven upper bound %d, status %s" % (best, ub, s.StatusName(st)))
print("tau(17): achieved %d, ceiling %d   (record 5730)" % (5346 + 2 * best, 5346 + 2 * ub))
if best > 192:
    print("*** NEW RECORD IN DIMENSION 17 ***")
    np.save(os.path.join(H, '..', 'data', 'd17_B_new.npy'),
            np.array([Kl[i] for i in range(N) if s.Value(x[i])], dtype=np.int64))
elif ub <= 192:
    print("*** 192 IS OPTIMAL -- dimension 17 CLOSED at tau(17) = 5730 ***")
