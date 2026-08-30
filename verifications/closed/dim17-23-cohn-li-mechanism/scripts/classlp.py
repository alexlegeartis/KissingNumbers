"""A RIGOROUS upper bound on alpha(G_0), via the Clebsch classes and sumset growth.

If this bound comes out at 160, dimension 19 is CLOSED at tau(19) = 11948, because
maxbip(G_0) <= 2*alpha(G_0).

Setup (section 71).  K_0 = N + W with N the common kernel of the five (-12)-eigenvalue
characters (|N| = 32) and K_0/N the Clebsch graph on 16 classes.  Fix a LINEAR section
r: K_0/N -> K_0.  An independent set A splits as A_c over the classes; writing
Atil_c = A_c + r_c inside N, the independence condition on the Clebsch edge (c, c') with
label i = c + c' is exactly

    Atil_c  n  ( Atil_{c'} + Ghat_i )  =  empty,      Ghat_i = G_i + r_{v_i} in N, |Ghat_i| = 4,

(using linearity of r), and since both sides live in N,

    |Atil_c|  +  |Atil_{c'} + Ghat_i|  <=  32.

So the whole problem is controlled by SUMSET GROWTH in N = F_2^5.  Define

    f_i(m) = min { |X + Ghat_i| : X in N, |X| = m },

computed here EXACTLY by CP-SAT for every i and every m from 0 to 32 (a 64-variable model, so
this is cheap and not an estimate).  Then every independent set satisfies the integer program

    maximise  sum_c a_c    subject to    a_c + f_i(a_{c'}) <= 32   for every Clebsch edge,

whose optimum is a valid upper bound on alpha(G_0).  Ho's 320 corresponds to a_c = 32 on an
independent 5-set of Clebsch classes and 0 elsewhere, which is feasible here, so the bound is
at least 160 by construction; the question is whether it is more.

    python classlp.py
"""
import sys, os, itertools, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'common'))
import numpy as np
from capcheck import golay
from ortools.sat.python import cp_model

GOL = np.array(sorted(golay()), dtype=np.uint8)[:, 5:]
W_ = GOL.sum(1)
codes = [int("".join(map(str, r)), 2) for r in GOL]
Sraw = sorted(codes[i] for i in range(len(codes)) if 0 < W_[i] <= 4)
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
print(__doc__)

chis = [c for c in range(NV)
        if sum(-1 if bin(c & s).count('1') & 1 else 1 for s in S0) == -12]
cls = lambda v: sum(((bin(v & c).count('1') & 1) << i) for i, c in enumerate(chis))
N = sorted(v for v in range(NV) if cls(v) == 0)
Nidx = {v: i for i, v in enumerate(N)}
IMG = sorted({cls(v) for v in range(NV)})
print("|N| = %d ; classes = %d" % (len(N), len(IMG)))

# linear section r: pick preimages of a basis of the image and extend linearly
qb, qsp, pre = [], {0}, {}
for v in range(NV):
    c = cls(v)
    if c not in qsp:
        qb.append(c); pre[c] = v
        qsp = {x ^ y for x in qsp for y in (0, c)}
R = {0: 0}
for c in qb:
    R.update({k ^ c: R[k] ^ pre[c] for k in list(R)})
assert set(R) == set(IMG) and all(cls(R[c]) == c for c in IMG)
assert all(R[a] ^ R[b] == R[a ^ b] for a in IMG for b in IMG), "section is not linear"

blocks = collections.defaultdict(list)
for s in S0:
    blocks[cls(s)].append(s)
GH = {vi: sorted(Nidx[g ^ R[vi]] for g in blk) for vi, blk in blocks.items()}
print("the five Ghat_i inside N: %s" % GH)

XOR = [[Nidx[N[i] ^ N[j]] for j in range(32)] for i in range(32)]


def growth(gh):
    """f(m) = min |X + Ghat| over |X| = m, exactly."""
    out = {}
    for m in range(33):
        mod = cp_model.CpModel()
        x = [mod.NewBoolVar("x%d" % i) for i in range(32)]
        y = [mod.NewBoolVar("y%d" % i) for i in range(32)]
        mod.Add(sum(x) == m)
        for i in range(32):
            for g in gh:
                mod.AddImplication(x[i], y[XOR[i][g]])
        mod.Minimize(sum(y))
        s = cp_model.CpSolver()
        s.parameters.max_time_in_seconds = 10.0
        s.parameters.num_search_workers = 4
        st = s.Solve(mod)
        assert st == cp_model.OPTIMAL, (m, s.StatusName(st))
        out[m] = int(s.ObjectiveValue())
    return out


F = {}
for vi, gh in sorted(GH.items()):
    F[vi] = growth(gh)
    print("  f_%d = %s" % (vi, [F[vi][m] for m in range(33)]), flush=True)

mod = cp_model.CpModel()
a = {c: mod.NewIntVar(0, 32, "a%d" % c) for c in IMG}
sel = {c: [mod.NewBoolVar("s%d_%d" % (c, m)) for m in range(33)] for c in IMG}
for c in IMG:
    mod.AddExactlyOne(sel[c])
    mod.Add(a[c] == sum(m * sel[c][m] for m in range(33)))
edges = [(c, c ^ vi, vi) for c in IMG for vi in GH if (c ^ vi) in R and c < (c ^ vi)]
print("\nClebsch edges: %d" % len(edges))
for c, cp, vi in edges:
    for m in range(33):
        mod.Add(a[c] + F[vi][m] <= 32).OnlyEnforceIf(sel[cp][m])
        mod.Add(a[cp] + F[vi][m] <= 32).OnlyEnforceIf(sel[c][m])
mod.Maximize(sum(a.values()))
s = cp_model.CpSolver()
s.parameters.max_time_in_seconds = 120.0
s.parameters.num_search_workers = 8
st = s.Solve(mod)
bound = int(s.ObjectiveValue())
print("\nRIGOROUS upper bound: alpha(G_0) <= %d   (status %s)" % (bound, s.StatusName(st)))
print("profile: %s" % sorted((s.Value(a[c]) for c in IMG), reverse=True))
print("=> maxbip(G_0) <= %d ; Ho's record is 320" % (2 * bound))
if 2 * bound <= 320:
    print("*** DIMENSION 19 IS CLOSED: tau(19) = 11948 ***")
