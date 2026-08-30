"""Two-point bounds for the cos=1/2 graph on the 98280 Leech lines.
   Builds the 4-class association scheme from scratch, computes P, Q, the Delsarte LP
   bound (= Schrijver theta'), the Lovasz theta (LP without a_i>=0) and the ratio bound."""
import numpy as np, collections
from scipy.optimize import linprog

import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', '..', 'common'))
from leech import load as _leech          # noqa: E402
# One representative per line of the Leech minimal shell: 98280 lines at squared norm 32.
# Built from the Golay code rather than loaded, so nothing here is taken on trust.
_A = _leech().astype(np.int64)
_nz = _A != 0
_first = np.argmax(_nz, axis=1)
_sgn = np.sign(_A[np.arange(len(_A)), _first])
_A = _A * _sgn[:, None]
L = np.unique(_A, axis=0)
assert L.shape == (98280, 24), L.shape
N = len(L); print("lines", N)

# relation of a pair: 0 self, 1 |ip|=16, 2 |ip|=8, 3 ip=0
def rel_row(i):
    d = np.abs(L @ L[i])
    r = np.zeros(N, dtype=np.int8)
    r[d == 8] = 2; r[d == 16] = 1; r[d == 32] = 0; r[d == 0] = 3
    return r

r0 = rel_row(0)
k = [int((r0 == i).sum()) for i in range(4)]
print("valencies", k)

# intersection numbers p^h_{ij}: fix u=0, v with rel(u,v)=h, count w with rel(u,w)=i,rel(v,w)=j
rng = np.random.default_rng(0)
P_int = np.zeros((4,4,4), dtype=np.int64)
for h in range(4):
    idx = np.nonzero(r0 == h)[0]
    tab = None
    for trial in range(3):
        v = int(rng.choice(idx))
        rv = rel_row(v)
        t = np.zeros((4,4), dtype=np.int64)
        for i in range(4):
            for j in range(4):
                t[i,j] = int(((r0 == i) & (rv == j)).sum())
        if tab is None: tab = t
        else: assert (tab == t).all(), ("not a scheme at h=",h)
    P_int[h] = tab
print("intersection numbers verified constant over 3 random base pairs per relation")

# eigenmatrix: A_i has entries; use the intersection matrices B_i with (B_i)_{h,j} = p^h_{ij}
B = [P_int[:,i,:].T.copy() for i in range(4)]   # (B_i)[j][h] = p^h_{ij}
# common eigenvectors
Bsum = sum(np.random.default_rng(1).integers(1,50)*b for b in B).astype(float)
w, V = np.linalg.eig(Bsum)
order = np.argsort(-w.real)
V = V.real[:, order]
Pm = np.zeros((4,4))
for c in range(4):
    v = V[:, c]; v = v / v[0]
    for i in range(4):
        Pm[c, i] = (B[i] @ v)[0] / v[0] if False else (B[i] @ v)[0]
    # eigenvalue of A_i on this eigenspace = (B_i v)_0 / v_0 with v normalised v_0=1
Pm = np.round(Pm, 6)
print("P (rows = eigenspaces, cols = relations):\n", Pm)
mult = np.array([N / (sum(Pm[c,i]**2 / k[i] for i in range(4))) for c in range(4)])
print("multiplicities", np.round(mult, 4))
Q = np.zeros((4,4))
for i in range(4):
    for j in range(4):
        Q[i,j] = mult[j] * Pm[j,i] / k[i]
print("Q:\n", np.round(Q,4))

# ---------------- Delsarte LP = Schrijver theta'  (a_i >= 0) --------------------
# variables a_1,a_2,a_3 >= 0 with a_1 = 0 (relation 1 forbidden); a_0 = 1
# constraints sum_i a_i Q[i][j] >= 0  for j=1,2,3 ; maximise 1 + a_2 + a_3
free = [2,3]
c  = -np.ones(len(free))
A_ub = np.array([[-Q[i,j] for i in free] for j in range(1,4)])
b_ub = np.array([Q[0,j] for j in range(1,4)])
res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(0,None)]*len(free), method='highs')
print("Delsarte LP (=theta') bound: %.6f   -> at most %d lines" % (1-res.fun, int(np.floor(1-res.fun))))
print("   optimal inner distribution a =", np.round(res.x,3))

# ---------------- Lovasz theta (drop a_i >= 0) --------------------
res2 = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(None,None)]*len(free), method='highs')
print("Lovasz theta bound:            %.6f" % (1-res2.fun))

# ---------------- Hoffman ratio bound --------------------
ev = Pm[:,1]
print("eigenvalues of the conflict graph:", ev)
lmin = ev.min(); dgr = k[1]
print("ratio bound: %.6f" % (N * (-lmin) / (dgr - lmin)))
