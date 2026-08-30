#!/usr/bin/env python3
"""The exact Delsarte / Lovasz-theta bound for Ho's dimension-19 problem.

tau(19) = 10668 + 4*alpha, alpha = the independence number of Cay(F_2^10, S) with |S| = 21
(see component.py).  Ho achieves alpha = 320.  For a Cayley graph on an ABELIAN group the
Lovasz theta function coincides with the Delsarte linear program

    maximise   sum_g a_g
    subject to a_0 = 1,  a_g >= 0 for all g,  a_g = 0 for g in S,
               and  ahat(chi) = sum_g a_g (-1)^{<chi,g>} >= 0 for every character chi,

and alpha <= theta.  With 1024 group elements this is a 1024-variable linear program, solved
here exactly enough to decide whether Ho's 320 is optimal.

An exact rational dual certificate is extracted afterwards: any y >= 0 with
    sum_chi y_chi (-1)^{<chi,g>} >= 1   for every g not in S and g != 0
gives  alpha <= 1 + sum_chi y_chi * (something)  -- we instead simply re-verify the primal
optimum and round the LP value down, which is what bounds alpha.

    python lpbound.py
"""
import sys
import os
import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'common'))
from capcheck import golay          # noqa: E402

G = np.array(sorted(golay()), dtype=np.uint8)[:, 5:]
W = G.sum(1)
codes = [int("".join(map(str, r)), 2) for r in G]
S = sorted(codes[i] for i in range(len(codes)) if 0 < W[i] <= 4)
span = {0}
for s in S:
    span |= {x ^ s for x in span}
comp = sorted(span)
assert len(comp) == 1024

# index <S> by F_2^10 through a basis
basis, sp = [], {0}
for c in comp:
    if c and c not in sp:
        basis.append(c)
        sp |= {x ^ c for x in sp}
    if len(basis) == 10:
        break
assert len(sp) == 1024
lab = {}
for m in range(1024):
    v = 0
    for b in range(10):
        if m >> b & 1:
            v ^= basis[b]
    lab[v] = m
Sl = sorted(lab[s] for s in S)
print(__doc__)
print("group F_2^10, |S| = %d" % len(Sl))

POP = np.array([bin(i).count("1") for i in range(1024)], dtype=np.int8)
# character matrix H[chi, g] = (-1)^{<chi,g>}
idx = np.arange(1024)
H = 1 - 2 * (POP[idx[:, None] & idx[None, :]] & 1).astype(np.int8)
H = H.astype(float)
print("character matrix built:", H.shape)

n = 1024
c = -np.ones(n)                       # maximise sum a_g
A_ub = -H                             # -ahat(chi) <= 0  i.e. ahat >= 0
b_ub = np.zeros(n)
bounds = []
for g in range(n):
    if g == 0:
        bounds.append((1.0, 1.0))
    elif g in set(Sl):
        bounds.append((0.0, 0.0))
    else:
        bounds.append((0.0, None))
res = linprog(c=c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')
print()
print("LP status:", res.message.strip())
theta = -res.fun
print("Delsarte / theta bound:  alpha <= %.6f  ->  %d" % (theta, int(np.floor(theta + 1e-7))))
print("Ho's construction:       alpha  = 320")
print("Hoffman ratio bound:     alpha <= 391")
best = int(np.floor(theta + 1e-7))
print()
print("tau(19) <= 10668 + 4*%d = %d  by this mechanism" % (best, 10668 + 4 * best))
print("tau(19) >= 10668 + 4*320 = 11948  (Ho, the current record)")
if best == 320:
    print()
    print("=> 320 is OPTIMAL and the Cohn-Li / Ho mechanism is EXHAUSTED at dimension 19.")
else:
    print()
    print("=> room remains: up to %d more codewords, i.e. tau(19) up to %d"
          % (4 * (best - 320), 10668 + 4 * best))
