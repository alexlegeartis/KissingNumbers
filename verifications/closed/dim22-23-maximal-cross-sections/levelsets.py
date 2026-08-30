"""Idea (b): unions of Leech level sets projected into u^perp.
Norm^2=32 scaling.  v' = v - (C/32) u, |v'|^2 = 32 - C^2/32.
Pair (v,w) admissible iff  P - Cv*Cw/32 <= (1/2)sqrt((32-Cv^2/32)(32-Cw^2/32)).
"""
import numpy as np, collections, itertools, math
from base import leech
from fractions import Fraction as F

A = leech(); u = A[0].copy(); C = A @ u
LEV = {0:np.nonzero(C==0)[0], 8:np.nonzero(C==8)[0], 16:np.nonzero(C==16)[0],
       -8:np.nonzero(C==-8)[0], -16:np.nonzero(C==-16)[0]}
print("sizes:", {k:len(v) for k,v in LEV.items()})

def norm2(c): return F(32) - F(c*c,32)
def forbidden(a,b):
    """exact set of P values that violate, over P in {-32,-16,-8,0,8,16} (P=32 => same vec)"""
    bad=[]
    lhsc = F(a*b,32)
    r2 = norm2(a)*norm2(b)          # rhs = sqrt(r2)/2 ; compare (P-lhsc) with that
    for P in (-32,-16,-8,0,8,16,32):
        x = F(P) - lhsc
        if x <= 0: continue
        if 4*x*x > r2: bad.append(P)
    return bad

print("\nforbidden inner products by level pair (exact rational test):")
for a,b in itertools.combinations_with_replacement([0,8,16,-8,-16],2):
    f = forbidden(a,b)
    if f: print("   (%4d,%4d): P in %s"%(a,b,f))
print("   all other level pairs: no constraint")

# --- degrees of the conflict graph -------------------------------------------
def deg(src_lev, tgt_lev, P, nsample=None):
    S = A[LEV[src_lev]]; T = A[LEV[tgt_lev]]
    idx = range(len(S)) if nsample is None else np.random.default_rng(0).choice(len(S), nsample, replace=False)
    ds = []
    for i in idx:
        ds.append(int((T @ S[i] == P).sum()))
    return collections.Counter(ds)

print("\nconflict-graph degrees (P=16 edges):")
for (s,t) in [(0,8),(8,0),(0,16),(16,0)]:
    c = deg(s,t,16,nsample=40)
    print("   level %3d -> level %3d : %s"%(s,t,dict(c)))
