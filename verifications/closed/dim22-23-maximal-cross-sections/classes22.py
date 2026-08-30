"""(c) Leech minimal vectors classified by (<u1,v>,<u2,v>) for a 60-degree pair,
projected into the 22-dim orthogonal complement.  Exact rational admissibility."""
import numpy as np, collections, itertools
from fractions import Fraction as F
from base import leech

A = leech()
u1 = A[0].copy()
u2 = A[np.nonzero(A@u1==16)[0][0]].copy()
assert u1@u1==32 and u2@u2==32 and u1@u2==16
c1 = A@u1; c2 = A@u2
cls = collections.Counter(zip((c1//8).tolist(), (c2//8).tolist()))
print("class sizes (a1,a2 = inner products / 8):")
tot=0
Qof={}
for (a,b),n in sorted(cls.items(), key=lambda kv:(-kv[1],kv[0])):
    Q = a*a+b*b-a*b
    Qof[(a,b)]=Q
    tot+=n
    print("   a=(%2d,%2d)  Q=%2d  |v'|^2=%s   n=%6d"%(a,b,Q,F(32)-F(8*Q,3),n))
print("total", tot)

# --- exact pair admissibility between classes -------------------------------
def Bform(a,bb):
    return F(8,3)*(a[0]*bb[0]+a[1]*bb[1]) - F(4,3)*(a[0]*bb[1]+a[1]*bb[0])
def n2(a): return F(32) - F(8,3)*(a[0]**2+a[1]**2-a[0]*a[1])

CL = [k for k in cls if n2(k) > 0]
CL.sort(key=lambda k:(Qof[k],k))
print("\nout-of-plane classes:", len(CL), " (in-plane:", [k for k in cls if n2(k)==0], ")")

def badP(a,bb, same):
    """inner products P (of the 24-dim vectors) that VIOLATE cos<=1/2"""
    bad=[]
    r2 = n2(a)*n2(bb)
    for P in (-32,-16,-8,0,8,16,32):
        if same and P==32: continue        # v=w
        if (not same) and P==32: continue  # impossible across classes
        x = F(P) - Bform(a,bb)
        if x<=0: continue
        if 4*x*x > r2: bad.append(P)
    return bad

rows={}
for a in CL:
    for bb in CL:
        rows[(a,bb)] = badP(a,bb, a==bb)
free = [ (a,bb) for a in CL for bb in CL if not rows[(a,bb)] ]
print("\nconflict-free ORDERED class pairs: %d of %d"%(len(free), len(CL)**2))

# build the "compatibility graph" on classes (a class is self-compatible or not)
selfok = {a: not rows[(a,a)] for a in CL}
print("classes internally conflict-free:", sum(selfok.values()), "of", len(CL))
import sys
# maximum-weight clique in compatibility graph over self-ok classes
good = [a for a in CL if selfok[a]]
adj = {a:set(b for b in good if b!=a and not rows[(a,b)] and not rows[(b,a)]) for a in good}
best=[0,None]
def expand(R, P_, w):
    if w + sum(cls[x] for x in P_) <= best[0]: return
    if not P_:
        if w>best[0]: best[0]=w; best[1]=list(R)
        return
    for i,v in enumerate(list(P_)):
        expand(R+[v], [x for x in P_[i+1:] if x in adj[v]], w+cls[v])
expand([], sorted(good, key=lambda a:-cls[a]), 0)
print("\nMAX conflict-free union of classes: %d   classes %s"%(best[0], best[1]))
print("Lambda_22 (class (0,0)) = %d"%cls[(0,0)])
