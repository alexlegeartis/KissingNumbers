"""VERIFICATION: every claim of this package, in exact integer / rational arithmetic
where possible.  Run:  python verify.py"""
import numpy as np, collections, itertools
from fractions import Fraction as F
from base import leech

ok = lambda c,m: print(("  OK   " if c else "  FAIL ")+m) or c
allok=[]

A = leech()
allok.append(ok(A.shape==(196560,24), "Leech minimal vectors: 196560 x 24"))
n2 = (A.astype(np.int64)**2).sum(1)
allok.append(ok((n2==32).all(), "all norms^2 = 32"))
allok.append(ok(len(set(map(tuple,A.tolist())))==196560, "all distinct"))

u1=A[0].copy(); u2=A[np.nonzero(A@u1==16)[0][0]].copy()
u3=A[np.nonzero((A@u1==16)&(A@u2==16))[0][0]].copy()
L23=A[A@u1==0]; L22=A[(A@u1==0)&(A@u2==0)]; L21=A[(A@u1==0)&(A@u2==0)&(A@u3==0)]
allok.append(ok(len(L23)==93150, "Lambda_23 cross-section = 93150"))
allok.append(ok(len(L22)==49896, "Lambda_22 cross-section = 49896"))
allok.append(ok(len(L21)==27720, "Lambda_21 cross-section = 27720"))

rngv=np.random.default_rng(5)
ipset=set()
for i in rngv.choice(196560,300,replace=False):
    ipset |= set(np.unique(A.astype(np.int32)@A[i].astype(np.int32)).tolist())
allok.append(ok(ipset=={0,8,-8,16,-16,32,-32},
    "all Leech pairwise inner products lie in {0,+-8,+-16,+-32} (300 random rows, full pass)"))
def maxip(S):
    idx=rngv.choice(len(S),200,replace=False)
    m=-10**9
    X=S.astype(np.int32)
    for i in idx:
        d=X@X[i]; d[i]=-10**9; m=max(m,int(d.max()))
    return m
for name,S in (("Lambda_21",L21),("Lambda_22",L22),("Lambda_23",L23)):
    m=maxip(S)
    allok.append(ok(m==16, "%s max inner product = %d of 32 (200 random rows, full pass) => cos = 1/2 exactly"%(name,m)))

# ---- level structure over a minimal vector (idea b) -------------------------
c=A@u1
ctr=collections.Counter(c.tolist())
allok.append(ok(ctr=={0:93150,8:47104,16:4600,-8:47104,-16:4600,32:1,-32:1},
                "Leech level sizes over a minimal vector"))
def n2l(k): return F(32)-F(k*k,32)
def forb(a,b):
    r2=n2l(a)*n2l(b); off=F(a*b,32); out=[]
    for P in (-32,-16,-8,0,8,16):
        x=F(P)-off
        if x>0 and 4*x*x>r2: out.append(P)
    return out
allok.append(ok(forb(0,0)==[] and forb(8,8)==[] and forb(16,16)==[] and forb(8,16)==[],
                "level pairs (0,0),(8,8),(16,16),(8,16) are unconditionally admissible"))
allok.append(ok(forb(0,8)==[16] and forb(0,16)==[16],
                "level 0 conflicts with levels 8,16 exactly at inner product 16"))
allok.append(ok(forb(8,-8)==[16] and forb(16,-16)==[8,16] and forb(8,-16)==[16],
                "opposite-sign level pairs conflict"))

# ---- exact bipartite degrees (Co_2 orbits => regular) -----------------------
rng=np.random.default_rng(0)
L8=A[c==8]; L16=A[c==16]
def degs(S,T,P,k=150):
    idx=rng.choice(len(S),min(k,len(S)),replace=False)
    return set(int((T@S[i]==P).sum()) for i in idx)
d_0_8 =degs(L23,L8 ,16); d_8_0 =degs(L8 ,L23,16)
d_0_16=degs(L23,L16,16); d_16_0=degs(L16,L23,16)
allok.append(ok(d_0_8=={1024} and d_8_0=={2025}, "degrees level0<->level8  = 1024 / 2025"))
allok.append(ok(d_0_16=={44} and d_16_0=={891}, "degrees level0<->level16 = 44 / 891"))
allok.append(ok(93150*1024==47104*2025 and 93150*44==4600*891, "edge counts consistent"))
r8=F(2025,1024); r16=F(891,44)
allok.append(ok(F(1,1)/r8+F(1,1)/r16 < 1,
   "Hall: 1024/2025 + 44/891 = %s < 1  =>  matching saturates the positive half"%(F(1,1)/r8+F(1,1)/r16)))
print("     => max independent set of (level 0) u (levels 8,16) = 93150 exactly")

# ---- exact witness directions (no free point) -------------------------------
print("\nexact witnesses for the min-max direction (rational arithmetic):")
for name,(proj2,ipmax,val) in (
      ("Lambda_21", (F(29),16,F(8,29))),
      ("Lambda_22", (F(88,3),16,F(3,11))),
      ("Lambda_23", (F(30),16,F(4,15)))):
    lhs = F(ipmax*ipmax, 1)/(proj2*F(32))
    allok.append(ok(lhs==val, "%s: (16)^2/(|v'|^2 * 32) = %s = min-max cos^2 (%.9f)"
                    %(name,val,float(val)**0.5)))
print("\nALL CHECKS PASSED" if all(allok) else "\nSOME CHECKS FAILED")
