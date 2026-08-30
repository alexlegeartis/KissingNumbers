"""Exact ceilings for the Cohn-Li mechanism in dimensions 18, 19, 21.

Their extra vectors need a code of minimum distance >= n/4 inside the dual of the
constant-weight code, and tau(n) = base(n) + (code size) with base from their Table 2.1:
    n     18     19     20     21
    base  6500   10668  17400  27720
    record 7654  11948  19448  29768
So the code must beat 1154, 1280, 2048, 2048 respectively.

The dual is a punctured extended Golay code, so this is again a maximum independent set in a
Cayley graph on F_2^12 with S = the nonzero words of punctured weight < n/4.  For an abelian
Cayley graph the Delsarte LP equals the Lovasz theta, and alpha <= theta."""
import sys,os,numpy as np
from scipy.optimize import linprog
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'common'))
from capcheck import golay
import math
GOL=np.array(sorted(golay()),dtype=np.uint8)
BASE={18:6500,19:10668,20:17400,21:27720}
REC={18:7654,19:11948,20:19448,21:29768}
def run(n):
    p=24-n; d=math.ceil(n/4); thr=d-1
    Q=GOL[:,p:]; W=Q.sum(1)
    codes=[int("".join(map(str,r)),2) for r in Q]
    assert len(set(codes))==4096
    S=sorted(codes[i] for i in range(len(codes)) if 0<W[i]<=thr)
    span={0}
    for s in S: span|={x^s for x in span}
    comp=sorted(span); m=len(comp); ncomp=4096//m
    basis,sp=[],{0}
    for c in comp:
        if c and c not in sp:
            basis.append(c); sp|={x^c for x in sp}
    r=len(basis); lab={}
    for k in range(1<<r):
        v=0
        for b in range(r):
            if k>>b&1: v^=basis[b]
        lab[v]=k
    Ss=set(lab[s] for s in S)
    POP=np.array([bin(i).count("1") for i in range(1<<r)],dtype=np.int8)
    idx=np.arange(1<<r)
    H=(1-2*(POP[idx[:,None]&idx[None,:]]&1)).astype(float)
    N=1<<r
    bnds=[(1.,1.) if g==0 else ((0.,0.) if g in Ss else (0.,None)) for g in range(N)]
    res=linprog(c=-np.ones(N),A_ub=-H,b_ub=np.zeros(N),bounds=bnds,method='highs')
    th=-res.fun
    need=REC[n]-BASE[n]+1
    tot=ncomp*th
    print("  n=%2d  d>=%d  |S|=%3d  rank %2d  %d component(s) of %4d"%(n,d,len(S),r,ncomp,m))
    print("        theta(component) = %.4f   ceiling on the code = %.1f"%(th,tot))
    print("        need >= %d to beat the record %d   ->  %s"%(need,REC[n],
          "OPEN, room = %.0f"%(tot-need+1) if tot>=need else "DEAD"))
    return tot
print(__doc__); print()
for n in (21,19,18):
    run(n); print()
