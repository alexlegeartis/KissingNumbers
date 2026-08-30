#!/usr/bin/env python3
"""The cap-direction line systems for dimensions 57-63, extracted and verified exactly.

At cap level t = 3/4 a part is an antipodal pair and distinct parts are at <z,z'> <= 1/2, so
the direction set W is a union of antipodal pairs whose LINES are pairwise at >= 60 degrees.
Writing

    lambda(k) = |W| / 2 = the number of parts = the number of classes consumed,

the gain is 2 * (sum of the lambda(k) largest class sizes), strictly increasing in lambda(k).

The package used lambda(k) = tau_lat(k)/2, the lines of the best k-dimensional LATTICE.  That
is not the binding quantity: ANY antipodally symmetric 60-degree code of size 2m gives m lines.
Cohn's record configurations are antipodally symmetric for k = 9, 10, 11, 13, 14, 15 -- checked
here in exact arithmetic -- so lambda(k) = tau(k)/2 is available there instead.

    k         9    10    11    12    13    14    15
    tau_lat/2 136  168   219   378   459   711   1170
    tau/2     153  255   302   (420) 577   966   1282
    used      153  255   302   378   577   966   1282

(k = 12 is excluded: tau(12) = 841 is ODD, so that configuration cannot be antipodally
symmetric, and the best antipodal system known in R^12 remains K_12's 378 lines.)

Stores each configuration as integers (A + B sqrt(R))/den with the inner-product coefficients."""
import numpy as np, re, os, sys
from fractions import Fraction as F
from math import gcd
HERE=os.path.dirname(os.path.abspath(__file__))
def parse(tok):
    t=tok.replace(' ','')
    m=re.fullmatch(r'(-?)(\d*)\*?sqrt\((\d+)\)(?:/(\d+))?',t)
    if m:
        s=-1 if m.group(1)=='-' else 1
        num=int(m.group(2)) if m.group(2) else 1
        return (F(0),F(s*num,int(m.group(4)) if m.group(4) else 1),int(m.group(3)))
    return (F(t),F(0),None)
def load(dim,path):
    pts=[];coef=None;want=False
    for line in open(path):
        line=line.strip()
        if line.startswith('Dimension:'):
            if want and pts: break
            want=(int(line.split(':')[1])==dim); pts=[] if want else pts
        elif line.startswith('Inner product'):
            if want: coef=[int(x) for x in line.split(':')[1].split(',')]
        elif line.startswith('Number') or line=='Points:' or not line: continue
        elif want: pts.append([parse(x) for x in line.split(',')])
    return pts,np.array(coef,np.int64)
SRC=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'..','..','..','..','..','dimensions1-24.txt')
print(__doc__)
print("   k   |W|   sqrt   den   antipodal   60-degree   lambda(k)")
for k in (9,10,11,13,14,15):
    pts,c=load(k,SRC); n=len(pts)
    rads={p[2] for row in pts for p in row if p[2] is not None}
    assert len(rads)<=1
    R=rads.pop() if rads else 0
    den=1
    for row in pts:
        for a,b,_ in row:
            den=den*a.denominator//gcd(den,a.denominator)
            den=den*b.denominator//gcd(den,b.denominator)
    A=np.array([[int(p[0]*den) for p in row] for row in pts],np.int64)
    B=np.array([[int(p[1]*den) for p in row] for row in pts],np.int64)
    P=(A*c)@A.T + R*((B*c)@B.T)
    Q=(A*c)@B.T + (B*c)@A.T
    N=int(P[0,0])
    assert (np.diag(P)==N).all() and (np.diag(Q)==0).all()
    key={(A[i].tobytes(),B[i].tobytes()) for i in range(n)}
    sym=all(((-A[i]).tobytes(),(-B[i]).tobytes()) in key for i in range(n))
    def le(L,Rr):
        if R==0: return L<=0
        if Rr>=0 and L<=0: return True
        if Rr<0 and L>0: return False
        return (L*L<=R*Rr*Rr) if L>0 else (L*L>=R*Rr*Rr)
    bad=0
    for i in range(n):
        for j in np.nonzero(2*P[i]-N>0)[0]:
            if j!=i and not le(2*int(P[i,j])-N,-2*int(Q[i,j])): bad+=1
    assert sym and bad==0
    np.savez_compressed(os.path.join(HERE,'..','data','capdirs_%d.npz'%k),A=A,B=B,c=c,
                        R=np.array([R]),den=np.array([den]),norm=np.array([N]))
    print("   %2d %5d  %5d %5d   %-9s   %-9s   %d"%(k,n,R,den,sym,"yes",n//2))
print()
print("saved capdirs_<k>.npz for k = 9, 10, 11, 13, 14, 15")
