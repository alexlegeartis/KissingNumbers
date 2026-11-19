"""Independent integer discriminant replay for the 36D rebuilt cell layer."""
from collections import Counter
from fractions import Fraction as F
from math import isqrt,lcm
from pathlib import Path
import numpy as np,json
P=Path(__file__).resolve().parent
a=json.loads((P/'certificate.json').read_text())
G=np.asarray(a['gram_12'],dtype=np.int64)
V=np.asarray(a['sublattice_matrix'],dtype=np.int64)
I,J=np.array(a['label_first']),np.array(a['label_second'])
assert G.shape==V.shape==(12,12) and np.array_equal(G,G.T) and (G.diagonal()%2==0).all()
assert np.array_equal(V@G@V.T,3*G)
assert np.array_equal(np.sort(np.r_[I,J]),np.arange(12))
assert ((V[:,I]+V[:,J])%3==0).all()
assert ((G[I]-G[J])%3==0).all()
# Complete squares with the opposite variable ordering to the source
# rational LDL routine; enumerate solely via bounded integer squares.
R=[[F(int(x)) for x in row] for row in G[::-1,::-1]]
eq=[]
for i in range(12):
    d=R[i][i]
    assert d>0
    row=[F(1)]+[R[i][j]/d for j in range(i+1,12)]
    den=lcm(*(x.denominator for x in row))
    eq.append((d/F(den**2),den,[int(den*x) for x in row[1:]]))
    for j in range(i+1,12):
        for k in range(i+1,12):
            R[j][k]-=R[i][j]*R[i][k]/d
L=lcm(*(w.denominator for w,_,_ in eq))
terms=[(int(L*w),den,z) for w,den,z in eq]
co=np.array(a['outer_cosets'],dtype=np.int64)
assert co.ndim==2 and co.shape[1]==3 and len(np.unique(co,axis=0))==len(co)
assert (0<=co).all() and (co<729).all() and (co[0]==0).all()
base=3**np.arange(6,dtype=np.int64)
u=[0]*12
s=[Counter() for _ in range(729)]
def rec(k,left):
    if k<0:
        assert (8*L-left)%L==0
        v=u[::-1]
        lab=sum(((v[i]+v[j])%3)*p for i,j,p in zip(I,J,base.tolist()))
        s[lab][(8*L-left)//L]+=1
        return
    w,den,z=terms[k]; rad=isqrt(left//w); shift=sum(t*v for t,v in zip(z,u[k+1:]))
    lo=-((rad+shift)//den);hi=(rad-shift)//den
    for x in range(lo,hi+1):
        u[k]=x
        rec(k-1,left-w*(den*x+shift)**2)
rec(11,8*L)
hist=sum(s,Counter())
assert hist=={0:1,4:756,6:4032,8:20412}
X=(np.arange(729)[:,None]//base)%3
Q=2*G[np.ix_(I,I)]%3
q=np.sum((X@Q)*X,axis=1)%3
M=np.array([6,8,4],dtype=np.int64)[q];M[0]=0
assert all(min(shell)==m for shell,m in zip(s,M))
weights=np.array([shell[m] for shell,m in zip(s,M)])
assert Counter(zip(M,weights))=={(0,1):1,(4,3):252,(6,18):224,(8,81):252}
assert np.all(M[co[1:]].sum(axis=1)==12)
tab=M[(((X[:,None,:]-X[None,:,:])%3)*base).sum(axis=2)]
# Independent range of the label map: its kernel has index 729. Similarity
# gives index 3^6=729 (determinant via determinants of the nonsingular
# Grams), so the sublattice checked above is exactly the label kernel.
B=len(co);pairs=0;mi=99
for i,(x,y,z) in enumerate(co):
    if i:
        a0,b0,c0=co[:i].T
        v=tab[x,a0]+tab[y,b0]+tab[z,c0]
        mi=min(mi,int(v.min()));pairs+=len(v)
assert pairs==B*(B-1)//2 and mi==12
total=3*hist[4]+int(weights[co[1:]].prod(axis=1).sum())
assert total==a['points']
lift=np.array(a['lift_cosets'],dtype=np.int64)
assert a['lift_heights']==[-2,2]
assert lift.ndim==2 and lift.shape[1]==3 and len(np.unique(lift,axis=0))==len(lift)
assert np.all((0<=lift)&(lift<729))
assert np.array_equal(np.sort(M[lift],axis=1),np.tile([0,4,4],(len(lift),1)))
L0,L1,L2=lift.T
D0,D1,D2=co.T
LL=tab[L0[:,None],L0]+tab[L1[:,None],L1]+tab[L2[:,None],L2]
LE=tab[L0[:,None],D0]+tab[L1[:,None],D1]+tab[L2[:,None],D2]
p0=int(LL[np.triu_indices(len(lift),1)].min());p1=int(LE.min())
assert (p0,p1)==(12,8)
total37=total+2*int(weights[lift].prod(axis=1).sum())
assert total37==a['points_37']
print(json.dumps({'status':'PASS','points_36':total,'points_37':total37,'outer_pairs_checked':pairs,
                  'lift_lift':len(lift)*(len(lift)-1)//2,'lift_equator':len(LE)*len(co),
                  'minimum_distances':[mi,p0,p1],'shell_counts':dict(hist)},sort_keys=True))
