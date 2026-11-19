"""Independent exact short-shell and block-distance check for the two witnesses."""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
from collections import Counter
import hashlib,json
import numpy as np
P=Path(__file__).resolve().parent
c=json.loads((P/'certificate.json').read_text())
Gn=c['gram_12']; Mn=c['sublattice_matrix']
assert len(Gn)==len(Mn)==12 and all(len(r)==12 for r in Gn+Mn)
assert all(type(x) is int and abs(x)<1024 for r in Gn+Mn for x in r)
G=np.asarray(Gn,dtype=np.int64); M=np.asarray(Mn,dtype=np.int64)
assert np.array_equal(G,G.T) and np.array_equal(M@G@M.T,3*G)
# This identity and positive definite G imply |det M|=3^6, with no numerical determinant.
u=c['label_first'];v=c['label_second']
assert sorted(u+v)==list(range(12)) and len(u)==len(v)==6
assert np.all((M[:,u]+M[:,v])%3==0)
# Positive pivots and the exact completed-square identities from symmetric
# elimination. The remaining suffix quadratic stays nonnegative.
A=[[F(x) for x in row] for row in Gn]; stages=[]
while A:
    a=A[0][0]; assert a>0
    h=[x/a for x in A[0][1:]]
    stages.append((a,h))
    A=[[A[i+1][j+1]-a*x*y for j,y in enumerate(h)] for i,x in enumerate(h)]
X=[0]*12
hist=Counter(); w=[500]*729; mult=[0]*729
w[0]=0;mult[0]=1
pow3=[3**i for i in range(6)]
def search(j,R):
    if j<0:
        n=sum(X[a]*Gn[a][b]*X[b] for a in range(12) for b in range(12))
        assert n==8-R
        if n!=0:
            hist[n]+=1
            label=sum(((X[a]+X[b])%3)*p for a,b,p in zip(u,v,pow3))
            if n<w[label]:w[label]=n;mult[label]=0
            if n==w[label]:mult[label]+=1
        return
    a,h=stages[j]
    S=sum((b*x for b,x in zip(h,X[j+1:])),F(0))
    # If denominator(S)=d, then the integer d*X_j+numerator(S)
    # must have absolute value at most floor(sqrt(R*d*d/a)).
    den=S.denominator
    val=R*den**2/a
    rad=isqrt(val.numerator//val.denominator)
    lo=-((rad+S.numerator)//den)
    hi=(rad-S.numerator)//den
    for x in range(lo,hi+1):
        X[j]=x
        rem=R-a*(x+S)**2
        assert rem>=0
        search(j-1,rem)
search(11,F(8))
assert hist=={4:756,6:4032,8:20412}
assert Counter(zip(w,mult))=={(0,1):1,(4,3):252,(6,18):224,(8,81):252}
def array(name):
    V=c[name]
    assert all(isinstance(x,list) and len(x)==3 and all(type(y) is int and 0<=y<729 for y in x) for x in V)
    assert len(set(map(tuple,V)))==len(V)
    return np.asarray(V,dtype=np.int64)
E=array('outer_cosets')
assert np.array_equal(E[0],[0,0,0])
w=np.asarray(w,dtype=np.int64);mult=np.asarray(mult,dtype=np.int64)
assert np.all(w[E[1:]].sum(1)==12)
N36=3*hist[4]+np.prod(mult[E[1:]],axis=1).sum()
assert N36==c['points']
digs=(np.arange(729)[:,None]//np.asarray(pow3))%3
W=w[((digs[:,None,:]-digs[None,:,:])%3)@np.asarray(pow3)]
def min_dist(B,C,same=False):
    m=100
    for j in range(0,len(B),96):
        v=B[j:j+96]
        a=sum(W[v[:,None,i],C[None,:,i]] for i in range(3))
        if same:a[np.arange(len(v)),np.arange(j,j+len(v))]=100
        m=min(m,int(a.min()))
    return m
minimum=min_dist(E,E,True)
assert minimum>=12
print(json.dumps({'status':'PASS','points_36':int(N36),'distance':minimum,
                  'certificate_sha256':hashlib.sha256((P/'certificate.json').read_bytes()).hexdigest()}))
