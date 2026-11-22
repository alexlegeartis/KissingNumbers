"""Alternate arithmetic for the local 60D augmentation. Run with verify.py arguments."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product, combinations
import runpy
import numpy as np
d=runpy.run_path(str(Path(__file__).with_name("verify.py")))
G,W,Z,I,D,H,S=(d[k] for k in ('G','W','Z','I','D','H','S'))
A=W@G
def two_radicals(a,b):
    gap=F(1,4)-a-b
    return gap>=0 and gap*gap>=4*a*b

new_W=np.array([s*W[j] for j in I for s in (-1,1)],dtype='i8')
new_Z=np.repeat(Z,2,axis=0)
Gram=new_W@G@new_W.T+2*(new_Z@new_Z.T)
assert np.all(np.diag(Gram)==32)
assert np.max(Gram[np.triu_indices(len(Gram),1)])<=16
assert np.max(np.abs(new_Z@D.T))<=4
assert not np.any(new_W@G@H.T)
assert np.max(np.sum(abs(new_Z),axis=1))<=8

sparse=0
for j in (0,6):
    for q in range(8):
        supp=[j+2*k+((q>>k)&1) for k in range(3)]
        modes=[]
        for s in product((-1,1),repeat=3):
            p=np.zeros(12,dtype='i8');p[supp]=s;modes.append(p)
        for z,i in zip(Z,I,strict=True):
            cap=max(abs(int(x)) for x in z@np.array(modes).T)
            for v in S[min(q,q^7)]:
                h=abs(int(A[i]@v))
                assert two_radicals(F(h*h,512),F(5*cap*cap,384))
                sparse+=1
heavy=0
for row in d['root']['upgrade']:
    left=set(d['balls'][row['color']])
    for heads in row['tetrahedra']:
        c=min(left);words={c^t for t in d['L']};left-=words
        if heads is None: continue
        X=np.array([[1-2*(w>>j&1) for j in range(12)]for w in sorted(words)],dtype='i8')
        for z,i in zip(Z,I,strict=True):
            cap=max(abs(int(x)) for x in z@X.T)
            for v in np.array(heads,dtype='i8'):
                h=abs(int(A[i]@v))
                assert two_radicals(F(h*h,480),F(cap*cap,320))
                heavy+=1
    assert not left

# Replay the blocker union with actual integer matrix multiplication.
# All x and transform images are overinclusive as in verify.py, and
# every bound below is checked before the products.
T,V,E,O,M,Y,d0,perm=(d[k] for k in ('T','V','E','O','M','Y','d0','perm'))
cost=seen=0
for j,t in enumerate(T):
    cut=np.array([max(abs(int(q)) for q in block@D[perm[j]]) for block in [Z[I==i]for i in range(len(W))]])
    x=np.vstack((V[M[j]], E[O[j]:O[j+1]], np.array(d0.get(j,[]),dtype='i8').reshape(-1,48),Y[j]))
    # Columns where |tail product|<=2 cannot block, since mother minimum gives |w.x|<=6.
    jx=np.flatnonzero(cut>2)
    if len(jx):
        y=t@(G@W[jx].T)
        assert 48*max(abs(int(x.min())),abs(int(x.max())))*max(abs(int(y.min())),abs(int(y.max())))<2**62
        bad=np.any(np.abs(x.astype('i8')@y.astype('i8'))>8-cut[jx],axis=1)
        assert int(bad.sum())==d['class_costs'][j]
        cost+=int(bad.sum())
    else:
        assert d['class_costs'][j]==0
    seen+=len(x)
assert cost==d['cost']
print(f'PASS independent literal checks: {sparse} sparse, {heavy} heavy, {seen} superset owner rows, {cost} charged blockers.')
