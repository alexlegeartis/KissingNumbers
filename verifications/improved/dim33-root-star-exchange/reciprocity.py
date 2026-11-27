"""Exact finite inputs for the star/roots reciprocal-charge upper bound."""
import itertools
from collections import Counter
from pathlib import Path
import numpy as np

W=np.load(Path(__file__).resolve().parent/'witness.npz',allow_pickle=False)
S=W['supports'].astype(np.uint64)
sparse=W['sparse'].astype(np.int32)
octs=sparse[np.count_nonzero(sparse,axis=1)==8]
assert len(octs)==2216
Q=np.array(list(itertools.product((-1,1),repeat=5)),np.int32)*[4,2,2,2,2]
U=[]
for h in range(33):
    S_h=S[((S>>np.uint64(h))&1)!=0]
    leaf=np.array(list(itertools.combinations([i for i in range(33) if i!=h],4)),np.uint64)
    masks=np.bitwise_or.reduce(np.uint64(1)<<leaf,axis=1)
    pass0=[]
    for i in range(0,len(leaf),2048):
        sel=(np.bitwise_count(masks[i:i+2048,None]&S_h[None,:])<=2).all(axis=1)
        pass0.extend(leaf[i:i+2048][sel].tolist())
    for tup in pass0:
        a=np.zeros((32,33),np.int32)
        a[:,[h,*tup]]=Q
        U.extend(a[(a@octs.T<=16).all(axis=1)])
A=np.asarray(U,np.int32)
assert A.shape==(1486,33) and len(np.unique(A,axis=0))==1486

def head(r):
    h=int(np.argmax(abs(r)))
    return h, int(r[h])//4

def leaves(r):
    return {(int(i),int(r[i])//2) for i in np.flatnonzero(abs(r)==2)}

uses={(head(x),v) for x in A for v in leaves(x)}
both={tuple(sorted((u,v))) for u,v in uses if (v,u) in uses}
exc={(20,-1),(20,1),(30,-1),(30,1)}
base={r for r in both if all(x not in exc for x in r)}
assert len(both)==34 and len(base)==28

local={}
for u in sorted(exc):
    X=[x for x in A if head(x)==u]
    deficit=[]
    for k in range(1,len(X)+1):
        for rows in itertools.combinations(X,k):
            g=np.asarray(rows)@np.asarray(rows).T
            np.fill_diagonal(g,-32)
            if g.max()<=16:
                deficit.append(len(set.union(*(leaves(x) for x in rows)))-k)
    d=min(deficit)
    q=sum(u in r for r in both)
    assert d==3 and q<=2
    local[u]=(len(X),q,d)
print('eligible stars',len(A),'possible reciprocal roots',len(both))
print('outside exceptional heads',len(base),'; exceptional (candidate count, possible credit, minimum deficit):',local)
print('upper bound with fixed octads: 131072 + 1837*128 + 2216 + 2112 + 28 =',131072+1837*128+2216+2112+28)
