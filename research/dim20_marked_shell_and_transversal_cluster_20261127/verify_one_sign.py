#!/usr/bin/env python3
"""Exact Golay/tetrad one-coordinate cut; Python 3 + NumPy, no solver."""
from itertools import combinations,product
from collections import Counter,defaultdict
import numpy as np

n=20; marking=0; removed=(0,1,2,3)
G=np.zeros((12,24),np.int8)
for i in range(12):
    for x in (0,2,4,5,6,10,11): G[i,(i+x)%23]=1
G[:,23]=G[:,:23].sum(1)%2
C=((np.arange(4096)[:,None]>>np.arange(12)&1)@G)%2
assert len(np.unique(C,axis=0))==4096
assert Counter(C.sum(1))=={0:1,8:759,12:2576,16:759,24:1}
O=C[(C[:,:4]==0).all(1)][:,4:]; O=O[O.sum(1)==8]
tetrads=[np.flatnonzero(u[4:]) for u in C[(C[:,:4]==1).all(1)&(C.sum(1)==8)]]
assert len(tetrads)==5 and len(set(np.concatenate(tetrads)))==20
P=C[:,4:]; assert len(np.unique(P,axis=0))==4096
s=np.array([t[0] for t in tetrads])
F=P[P[:,s].sum(1)%2==0]
assert len(F)==2048 and (F.sum(0)==1024).all()
assert len(O)==130 and (O.sum(0)==52).all()
assert np.all((F@O.T)%2==0)
for marking in range(n):
    u=np.array(list(combinations([i for i in range(n) if i!=marking],7)),np.int8)
    b=np.zeros(len(u),bool)
    for o in O[O[:,marking]==0]: b|=o[u].sum(1)>=5
    for o in O[O[:,marking]==1]: b|=o[u].sum(1)==6
    M=u[~b]; assert len(M)==260
    mask=[sum(1<<int(x) for x in a) for a in M]
    fs=defaultdict(list)
    for j,a in enumerate(M):
        for t in combinations(a.tolist(),4): fs[t].append(j)
    groups=[]
    for t, ix in fs.items():
        for v in combinations(ix,5):
            if all((mask[i]&mask[j]).bit_count()==4 for i,j in combinations(v,2)):
                groups.append(v)
    deg=Counter(x for v in groups for x in v)
    assert len(groups)==256 and Counter(deg.values())=={4:240,16:20}
    sp={x for x in deg if deg[x]==16}
    assert Counter(len(sp.intersection(v)) for v in groups)=={1:240,5:16}
    chosen=[v for v in groups if len(sp.intersection(v))==1]
    w=Counter(x for v in chosen for x in v)
    assert len(w)==260 and min(w.values())==4
print('PASS: every marked coordinate; 260 support shapes; exact quarter-weight clique cover bounds all signed refills by 3,840, and full trade by 18,936.')
