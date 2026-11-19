"""Explicit dense–dense dot-product replay for 38D (NumPy required)."""
from itertools import combinations
from pathlib import Path
from collections import Counter
from math import isqrt
import argparse,json
import numpy as np
from verify import load,HEADS
P=Path(__file__).resolve().parent
S=P.parent/'dim37-38-qr-merge'

def replay(full_dense=False):
    K,Q=load()
    H=np.asarray(K,dtype=np.uint64)
    T=np.array([k for i,q in enumerate(Q) for k in [HEADS[i]]*len(q)],dtype=np.int64)
    A=np.array([v for g in Q for v in g],dtype=np.uint64)
    # Independent vectorized Hamming count on the actual 37-column values,
    # using 5 byte chunks with a byte-popcount lookup table.
    tab=np.array([i.bit_count() for i in range(256)],dtype=np.uint8)
    pair_bounds={}
    for i,a in enumerate(A):
        data=(A[:i,None]^a^H[None,:])
        count=np.zeros(data.shape,dtype=np.uint8)
        for k in range(5):
            count+=tab[(data>>(8*k))&255]
        mi=count.min(axis=1)
        for j,m in enumerate(mi):
            dot=148-8*int(m)+int(T[i]*T[j])
            assert dot<=0 or 4*dot*dot<=int((148+T[i]**2)*(148+T[j]**2))
            key=tuple(sorted([int(T[i]),int(T[j])]))
            pair_bounds[key]=max(pair_bounds.get(key,0),dot)
    words=(A[:,None]^H).ravel()
    t=np.repeat(T,len(H))
    assert len(words)==253952 and all(np.unique(words[t==k]).size==(t==k).sum() for k in HEADS)
    h={(int(x),int(k)) for x,k in zip(words,t)}
    assert len(h)==253952
    # A second combinatorial proof of the signed octad conditions uses
    # disjoint five-subsets; two octads intersect >=5 iff they duplicate one.
    seen=set()
    supports=(S/'supports_3088.txt').read_text().splitlines()
    for row in supports:
        ix=[i for i,c in enumerate(row) if c=='1']
        assert len(row)==38 and len(ix)==8
        for k in combinations(ix,5):
            assert k not in seen
            seen.add(k)
    assert len(seen)==3088*56
    assert (4*34**2<=32*157 and 4*32**2<=32*148 and 4*20**2<=32*148)
    record=dict(dimension=38,dense=len(h),total=len(h)+128*len(supports)+2*38*37,status='PASS',
                product_extrema=[list(k)+[v] for k,v in sorted(pair_bounds.items())])
    if full_dense:
        Y=(1-2*((words[:,None]>>np.arange(37,dtype=np.uint64))&1).astype(np.int8))*2
        Z=np.c_[Y,t].astype(np.float32)
        # Inputs/products and all partial sums are small integers in exact
        # float32 range, independent of association of the BLAS additions.
        assert 38*int(max(abs(Z.max()),abs(Z.min())))**2<2**24
        N=(Z*Z).sum(1).astype(np.int64)
        total=0
        for i in range(0,len(Z),512):
            B=Z[i:i+512];C=Z[i:]
            P0=B@C.T
            np.fill_diagonal(P0[:,:len(B)],0)
            for n in (148,157):
                idx=N[i:i+len(B)]==n
                th=np.where(N[i:]==148,isqrt(n*148//4),isqrt(n*157//4))
                assert np.all(P0[idx]<=th)
            total+=len(B)*(len(C)-len(B))+len(B)*(len(B)-1)//2
            if i and i%65536==0: print('dense rows',i,flush=True)
        assert total==len(h)*(len(h)-1)//2
        record['all_dense_pairs']=total
    return record

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--full-dense',action='store_true')
    z=p.parse_args()
    print(json.dumps(replay(z.full_dense),sort_keys=True),flush=True)
