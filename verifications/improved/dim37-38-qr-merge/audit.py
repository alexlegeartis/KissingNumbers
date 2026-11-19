"""Independent checks for the projected QR certificates (NumPy)."""
import argparse
import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from math import isqrt
import numpy as np
from supports37 import supports

P=Path(__file__).resolve().parent
PALETTES={37:(4,1,-4,-1),38:(3,None,0,-3)}

def load_cube(d):
    code=json.loads((P/f'dense_parent{d+1}.json').read_text())
    def row(k):return np.array([int(x,16) for x in code[k]],dtype=np.uint64)
    G=row('generators_hex')
    R=row('representatives_hex')
    n=1<<len(G)
    words=np.zeros(n,dtype=np.uint64)
    for i,g in enumerate(G):
        words[1<<i:1<<(i+1)]=words[:1<<i]^g
    return code,words,R

def verify(d,full):
    code,H,R=load_cube(d)
    f=PALETTES[d];ordinary=d-1
    sizes=[0 if x is None else 4*ordinary+x*x for x in f]
    assert d+1==code['length']
    # Independent, exact reduced computation. For each pair of affine
    # offsets, enumerate every full xor mask from the kernel. Its high
    # part gives the projected body product; its low bits give the xor
    # of the merged patterns. The calculation does not depend on a bound
    # imported from the original length-38/39 code paper.
    floor=Counter(); extrema={}
    for i,a in enumerate(R):
        for j,b in enumerate(R[:i+1]):
            D=H^a^b
            if i==j:
                assert D[0]==0 and (D[1:]!=0).all()
                D=D[1:]
            else:assert (D!=0).all()
            words=D.tolist()
            floor[min(x.bit_count() for x in words)]+=1
            mx=[-1]*4
            for x in words:
                h=(x>>2).bit_count()
                mx[x&3]=h if mx[x&3]<0 else min(h,mx[x&3])
            for a0 in range(4):
                if f[a0] is None:continue
                for b0 in range(4):
                    if f[b0] is None or mx[a0^b0]<0:continue
                    v=4*(ordinary-2*mx[a0^b0])+f[a0]*f[b0]
                    assert v<=0 or 4*v*v<=sizes[a0]*sizes[b0]
                    k=tuple(sorted((sizes[a0],sizes[b0])))
                    extrema[k]=max(extrema.get(k,0),v)
    assert floor==({10:55,12:11} if d==37 else {10:45,12:10})
    Q=np.ravel(H[None,:]^R[:,None]);Q=Q[np.isin(Q&3,[i for i in range(4) if f[i] is not None])]
    assert len(np.unique(Q))==len(Q)
    S=supports() if d==37 else [int(x[::-1],2) for x in (P/'supports_3088.txt').read_text().split()]
    B=np.array(S,dtype=np.uint64)
    # Explicit middle sparsity checking through unique five-subsets.
    seen=set()
    for x in S:
        ids=[j for j in range(d) if x>>j&1]
        assert len(ids)==8
        for k in combinations(ids,5):
            assert k not in seen
            seen.add(k)
    if d==37:assert all(x<1<<(d-1) for x in S)
    for i,s in enumerate(f):
        if s is None:continue
        dot=32 if d==37 else 2*(14+max(2,abs(s)))
        assert 4*dot*dot<=32*sizes[i]
        r=4*(2+max(2,abs(s)))
        assert 4*r*r<=32*sizes[i]
    r=d*(d-1)*2
    n=len(Q)+128*len(B)+r
    assert n=={37:533864,38:643836}[d]
    result={'dimension':d,'points':n,'dense':len(Q),'sparse_supports':len(B),'roots':r,
      'dense_extrema':[list(k)+[v] for k,v in sorted(extrema.items())],
      'checked_five_subsets':len(seen),'status':'PASS'}
    if full:
        t=(1-2*((Q[:,None]>>np.arange(2,d+1,dtype=np.uint64))&1).astype(np.int8))*2
        tail=np.array([0 if x is None else x for x in f],dtype=np.int8)[Q&3]
        X=np.c_[t,tail].astype(np.float32)
        assert d*int(np.max(np.abs(X)))**2<2**24
        N=(X*X).sum(axis=1).astype(np.int64)
        assert set(N.tolist())=={n for n in sizes if n}
        kinds=np.unique(N)
        pairs=0
        for i in range(0,len(X),512):
            L=X[i:i+512]
            V=X[i:]
            b=L@V.T
            np.fill_diagonal(b[:,:len(L)],0)
            for a in kinds:
                I=N[i:i+len(L)]==a
                if not I.any():continue
                check=np.full(len(V),-1,dtype=np.float32)
                for v in kinds:
                    check[N[i:]==v]=isqrt(int(a*v)//4)
                assert not np.any(b[I]>check[None,:])
            pairs+=len(L)*(len(V)-len(L))+len(L)*(len(L)-1)//2
            if i and i%65536==0:print(d,'dense rows',i,flush=True)
        assert pairs==len(X)*(len(X)-1)//2
        result['all_dense_pairs']=pairs
    print(json.dumps(result,sort_keys=True),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dimension',type=int,choices=[37,38])
    p.add_argument('--full-dense',action='store_true')
    a=p.parse_args()
    for d in ([a.dimension] if a.dimension else [37,38]):
        verify(d,a.full_dense)
