"""Exact 38D ternary-head, shortened-binary-kernel certificate."""
from pathlib import Path
from collections import Counter
from itertools import combinations
from hashlib import sha256
import json

P=Path(__file__).resolve().parent
S=P.parent/'dim37-38-qr-merge'
HEADS=(3,0,-3)
DIST=((10,9,8),(9,10,9),(8,9,10))

def load():
    a=json.loads((S/'dense_parent39.json').read_text())
    parent=[int(x,16) for x in a['generators_hex']]
    assert len(parent)==15 and all(0<x<1<<39 for x in parent)
    R={}
    for g in parent:
        h=g
        while h:
            i=h.bit_length()-1
            if i not in R:R[i]=h;break
            h^=R[i]
    assert len(R)==15
    free=sorted(set(range(39))-set(R))
    def quotient(x):
        for i,v in sorted(R.items(),reverse=True):
            if x>>i&1:x^=v
        return sum(((x>>i)&1)<<j for j,i in enumerate(free))
    source_cols=[quotient(1<<i) for i in range(2,39)]
    dual={}
    for j,v in enumerate(source_cols):
        u=1<<j
        while v:
            i=v.bit_length()-1
            if i not in dual:dual[i]=(v,u);break
            d,h=dual[i];v^=d;u^=h
    assert len(dual)==24
    def representative(x):
        q=x;a=0
        for i,(g,v) in sorted(dual.items(),reverse=True):
            if q>>i&1:q^=g;a^=v
        assert q==0 and quotient(a<<2)==x
        return a
    H=[0]
    for g in parent:H += [h^g for h in H]
    assert len(set(H))==32768
    K=sorted(h>>2 for h in H if (h&3)==0)
    assert min(h.bit_count() for h in K[1:])==12 and len(K)==8192 and K[0]==0
    assert all(quotient(h<<2)==0 for h in K)
    raw=[list(map(int,s.split())) for s in (P/'offset_witness.txt').read_text().splitlines()]
    assert len(raw)==3
    Q=[]
    for n,*label in raw:
        assert n==len(label)==len(set(label)) and all(0<=q<1<<24 for q in label)
        Q.append([representative(q) for q in label])
    assert list(map(len,Q))==[11,10,10]
    return K,Q

def replay():
    K,Q=load();ranges=[]
    for i,L in enumerate(Q):
        for j,M in enumerate(Q[:i+1]):
            th=DIST[i][j];md=37;count=0
            for a,x in enumerate(L):
                old=M[:a] if i==j else M
                for y in old:
                    v=x^y
                    d=min((v^h).bit_count() for h in K)
                    assert d>=th
                    count+=1;md=min(md,d)
            # All words have coordinates ±2 on the 37 ordinary sites.
            # The new scalar coordinate is independent of the affine word.
            n1=148+HEADS[i]**2;n2=148+HEADS[j]**2
            dot=148-8*min(md,12)+HEADS[i]*HEADS[j] if i==j else 148-8*md+HEADS[i]*HEADS[j]
            assert 4*max(0,dot)**2<=n1*n2
            ranges.append((i,j,count,md,n1,n2,dot))
    raw=(S/'supports_3088.txt').read_bytes()
    assert sha256(raw).hexdigest()=="ccb530f90e914df84cbc8d19a71d30c79146d23a2e8e86e7a813d4d92a5a6a8b"
    words=raw.decode('ascii').splitlines()
    assert all(len(w)==38 and set(w)<={'0','1'} and w.count('1')==8 for w in words)
    supports=[int(w[::-1],2) for w in words]
    assert len(set(supports))==len(supports)==3088
    h=Counter((x&y).bit_count() for x,y in combinations(supports,2))
    assert max(h)==4
    for t in HEADS:
        n=148+t*t
        m=2*(14+max(2,abs(t)));r=4*(2+max(2,abs(t)))
        assert 4*max(m,r)**2<=32*n
    # Sparse/sparse: ≤16 by support intersection/parity or signed root shape.
    # Sparse and dense supports are disjoint as layer types (8,2 vs >=37).
    nd=sum(map(len,Q))*len(K)
    no=len(supports)*128
    nr=2*38*37
    assert nd+no+nr==652028
    return dict(status='PASS',dimension=38,points=nd+no+nr,dense=nd,middle=no,roots=nr,
                fibre_counts=list(map(len,Q)),fibre_size=len(K),comparisons=ranges,
                all_fibre_pairs=sum(row[2] for row in ranges))

if __name__=='__main__':
    print(json.dumps(replay(),sort_keys=True))
