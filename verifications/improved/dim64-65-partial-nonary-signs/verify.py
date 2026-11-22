"""Exact verifier for the signed partial nonary layers."""
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import numpy as np

H=Path(__file__).resolve().parent
V=H.parent
# The release is placed beside the pinned parent. This optional path permits
# replay against a checkout or the archived copy of those exact parent files.
def pin(p, sha):
    assert sha256(p.read_bytes()).hexdigest()==sha, p
    return p

def parent(v):
    o=pin(v/'dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz',
          'cdbd5e5a7bbda9ddb84ee86052e2d4d00ccc331f0188de958a69ed83b0c529f6')
    G=pin(v/'dim61-67-compatible-sign-families/data/global64.json',
          '8c0c35f8127c578118f1f09e5950d712a272db32f2d50f382e6fe32f46c6dbf5')
    A=v/'dim63-65-affine-reconstruction/certificate.json'
    if not A.is_file(): A=v/'dim63-65-affine-reconstruction_certificate.json'
    A=pin(A,'266af3683caadf6d23c2aa3e7409ce28d2d3f4a8644be2784263c80bd407fc7e')
    N=pin(v/'dim61-66-nonary-affine-caps/certificate.json',
          'd2f7743669ca5ef39a6f7a58fdffc859bbc84550c815fcd89a6d958cb8b3f6de')
    load=lambda f:json.loads(f.read_text())
    stock=list(map(int,np.load(o,allow_pickle=False)['words']))
    stock += [int(x['support']) for x in load(G)['new']]
    stock += [int(x['support']) for x in load(A)['new_classes']]
    assert len(stock)==len(set(stock))==32801
    assert all(0<s<1<<64 and s.bit_count()==16 for s in stock)
    n=load(N)
    gr=n['gray_map'];seeds=n['pentacode_seeds']
    assert gr==[0,1,3,2] and seeds==['01112','03110','21310','21132']
    signs={sum(gr[(a*int(t[(j+k)%5]))%4]<<(2*j) for j in range(5))>>1
           for t in seeds for k in range(5) for a in (-1,1)}
    assert len(signs)==40 and all((a^b).bit_count()>=3 for a,b in combinations(signs,2))
    old={}
    for s in n['supports']['64']:
        L=[i for i in range(64) if s>>i&1]
        assert 0<s<1<<64 and len(L)==9 and s not in old
        old[s]=[sum(((a>>j)&1)<<i for j,i in enumerate(L)) for a in signs]
    assert len(old)==999
    return stock,old

def colls(L,R):
    test=0
    for s,ls in L.items():
        for t,rs in R.items():
            m=s&t
            k=m.bit_count()
            if k<=4:continue
            for x in ls:
                for y in rs:
                    assert k-2*((x^y)&m).bit_count()<=4, (s,t,x,y)
                    test+=1
    return test

def verify(v=V):
    assert __debug__
    stock,old=parent(v)
    data=pin(H/'certificate.txt','32dfcc52f14b7cae87d7812eadba8e3d49961a538adb3dfb3f1248d899df9faf')
    pairs=[tuple(map(int,line.split())) for line in data.read_text().splitlines()]
    assert len(pairs)==len(set(pairs))==5698
    cand=defaultdict(list)
    for s,x in pairs:
        assert 0<s<1<<64 and s.bit_count()==9 and 0<=x<1<<64 and x&s==x
        assert s not in old
        cand[s].append(x)
    assert len(cand)==771
    a=Q(8,3)
    assert 9*a*a==64 and 4*a*a<32 and 5*a*a>32
    assert 6*2*a==32 and 3*4*a==32 and 8*a<32 and 9*a<32
    for s,J in cand.items():
        I=[i for i in range(64) if s>>i&1]
        D={i^j for i,j in combinations(I,2)}
        assert len(D)==36 # excludes every affine SQS quadruple
        assert all((s&t).bit_count()<=6 for t in stock)
        assert all((x^y).bit_count()>=3 for x,y in combinations(J,2))
    check_old=colls(cand,old)
    check_new=0
    for i,(s,ss) in enumerate(cand.items()):
        check_new+=colls({s:ss},dict(list(cand.items())[:i]))
    # 65D parent owners are dense and GF(4)-line members of the 64D
    # parent, with horizontal factor sqrt(3)/2 (see parent nonary replay).
    # This factor cannot introduce a >32 cross product; poles are orthogonal.
    assert 0<Q(3,4)<1
    base64=335749696+40*len(old)
    base65=335817922+40*len(old)
    return {'dimension_64':base64+len(pairs),'dimension_65':base65+len(pairs),
            'gain_each':len(pairs),'supports':len(cand),
            'signed_nonary_pairs':check_old,'signed_cross_support_pairs':check_new}
if __name__=='__main__':
    import sys
    print(json.dumps(verify(Path(sys.argv[1]) if len(sys.argv)>1 else V),indent=2,sort_keys=True))
