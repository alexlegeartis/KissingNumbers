"""Exact partial-sign nonary refill of the pinned 62D weighted construction."""
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
import hashlib,importlib.util,json,struct,sys

P=Path(__file__).resolve().parent
V=Path(sys.argv[sys.argv.index('--registry')+1]).resolve() if '--registry' in sys.argv else P.parent
def pin(name,h):
    p=V/name
    assert hashlib.sha256(p.read_bytes()).hexdigest()==h,p
    return p
def load(p):
    spec=importlib.util.spec_from_file_location(p.parent.name,p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m
def stock():
    S=load(pin('dim61-66-nonary-affine-caps/verify.py','910f9bea7ac35b5995bb1f21cbf4023c65ebbea15bd02b4510415a170d4e59df'))
    old=pin('dim61-66-nonary-affine-caps/certificate.json','d2f7743669ca5ef39a6f7a58fdffc859bbc84550c815fcd89a6d958cb8b3f6de')
    C=json.loads(old.read_text())['supports']['62']
    O,*_=S.setup()
    pin('dim61-62-omitted-fibre-spikes/verify.py','29749fd2a4019dacfd67c0ce89e744e70e627c4a270aee6be08c487505a4405f')
    src=pin('dim62-coordinate-weighted-distance/verify.py','797a822e083e5ac661da29f831d163d379910d1f61d0c9b9c749667ea3004cee')
    if '--stock-only' not in sys.argv:
        z=load(src).verify()
        assert z['dimension']==62 and z['lower_bound']==126570636 and z['removed_poles']==120,z
    return O,C
def geometry(O,C):
    assert len(O)==10635 and len(set(O))==len(O) and all(w.bit_count()==16 and 0<w<1<<60 for w in O)
    assert len(C)==2891 and len(set(C))==len(C) and all(w.bit_count()==9 and 0<w<1<<62 for w in C)
    seeds=('01112','03110','21310','21132')
    S=sorted({sum([0,1,3,2][(a*int(t[(j+k)%5]))%4]<<(2*j) for j in range(5))>>1
              for t in seeds for k in range(5) for a in (-1,1)})
    assert len(S)==40
    def supp(w):return [i for i in range(62) if w>>i&1]
    def emb(w,s):return sum((s>>j&1)<<i for j,i in enumerate(supp(w)))
    A=[(w,[emb(w,s)for s in S])for w in C]
    p=Q(8,3)
    assert 9*p*p==64 and 6*2*p==32 and 3*4*p==32 and 8*p<32
    assert p*11<32 and (p*Q(25,2))**2*Q(256,317)<32**2
    new=[]; added=0; paid=set()
    W=json.loads((P/'winner.json').read_text())
    assert len(W)==60 and len({w for w,s,b in W})==len(W)
    for w,S,D in W:
        assert 0<w<1<<62 and w.bit_count()==9 and w not in C
        assert len({i^j for i,j in combinations(supp(w),2)})==36
        assert all((w&t).bit_count()<=6 for t in O)
        assert 0<len(S)<=40 and len(set(S))==len(S) and all(0<=v<512 for v in S)
        assert all((u^v).bit_count()>=3 for u,v in combinations(S,2))
        X=[emb(w,s)for s in S]
        for u,Y in new:
            t=w&u;k=t.bit_count()
            if k>4:
                assert all(k-2*((x^y)&t).bit_count()<=4 for x in X for y in Y)
        B=set()
        for j,(u,Y) in enumerate(A):
            t=w&u;k=t.bit_count()
            if k<=4:continue
            for i,y in enumerate(Y):
                if any(k-2*((x^y)&t).bit_count()>4 for x in X):
                    B.add((j,i))
        assert B==set(map(tuple,D))
        paid|=B; new.append((w,X));added+=len(S)
    result={'dimension':62,'previous':126570636,'lower_bound':126570636+added-len(paid),
            'gain':added-len(paid),'added':added,'removed':len(paid)}
    assert result==json.loads((P/'claims.json').read_text()) and result['gain']>0
    print('PASS exact 62D complete exchange:',result)
if __name__=='__main__':
    geometry(*stock())
