#!/usr/bin/env python3
"""Exact 65D transfer of the mixed-amplitude clone construction."""
import hashlib,importlib.util,json
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
H=Path(__file__).resolve().parent
V=H.parent
ONE=(1<<64)-1
X0=[23,55,61,15,31,1,37,26,4,36,9,60,62,39,17,21,29,7,53,6,2,19,32,49,20,10,24,13,5,3,18,63,42,45,35,48,47,54,16,40,25,14,38,34,41,50,43,56,27,11,59,22,30,8,33,58,51,0,52,57,44,12,28,46]
def load(rel,h):
    p=V/rel;b=p.read_bytes()
    assert hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()==h
    return p
def module(path,n):
    s=importlib.util.spec_from_file_location(n,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    return m
def verify():
    assert __debug__
    prior=module(load('dim63-64-dense-owner-clones/verify.py',
                      'a22e5aad48fd5be5623fc844cb6348837f779aa9'),'clones')
    gf8=load('dim63-67-gf8-norm-poles/claims.json','45a826ff2f34261540adf71e299c4c4c58a36c27')
    record=next(x for x in json.loads(gf8.read_text())['dimensions'] if x['dimension']==65)
    assert record['lower_bound']==336257692 and record['previous']==336256796
    G,J,R=prior.stock()
    B=module(load('dim64-65-signed-parity/construction.py',
                   '5657a1072801fd2ded5e9fa02124976904e44e3e'),'bch')
    P,g=B.extended_bch(13)
    E,e=B.extended_bch(21)
    assert len(P)==30 and g==0x37cd0eb67 and len(E)==18 and e==0x2f30b529d3d5
    checks=[sum(((w&p).bit_count()%2)<<i for i,p in enumerate(P)) for w in P]
    X=[1]
    for i in range(62):X.append(B.gf64_multiply(X[-1],2))
    X.append(0)
    assert sorted(X)==sorted(X0)==list(range(64))
    perm=[X0.index(x) for x in X]
    def phys(w):return sum(1<<perm[j] for j in prior.bits(w))
    D0=[]
    for m in B.kernel(checks,30):
        w=0
        for j in prior.bits(m):w^=P[j]
        D0.append(w)
    dr,di=B.rref([phys(w) for w in D0],64)
    er,ei=B.rref([phys(w) for w in E],64)
    assert len(dr)==28 and all(w.bit_count()%4==0 for w in dr)
    # This affine-code relabeling acts on the BCH dense equator, E cap bodies
    # and GF8 heavy spikes simultaneously, and preserves all their pair
    # bounds. Interactions of these layers with old sparse stock depended
    # only on magnitudes (W16/W4/W9 with zero or the stated 65th coordinate).
    # No coordinate permutation of the inherited sparse stock is needed.
    Ns={h:63*4+h*h for h in (11,14)}
    for h,N in Ns.items():
        # W4/W9 and W16 without mark, plus with mark and mismatch.
        for v in (Q(4*(h+6)),Q(8,3)*(h+16),Q(64),Q(2*(h+26)),Q(16)):
            assert v*v<=16*N
        assert (h+66)**2<=16*N # other equatorial dense, cyclic d >=16
        assert Q(3,4)*(h+74)**2<=16*N # E-cap height4, d >=14 in P
        # An unretired GF8 spike has a different mark and its exterior
        # sign differs from the unflipped owner on >=14 sites, so, after
        # the one flip, distance >=13. In integer magnitudes, the common
        # ordinary contribution is <=4*38 and weights increase this
        # by 2(h-2)+32; GF8 norm squared in these units is 576.
        raw=4*38+2*(h-2)+32
        assert 4*raw*raw<=N*576
    W=set();used=0;S=[];slow=0;fast=0
    for ln in (H/'win65.txt').read_text().splitlines():
        w,M,L=map(int,ln.split())
        assert 0<=w<=ONE and 0<M<=ONE and 0<=L<=ONE and not M&L and L.bit_count()<=1
        assert w not in W and (M|L).bit_count()>=2
        assert B.in_span(w,G,J) and B.in_span(w,dr,di) and not B.in_span(w,er,ei)
        W.add(w);used|=M|L;slow|=L;fast|=M
        for mask,h in ((L,11),(M,14)):
            for a in prior.bits(mask):
                x=w^(1<<a)
                assert all(any((v&x).bit_count()%2!=t for v,t in c) for c in R[a])
                S.append((a,x,h))
    assert slow&fast==0 and len(S)==len({(a,w) for a,w,h in S})
    for (a,x,h),(b,y,k) in combinations(S,2):
        z=x^y;d=z.bit_count()
        ja=z>>a&1;jb=z>>b&1
        if a==b:
            t=4*(63-2*(d-ja))+h*k*(1-2*ja)
        else:
            t=4*(62-2*(d-ja-jb))+2*h*(1-2*ja)+2*k*(1-2*jb)
        assert t<=0 or 4*t*t<=Ns[h]*Ns[k]
    gain=len(S)-len(W)-16*used.bit_count()
    return dict(dimension=65,previous=record['lower_bound'],
                lower_bound=record['lower_bound']+gain,
                gain=gain,inserted=len(S),retired_dense_owners=len(W),
                retired_gf8_spikes=16*used.bit_count(),
                retained_gf8_spikes=16*(64-used.bit_count()),
                exact_new_pairs=len(S)*(len(S)-1)//2)
if __name__=='__main__':print(json.dumps(verify(),indent=2))
