#!/usr/bin/env python3
"""Exact 65D transfer of the mixed-amplitude clone construction."""
import gzip,hashlib,importlib.util,json
import numpy as np
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
H=Path(__file__).resolve().parent
V=H.parent
ONE=(1<<64)-1
X0=[23,55,61,15,31,1,37,26,4,36,9,60,62,39,17,21,29,7,53,6,2,19,32,49,20,10,24,13,5,3,18,63,42,45,35,48,47,54,16,40,25,14,38,34,41,50,43,56,27,11,59,22,30,8,33,58,51,0,52,57,44,12,28,46]
def load(rel,h):
    p=V/rel;b=p.read_bytes()
    sig=hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest() if len(h)==40 else hashlib.sha256(b).hexdigest()
    assert sig==h
    return p
def module(path,n):
    s=importlib.util.spec_from_file_location(n,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    return m
def verify(d):
    assert __debug__
    prior=module(load('dim63-64-dense-owner-clones/verify.py',
                      'a22e5aad48fd5be5623fc844cb6348837f779aa9'),'clones')
    gf8=load('dim63-67-gf8-norm-poles/claims.json','45a826ff2f34261540adf71e299c4c4c58a36c27')
    record=next(x for x in json.loads(gf8.read_text())['dimensions'] if x['dimension']==d)
    assert record['lower_bound']=={66:342390470,67:355196300}[d]
    G,J,_=prior.stock()
    O=load('dim61-68-95-96-code-layers/source/dim61_96/code64_supports.npz','cdbd5e5a7bbda9ddb84ee86052e2d4d00ccc331f0188de958a69ed83b0c529f6')
    original=set(map(int,np.load(O,allow_pickle=False)['words']))
    assert len(original)==30828 and all((w&g).bit_count()%2==0 for w in original for g in G)
    if d==66:
        c=json.loads(load('dim61-62-66-affine-reconstruction/dim66/certificate.json',
                          '9eb80d0529de4df4eb86f4366e2ecabcb779ecf0').read_text())
        g=c['family']
    else:
        c=json.loads(load('dim61-67-compatible-sign-families/data/subcode67.json',
                          'fcfe182b0d6e2becac726f11f06560b24f196197').read_text())
        f=load('dim61-68-95-96-code-layers/source/shared_parity/expanded67_68.npz',
               '92b815f8cea3197578fae34336095f418abd754a')
        z=np.load(f,allow_pickle=False)
        g={kind:[dict(support=int(s),checks=list(map(int,h)),shift=int(t)) for s,h,t in
            zip(z[kind+'_supports'],z[kind+'_duals'],z[kind+'_shifts'])] for kind in ('old','new')}
    # Modified parent phases replace the corresponding earlier phase.
    q={r['support']:r for r in g['old']+g['new']+c['modified_parent_classes']+c['new_classes']}
    assert all(w.bit_count()==16 and (w>>64).bit_count() in (0,2) for w in q)
    q={w:r for w,r in q.items() if w>>64==0}
    assert len(original|q.keys())=={66:32794,67:32669}[d]
    R=[[] for _ in range(64)]
    for r in q.values():
        A=list(prior.bits(r['support']))
        h=[sum(1<<A[j] for j in prior.bits(i)) for i in r['checks']]
        t=[(i&r['shift']).bit_count()%2 for i in r['checks']]
        for a in A:R[a].append(tuple(zip(h,t)))
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
    assert len(dr)==28 and all(w.bit_count()%4==0 for w in dr)
    pr,pi=B.rref([phys(w) for w in P],64)
    shift=phys(0x800000037cd0eb67)
    assert B.in_span(shift,pr,pi) and not B.in_span(shift,dr,di)
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
        assert Q(37,64)*(h+74)**2<=16*N # Kerdock cap, exterior P word
        # An unretired GF8 spike has a different mark and its exterior
        # sign differs from the unflipped owner on >=14 sites, so, after
        # the one flip, distance >=13. In integer magnitudes, the common
        # ordinary contribution is <=4*38 and weights increase this
        # by 2(h-2)+32; GF8 norm squared in these units is 576.
        raw=4*38+2*(h-2)+32
        assert 4*raw*raw<=N*576
    W=set();used=0;S=[];slow=0;fast=0
    for ln in gzip.decompress((H/f'win{d}.txt.gz').read_bytes()).decode('ascii').splitlines():
        w,M,L=map(int,ln.split())
        assert 0<=w<=ONE and 0<M<=ONE and 0<=L<=ONE and not M&L and L.bit_count()<=1
        assert w not in W and (M|L).bit_count()>=2
        assert B.in_span(w,G,J) and B.in_span(w,dr,di) 
        W.add(w);used|=M|L;slow|=L;fast|=M
        for mask,h in ((L,11),(M,14)):
            for a in prior.bits(mask):
                x=w^(1<<a)
                assert all(any((v&x).bit_count()%2!=t for v,t in c) for c in R[a])
                S.append((a,x,h))
    assert slow&fast==0 and len(S)==len({(a,w) for a,w,h in S})
    for h in Ns:
        for k in Ns:
            assert 4*(136+2*h+2*k)**2<=Ns[h]*Ns[k]
    for h,k in ((11,14),(14,14)):
        assert 4*(248-2*h-2*k)**2<=Ns[h]*Ns[k]
    grouped=[[] for _ in range(64)]
    for a,x,h in S:grouped[a].append((x,h))
    for a,points in enumerate(grouped):
        for (x,h),(y,k) in combinations(points,2):
            z=x^y;da=z>>a&1
            t=4*(63-2*(z.bit_count()-da))+h*k*(1-2*da)
            assert t<=0 or 4*t*t<=Ns[h]*Ns[k]
    gain=len(S)-len(W)-16*used.bit_count()
    return dict(dimension=d,previous=record['lower_bound'],
                lower_bound=record['lower_bound']+gain,
                gain=gain,inserted=len(S),retired_dense_owners=len(W),
                retired_gf8_spikes=16*used.bit_count(),
                retained_gf8_spikes=16*(64-used.bit_count()),
                exact_new_pairs=len(S)*(len(S)-1)//2)
if __name__=='__main__':print(json.dumps({'dimensions':[verify(d) for d in (66,67) if (H/f'win{d}.txt.gz').exists()]},indent=2))
