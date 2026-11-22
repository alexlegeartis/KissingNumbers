#!/usr/bin/env python3
"""Replay the 164-axis repacking of independent six-row modes in seven dimensions."""
from collections import Counter
from hashlib import sha1
from itertools import combinations, product
from pathlib import Path
import json

P = Path(__file__).resolve().parent
D89 = 'dim89-mixed-matching-rows'
D91 = 'dim91-93-triality-matching-transfer'
D94 = 'dim94-95-mixed-triality'
M94 = 'dim94-95-matching-row-append'
H96 = 'dim94-96-hexary-parity-packing'
D96 = 'dim96-mixed-triality-six-rows'
A96 = 'dim96-affine-phase-projections'
P96 = 'dim96-partial-triality-fibres'
PC = int.bit_count

def pin(folder, name, digest):
    raw = (P.parent/folder/name).read_bytes()
    assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest() == digest, (folder,name)
    return raw.decode()

def ints(folder, name, digest):
    return [tuple(map(int,l.split())) for l in pin(folder,name,digest).splitlines()]

def claim(folder, digest):
    return json.loads(pin(folder,'claims.json',digest))

def axes():
    avoid = (0,0,0,0,4,4)
    source = [(tuple((a+b*x+c*x*x+d*x*x*x)%5 for x in range(5)) + (d,))
              for a,b,c,d in product(range(5),repeat=4)]
    assert len(source) == len(set(source)) == 625
    C = [tuple(u-(u>v) for u,v in zip(W,avoid))
         for W in source if all(u!=v for u,v in zip(W,avoid))]
    assert len(C) == len(set(C)) == 164
    assert min(sum(a!=b for a,b in zip(u,v)) for u,v in combinations(C,2)) == 3
    E = [s for s in range(64) if PC(s)%2==0]
    assert len(E)==32 and min(PC(a^b) for a,b in combinations(E,2))==2
    # Same axis word: at least two signs reverse; distinct axis words:
    # at least three row pairs are orthogonal. All rows have squared norm 4.
    assert max(4*(6-4),4*(6-3)) == 12
    return len(C),len(E)

def validate_triples(T):
    assert len(T)==len(set(T))
    assert all(all(0<=x<1<<24 for x in (a,e,o)) and not(a&e or a&o or e&o)
               and PC(a|e|o)==6 for a,e,o in T)

def tags(T):
    return [(a,e,o,a|e|o) for a,e,o in T]

def active(T):
    D = [set() for _ in T]
    for (i,(a,e,o,m)),(j,(b,f,p,n)) in combinations(enumerate(tags(T)),2):
        q = PC(m&n)
        r = PC((a&b)|(e&f)|(o&p))
        if q+r>6:
            D[i].add(j)
            D[j].add(i)
    return D

def hexary(src, shift):
    out = [(0,x<<shift,(m^x)<<shift) for m,t in src for x in (t,m^t)]
    assert all(t&m==t for m,t in src)
    return out

def main():
    assert __debug__
    n,s = axes()
    extra = (n-128)*s
    # These immutable predecessor replay programs check every cross-mode
    # pair by row maxima alone for the cited seven-dimensional component
    # inventories, except at explicitly checked 96D interfaces below.
    for folder,dig in ((D89,'4f271af79e8c4f322b5b92aafd9f7996ad2b87d2'),
                       (D91,'64dce5a923599f375a5594d4deac0a77d5393f11'),
                       (D94,'c08f62e8a665a0001ff162e35bf10ed10306f5e4'),
                       (M94,'aeef5c71dbd4d0575ba2027fa36713408cef7ca6'),
                       (H96,'b19c5682d7279be148f722a87f103dd6e4d631c7'),
                       (D96,'6b5ad3e41d156fdf606f2acb49a1168fc3050dfc'),
                       (A96,'ccd9c4c67fe4342c98756c28f3d439ce52fe5642'),
                       (P96,'47406c61d2f22c7814c1b8ec6a0400bac9394eea')):
        pin(folder,'verify.py',dig)
    claims = {}
    c89 = claim(D89,'15c6db118b6d2209aaafe25112b5daf1d29a4b46')
    w89 = ints(D89,'witness.txt','bae0e9fa8bc696ef5c23e0cc6e81fdd1be856adb')
    assert len(w89) == len(set(w89))
    assert Counter(PC(m) for m,_ in w89) == {5:2,6:859,7:3,8:61}
    assert c89['dimension']==89 and c89['rays_per_pattern']['6']==4096
    claims[89] = (c89['lower_bound'],859)

    c91 = claim(D91,'ddfcafe618488c7e8cccdbd297ef356e020f436d')
    ix = json.loads(pin(D91,'selection.json','0e2f71d7a6f2acaddc5a7a100e977958cb38dd8b'))['append_matching']
    assert ix == sorted(set(ix)) and 0<=ix[0]<=ix[-1]<len(w89)
    assert Counter(PC(w89[i][0]) for i in ix) == {5:2,6:850,7:3,8:61}
    assert [r['dimension'] for r in c91['claims']]==[91,92,93]
    for r in c91['claims']:
        claims[r['dimension']] = (r['lower_bound'],850)

    c94=claim(D94,'af80970389d9b728c85dc29a9c53931467494ad3')
    mixed94=ints(D94,'witness.txt','5c45682ae2ea39b43ab4525d30ecaa8f6df855cd')
    h94=hexary(ints(H96,'best6_94.txt','5cfd7cb4793ff637bca3f5b03ef2c6cd852bb3b5'),2)
    mat94=ints(M94,'witness.txt','517698a0364c917753553dff7fde9e70f9820396')
    assert len(mixed94)==473 and len(h94)==304 and len(mat94)==len(set(mat94))==85
    assert not any(active(mixed94+h94))
    validate_triples(mixed94+h94)
    assert all((a|e|o)&3==0 for a,e,o in mixed94+h94)
    assert Counter(PC(m) for m,_ in mat94)=={4:1,6:55,7:29}
    cmat=claim(M94,'4d717bfb030caa966c8ea073129acfeb5032fbf0')
    assert c94['per_pattern']==cmat['per_pattern']['6']==4096
    assert [r['dimension'] for r in c94['claims']]==[94,95]
    assert [r['previous'] for r in c94['claims']]==[r['lower_bound'] for r in cmat['claims']]
    for r in c94['claims']:
        claims[r['dimension']] = (r['lower_bound'],len(mixed94)+len(h94)+55)

    m96=ints(D96,'witness.txt','a00a1bdd3f82f40ed930a8b894465a3feca824ee')
    m96+=ints(D96,'witness_extra.txt','7c2e5db39c7780c573356242bd8f12ff3a82762d')
    h96=hexary(ints(H96,'best6_96.txt','bef34cfc76b36934686b7714d20299e53c1ea4ca'),0)
    a96=[tuple(map(int,l.split()[:3])) for l in pin(A96,'witness.txt','358e8701a398155022005ae200382a8a493dff15').splitlines()]
    p96=[X[:3] for X in ints(P96,'witness.txt','0926b8dedda71178708a04f450edee3f60248997')]
    assert list(map(len,(m96,h96,a96,p96)))==[1170,438,103,60]
    T=m96+h96+a96+p96
    validate_triples(T)
    D=active(T)
    assert sum(len(x) for x in D)==2*317
    free=[i for i,row in enumerate(D) if not row]
    assert sum(i<1170 for i in free)==978 and sum(1170<=i<1608 for i in free)==356
    assert len(free)==1334 and all(i<1608 for i in free)
    c96=claim(P96,'57918208eeca0a6b1cbdc6f94c35e53228a67fdb')
    assert c96['dimension']==96 and c96['fibres']==60
    claims[96]=(c96['lower_bound'],len(free))
    out=[dict(dimension=d,previous=v,repacked_modes=k,per_mode=n*s,gain=k*extra,lower_bound=v+k*extra)
         for d,(v,k) in claims.items()]
    result=dict(axis_words=n,signs_per_word=s,claims=out)
    assert result == json.loads((P/"claims.json").read_text())
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
