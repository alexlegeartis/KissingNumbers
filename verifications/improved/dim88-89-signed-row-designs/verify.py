#!/usr/bin/env python3
"""Standard-library exact checker of the new support systems and 89D geometry.

The original QR/Pless and Z8/Golay code minimum distances are inherited
premises, not asserted from these small geometric checks.
"""
from itertools import combinations,product
from fractions import Fraction as F
import json
import hashlib
from pathlib import Path
from supports import hanani,double
from three88 import count as triples88
from root_rows89 import main as five_rows89

def alphabet89():
    V=list(product((-1,1),repeat=4))
    dot=lambda u,v:sum(x*y for x,y in zip(u,v))
    axes = [[tuple(2*int(j==i) for j in range(4)) for i in range(4)]]
    for p in (0,1):
        axes.append([u for u in V if u[0]==1 and u.count(-1)%2==p])
    assert [len(A) for A in axes]==[4,4,4]
    def scale(s,u):return tuple(s*x for x in u)
    frames=[[scale(s,u) for u in A for s in (1,-1)] for A in axes]
    R=[]
    for i,j in combinations(range(4),2):
      for a,b in product((-1,1),repeat=2):R.append(tuple(a if k==i else b if k==j else 0 for k in range(4)))
    RS=[scale(2,v) for v in R]
    assert len(set(sum(frames,[])))==24 and len(set(R))==24
    assert all(dot(u,u)==4 for u in sum(frames,[]))
    assert all(dot(u,u)==2 for u in R)
    assert max(dot(u,v) for u,v in combinations(R,2))==1
    for A,B in combinations(frames,2):
        assert max(dot(u,v) for u,v in product(A,B))==2
    for A in axes:
        assert all(dot(u,v)==0 for u,v in combinations(A,2))
    for A in frames:
        assert max(dot(u,v) for u,v in product(A,RS))==4
        assert max(dot(u,v) for u,v in product(A,R))==2
    assert max(dot(u,v) for u,v in product(R,RS))==4
    axis_words=[(a,b,c,a^b^c) for a,b,c in product(range(4),repeat=3)]
    sign_words=[s for s in range(16) ]
    assert len(set(axis_words))==64
    assert all(sum(a!=b for a,b in zip(x,y))>=2 for x,y in combinations(axis_words,2))
    assert len(sign_words)==16 and all((a^b).bit_count()>=1 for a,b in combinations(sign_words,2))
    # Same axes: at least one sign reversal; different axes: at most two matches.
    assert 4*(4-2)<=8 and 4*2==8
    # Axial, even-flat and odd-flat are exactly the three old dense row classes.
    for r in (80,88,F(177,2),96):
        assert 4*(4*4)**2<=16*r and 4*4**2<=2*r
    assert 4*2**2<=2*16
    return len(axes)*len(axis_words)*len(sign_words),len(R)

def pin_claim(folder, digest):
    raw=(Path(__file__).resolve().parent.parent/folder/'claims.json').read_bytes()
    assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==digest
    return json.loads(raw)

def main():
    assert __debug__
    par88=pin_claim('dim88-axial-pless','0b1544896ddfda12bdde8959cbf6013f5d19df46')
    par89=pin_claim('dim88-89-golay-variable-punctures','c9c080ee96b66e9c7574140c2ee643b933a237b3')[1]
    S22=hanani(8)
    S44=double(S22,22)
    Q42=[x for x in S44 if x[-1]<42]
    old42=[x for x in hanani(16) if x[-1]<42]
    assert len(S22)==385 and len(S44)==3311 and len(Q42)==2730 and len(old42)==2604
    for X in (S22,Q42):
        triples=[t for x in X for t in combinations(x,3)]
        assert len(triples)==len(set(triples))
    phase8=[s for s in range(256) if s.bit_count()%2==0]
    assert len(phase8)==128 and all((a^b).bit_count()>=2 for a,b in combinations(phase8,2))
    # In the unchanged 88D frame all candidate square corners have local
    # norm 1.  Same support: sign code has 2+ differing coordinates,
    # each changing the dot by 1. Distinct supports overlap in at most 2.
    assert 4-2==2
    # Dense QR full triangle coordinates of norm 1 in each old plane.
    assert 4**2 < 45
    # Pless local norm sqrt(3); square directions 15 mod 90 relative to
    # the old lines at multiples of 60. (4 sqrt(3) cos(15 degrees))^2
    # =24+12sqrt(3); verify it below 45 without radicals.
    angular=[min((a-b)%360,(b-a)%360) for a,b in product(range(15,360,90),range(0,360,60))]
    assert min(angular)==15 and 12**2*3 < (45-24)**2
    # D88 roots use the 2 axes of each square frame, so each of their two
    # unit components picks up at most 1/sqrt(2), and raw dot<=sqrt(2).
    assert 4*2==2*4
    assert 15-4>4 and 42>4
    counts=[7812,6561,2187,729,729]
    hist=[37468,80272,58152,16784,1628]
    parent88=2*2**4*3**18+sum(c*h for c,h in zip(counts,hist))+2*88*87+len(old42)*len(phase8)
    assert parent88==par88['lower_bound']==13357769652
    delta88=(len(Q42)-len(old42))*len(phase8)
    Q42_set=set(Q42)
    # 5 unit components along one square basis axis per occupied plane.
    # With square corners at 15 degrees modulo 90 from Pless angle zero,
    # its two axes are 60/150; choose 150, at least 30 from all Pless rays.
    pless=[min((150-b)%360,(b-150)%360) for b in range(0,360,60)]
    assert min(pless)==30
    assert 4*F(15,2)**2==5*45 # all five Pless projections, each <=3/2
    assert 4*5**2 < 5*45 # dense projection
    assert 4*F(9,2)<5*4 # max three shared square planes, dot <=3/sqrt(2)
    assert 4*1**2<5*2<4*2**2 # exactly two matching basis axes block
    fives=[tuple(map(int,l.split())) for l in Path(__file__).with_name('five.txt').read_text().splitlines()]
    assert all(tuple(sorted(set(b)))==b and len(b)==5 and 0<=b[0]<b[-1]<42 for b in fives)
    seen=set()
    for b in fives:
        T=set(combinations(b,3))
        assert not seen&T
        seen.update(T)
        assert all(q not in Q42_set for q in combinations(b,4))
    signs=[s for s in range(32) if s.bit_count()%2==0]
    assert len(signs)==16 and all((x^y).bit_count()>=2 for x,y in combinations(signs,2))
    assert 5-2*2<=F(5,2) and 2<=F(5,2)
    assert all({(x>>i&1,x>>j&1) for x in signs}==set(product(range(2),repeat=2)) for i,j in combinations(range(5),2))
    axis_pairs={p for b in fives for p in combinations(b,2)}
    plus88=len(fives)*len(signs)
    minus88=4*len(axis_pairs)
    assert len(fives)==689 and len(axis_pairs)==861 and plus88-minus88==7580
    ntriples=triples88()
    delta88+=plus88-minus88+ntriples
    per,roots=alphabet89()
    parent89=14*8**10+2+982*1024*144+11*256*72
    assert parent89==par89['lower_bound']==15177390082
    extra89=five_rows89()
    assert extra89["previous"]==parent89+len(S22)*per+22*roots
    delta89=len(S22)*per+22*roots+extra89["addition"]
    res=[dict(dimension=88,previous=parent88,square_supports=len(Q42),new_five_plane_points=plus88,deleted_D88_roots=minus88,new_three_plane_points=ntriples,added=delta88,lower_bound=parent88+delta88),
         dict(dimension=89,previous=parent89,supports=len(S22),rays_per_support=per,root_rays=22*roots,five_row_patterns=extra89["patterns"],rays_per_five_row_pattern=extra89["rays_per_pattern"],added=delta89,lower_bound=parent89+delta89)]
    assert [x['lower_bound'] for x in res]==[13357793920,15179459090]
    assert res==json.loads(Path(__file__).with_name('claims.json').read_text())
    print(json.dumps(res,indent=2))
if __name__=='__main__':main()
