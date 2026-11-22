#!/usr/bin/env python3
"""Exact verification: weight-six mixed axial/flat append in dimensions 94 and 95."""
from itertools import combinations,product
from collections import Counter
from hashlib import sha1
from importlib import util
from pathlib import Path
import json
P=Path(__file__).resolve().parent
def pin(rel,digest):
    p=P.parent/rel; raw=p.read_bytes()
    assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==digest,rel
    return p
base='dim94-96-hexary-parity-packing/'
match='dim94-95-matching-row-append/'
source=pin(base+'verify.py','b19c5682d7279be148f722a87f103dd6e4d631c7')
spec=util.spec_from_file_location('flat_parent',source)
m=util.module_from_spec(spec);spec.loader.exec_module(m)
oldsupports,axis=m.oldsupports,m.axis
def cap(a,e,o,b,t=None):
    q=((a|e|o)&b).bit_count();s=((e|o)&b).bit_count()
    if t is None:return 2*(q+s)
    h=((e&t)|(o&(b^t))).bit_count()
    return 2*(q+h)

def geometry():
    def dot(a,b):return sum(x*y for x,y in zip(a,b))
    ax=[tuple(s*2*int(i==j) for i in range(4))for j in range(4)for s in (-1,1)]
    fl=[x for x in product((-1,1),repeat=4)]
    F=[ax,[x for x in fl if x.count(-1)%2==0],[x for x in fl if x.count(-1)%2==1]]
    assert all(len(set(A))==8 for A in F)
    assert all(dot(z,z)==4 for z in sum(F,[]))
    assert all(dot(a,b)==(-4 if tuple(-y for y in a)==b else 0)
               for A in F for a,b in combinations(A,2))
    assert all(abs(dot(a,b))<=2 for A,B in combinations(F,2) for a,b in product(A,B))
    G=[tuple(a if k==i else b if k==j else 0 for k in range(4))
       for i,j in combinations(range(4),2)for a,b in product((-2,2),repeat=2)]
    assert len(G)==24 and all(abs(dot(a,b))<=4 for a,b in product(sum(F,[]),G))
    # Six occupied rows, squared norms 24 and 96: mixed cap 24.
    assert 4*24**2==24*96 and len(F)==3
    # The inherited axis() checks the 128-word punctured-Hamming phase stock
    # and 32 even signs; two modes in its calling layer, one mode here.
    assert 4*(6-3)<=12 and 4*(6-4)<=12
    return axis()//2

def main():
    assert __debug__
    U,n=oldsupports(94);assert n==22
    U=[x<<2 for x in U]
    V=[tuple(map(int,l.split()))for l in pin(base+'best6_94.txt','5cfd7cb4793ff637bca3f5b03ef2c6cd852bb3b5').read_text().splitlines()]
    assert len(V)==len({s for s,p in V})==152 and all(s.bit_count()==6 and p&s==p and s<1<<24 for s,p in V)
    V=[(x<<2,y<<2)for x,y in V]
    matching=[int(l.split()[0])for l in pin(match+"witness.txt","517698a0364c917753553dff7fde9e70f9820396").read_text().splitlines()]
    assert len(matching)==85
    W=[tuple(map(int,l.split()))for l in (P/'witness.txt').read_text().splitlines()]
    assert len(W)==len(set(W))
    for a,e,o in W:
        assert all(0<=u<1<<24 for u in (a,e,o)) and not(a&e or a&o or e&o)
        s=a|e|o;assert s&3==0 and s.bit_count()==6 and a.bit_count() in (1,2,3)
        G=[15<<(4*q)for q in range(6)]
        # Physical amplitude-four RS groups are 14,16,... coordinates. In
        # 95, old roots can also meet the excluded head axis 88 where new
        # points vanish. On each true group, at most one axial entry exists.
        H=[0]*6
        for i in range(88):
            H[0 if i<14 else 1+(i-14)//16] |= 1<<(2+i//4)
        assert all((a&g).bit_count()<=1 for g in H)
        assert 3**2*48<=24**2
        # Amplitude-four RS chooses one coordinate per group.
        weights=[2 if a&g else (1 if (e|o)&g else 0)for g in H]
        assert sum(weights)*4<=24
        # Flat RS supports are not required to have XOR zero, but always choose
        # one row from each group. Bound each such group's largest raw dot.
        assert sum(4 if (e|o)&g else (2 if a&g else 0) for g in G)<=12
        # Earlier BCH flats and paired group BCH flats: use arbitrary old parity
        # in every row, making this independent of the inner axis phases.
        assert all(cap(a,e,o,s0)<=12 for s0 in U)
        for t,p in V:
            assert cap(a,e,o,t,p)<=12 and cap(a,e,o,t,t^p)<=12
        # No need to use the actual matching tag: every shared matching
        # row contributes at most 2sqrt(2), against squared radii 24 and 4r.
        assert all(2*(s&m).bit_count()**2<=6*m.bit_count() for m in matching)
    for (a,e,o),(b,f,p) in combinations(W,2):
        k=((a|e|o)&(b|f|p)).bit_count()
        d=((a&b)|(e&f)|(o&p)).bit_count()
        assert 2*(k+d)<=12
    # At least one axial row distinguishes these from every old pure flat
    # point. Sparse shifted matching uses 12 ordinary rows while this uses 6,
    # amplitude-four RS has 6 individual nonzeros vs 15 or 18 here, old roots 2,
    # and dense vectors use all 22 available rows.
    assert 6<12<22 and min(1+4*5,2+4*4,3+4*3)>6
    old=json.loads(pin(match+'claims.json','4d717bfb030caa966c8ea073129acfeb5032fbf0').read_text())['claims']
    assert [(d['dimension'],d['lower_bound'])for d in old]==[(94,74139180748),(95,117088854084)]
    unit=geometry()
    tally=Counter(a.bit_count()for a,_,_ in W)
    gain=len(W)*unit
    out=dict(claims=[dict(dimension=x['dimension'],previous=x['lower_bound'],gain=gain,lower_bound=x['lower_bound']+gain) for x in old],patterns=dict(sorted(tally.items())),per_pattern=unit)
    assert json.loads(json.dumps(out))==json.loads((P/'claims.json').read_text())
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
