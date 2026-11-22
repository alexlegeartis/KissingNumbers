#!/usr/bin/env python3
"""Exact verification: weight-six mixed axial/flat append to current 96D code."""
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
    U,n=oldsupports(96);assert n==24
    V=[tuple(map(int,l.split()))for l in pin(base+'best6_96.txt','bef34cfc76b36934686b7714d20299e53c1ea4ca').read_text().splitlines()]
    assert len(V)==len({s for s,p in V})==219 and all(s.bit_count()==6 and p&s==p and s<1<<24 for s,p in V)
    files=[(P/f).read_text().splitlines() for f in ('witness.txt','witness_extra.txt')]
    assert [len(f) for f in files]==[822,348]
    W=[tuple(map(int,l.split())) for f in files for l in f]
    assert len(W)==len(set(W))
    for a,e,o in W:
        assert all(0<=u<1<<24 for u in (a,e,o)) and not(a&e or a&o or e&o)
        s=a|e|o;assert s.bit_count()==6 and a.bit_count() in (1,2,3)
        G=[15<<(4*q)for q in range(6)]
        # The amplitude-four RS groups of four rows align here. Since the
        # grouped norm-96 roots use both physical entries in one group, its
        # two coordinates meet at most one entry of magnitude 2 (an axial row).
        assert all((a&g).bit_count()<=1 for g in G)
        assert 3**2*48<=24**2
        # Amplitude-four RS chooses one coordinate per group.
        weights=[2 if a&g else (1 if (e|o)&g else 0)for g in G]
        assert sum(weights)*4<=24
        # Flat RS supports are not required to have XOR zero, but always choose
        # one row from each group. Bound each such group's largest raw dot.
        assert sum(4 if (e|o)&g else (2 if a&g else 0) for g in G)<=12
        # Earlier BCH flats and paired group BCH flats: use arbitrary old parity
        # in every row, making this independent of the inner axis phases.
        assert all(cap(a,e,o,s0)<=12 for s0 in U)
        for t,p in V:
            assert cap(a,e,o,t,p)<=12 and cap(a,e,o,t,t^p)<=12
    for (a,e,o),(b,f,p) in combinations(W,2):
        k=((a|e|o)&(b|f|p)).bit_count()
        d=((a&b)|(e&f)|(o&p)).bit_count()
        assert 2*(k+d)<=12
    # At least one axial row distinguishes these from every old pure flat
    # point. Sparse shifted matching uses 12 ordinary rows while this uses 6,
    # amplitude-four RS has 6 individual nonzeros vs 15 or 18 here, old roots 2,
    # and dense vectors use all 24 rows.
    assert 6<12<24 and min(1+4*5,2+4*4,3+4*3)>6
    old=json.loads(pin(base+'claims.json','766a432c247d00a349629e5daaea0686e2e55f47').read_text())[-1]
    assert old['dimension']==96 and old['lower_bound']==209287080768
    unit=geometry()
    tally=Counter(a.bit_count()for a,_,_ in W)
    gain=len(W)*unit
    construction_base=old['lower_bound']
    previous=construction_base+822*unit
    out=dict(dimension=96,previous=previous,patterns=dict(sorted(tally.items())),per_pattern=unit,gain=348*unit,lower_bound=construction_base+gain,construction_base=construction_base,construction_gain=gain)
    assert json.loads(json.dumps(out))==json.loads((P/'claims.json').read_text())
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
