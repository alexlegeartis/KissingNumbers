"""Exact BCH kernel, quotient phase, size and layer-compatibility audit."""
from itertools import combinations, product
from functools import reduce
from operator import xor
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent

# All multiplications in GF(32) reduced mod x^5+x^2+1.
def mul(x,y):
    z=0
    while y:
        if y&1:z^=x
        y>>=1;x<<=1
        if x&32:x^=37
    return z
assert all(len({mul(x,y)for y in range(32)})==32 for x in range(1,32))
H=[1|(x<<1)|(mul(mul(x,x),x)<<6) for x in range(24)]
B=[H[4*i:4*i+4] for i in range(6)]
def span(xs):
    T={0}
    for a in xs:T|={x^a for x in T}
    return T
assert len(span(H))==2**11
for w in range(1,6):
    assert all(reduce(xor,S)!=0 for S in combinations(H,w))
p=[0,4,913,1759,417,772]
q=[0,1025,1987,327,1632,1281]
assert max(p+q)<2**11

def phase(i,e):
    return (p[i] if e&1 else 0)^(q[i] if e&2 else 0)

# Erase the groups in which the support words differ; neither zero nor one
# flipped sign can meet the prescribed syndrome on the shared groups.
for i,j in combinations(range(6),2):
    E=span(B[i]+B[j])
    assert len(E)==128
    V={0}|{u for k in range(6) if k not in (i,j) for u in B[k]}
    forbidden={e^u for e in E for u in V}
    for d in (1,2,3):
        assert phase(i,d)^phase(j,d) not in forbidden

assert all(len(B[i])==4 for i in range(6))
DENSE=list(product((-1,1),repeat=4))+[
    tuple(s*2 if j==i else 0 for j in range(4)) for i in range(4) for s in (-1,1)]
MID=[tuple(2*a if k==i else 2*b if k==j else 0 for k in range(4))
     for i,j in combinations(range(4),2) for a,b in product((-1,1),repeat=2)]
def dot(x,y):return sum(a*b for a,b in zip(x,y))
assert all(abs(dot(x,y))<=4 for x,y in product(product((-1,1),repeat=4),DENSE+MID))
assert 6*4==24 and 4*24**2==24*96 and 4*48<=24**2
claims=json.loads((ROOT/'claims.json').read_text())
known={
    94:(22,512,2097152,74134063820),
    95:(22,512,2097152,117083737156),
    96:(24,1024,4194304,209278995264)
}
assert sorted(c['dimension'] for c in claims)==list(known)
for r in claims:
    n=r['dimension'];rows,ns,old,base=known[n]
    groups=[rows-20,4,4,4,4,4]
    assert sum(groups)==rows and rows>6
    S=[a for a in product(*map(range,groups)) if reduce(xor,a)==0]
    assert len(S)==ns
    # Positions of the first six symbols determine support uniquely.
    # Every H affine syndrome coset has 2^(24-11) entries.
    fat=len(S)*2**(24-11)
    assert (r['previous'],r['replaced_flat'],r['new_flat'],r['increase'],r['lower_bound'])==(base,old,fat,fat-old,base+fat-old)
    print(n,r['lower_bound'],r['increase'])
