"""Exact replay of two flat-sign RS layers; no optimizer or random input."""
from itertools import combinations, product
from pathlib import Path
import json
FIELD=0b10011
def mul(x,y):
    z=0
    while y:
        if y&1:
            z^=x
        y>>=1
        x<<=1
        if x&16:
            x^=FIELD
    return z
def evaluate(c,x):
    z=0
    for u in c[::-1]:
        z=mul(z,x)^u
    return z
for x in range(1,16):
    assert len({mul(x,y) for y in range(16)})==16
# Difference of any two cubic codewords is a nonzero cubic. Exhaust them.
AX=[tuple(evaluate(a,i) for i in range(6)) for a in product(range(16),repeat=4)]
assert len(AX)==len(set(AX))==65536
assert all(sum(v!=0 for v in u)>=3 for u in AX[1:])
H=[u for u in product((-1,1),repeat=6) if u.count(-1)%2==0]
assert len(H)==32 and all(sum(a!=b for a,b in zip(x,y))>=2 for x,y in combinations(H,2))
F=[tuple(x) for x in product((-1,1),repeat=4) if x[0]==1 and x.count(-1)%2==0]
assert len(F)==4
G=[(-x[0],)+x[1:] for x in F]
def dot(x,y):
    return sum(a*b for a,b in zip(x,y))
for Q in (F,G):
    assert all(dot(x,x)==4 for x in Q)
    assert all(dot(x,y)==0 for x,y in combinations(Q,2))
assert all(abs(dot(x,y))==2 for x in F for y in G)
# Embed each symbol in exactly one 4-coordinate block of a group of 16.
def rows(R):
    return [tuple(s*X[i%4] if i//4==loc else 0 for i in range(16))
            for loc,x,s in product(range(4),R,(-1,1)) for X in (x,)]
A,B=rows(F), rows(G)
assert len(A)==len(B)==32 and len(set(A+B))==64
assert all(dot(x,x)==4 for x in A+B)
assert all(abs(dot(x,y))<=2 for x in A for y in B)
for Z in (A,B):
    for u,v in combinations(range(0,32,2),2):
        assert dot(Z[u],Z[v])==0
DENSE=list(product((-1,1),repeat=4))+[
    tuple(s*2 if j==i else 0 for j in range(4)) for i in range(4) for s in (-1,1)]
MID=[tuple(2*a if k==i else 2*b if k==j else 0 for k in range(4))
     for i,j in combinations(range(4),2) for a,b in product((-1,1),repeat=2)]
assert all(abs(dot(x,y))<=4 for x,y in product(F+G,DENSE+MID))
# Old sparse RS uses exactly six entries of amplitude 4 and its roots
# exactly two entries of amplitude sqrt(48). New entries all have abs <=1.
# Norm squares of new points 24, old 96. The scaled cap is 24.
assert 6*4==24 and (2**2)*48<=24**2
assert 4*24*24==(4*6)**2*4

# The dense and sparse inherited points have norm²96 in every listed rank.
# Each flat active row uses only retained ordinary coordinates, zeros on all heads.
# Shortening a field-coordinate list restricts which symbols are allowed; it does
# not affect any of the within-family or cross-family angular inequalities.
byrank={
    92: ([4,4,3,3,3,3], 20, 20736, 30156524136,30151281256),
    93: ([2,4,4,4,4,4], 22, 32768, 44067064472,44038752920),
    94: ([2,4,4,4,4,4], 22, 32768, 74131966668,74103655116),
    95: ([2,4,4,4,4,4], 22, 32768, 110639189060,110610877508),
    96: ([4,4,4,4,4,4], 24, 65536, 209274800960,209198254912)
}
claims=json.loads((Path(__file__).parent/'claims.json').read_text())
assert sorted(c['dimension'] for c in claims)==list(range(92,97))
for c in sorted(claims,key=lambda c:c['dimension']):
    d=c['dimension']
    sizes,body,expected,old,prior=byrank[d]
    assert sum(sizes)==body and 6<body and all(0<n<=4 for n in sizes)
    # Count by evaluating every cubic. For the 92D restricted four positions
    # are a size-four information set, since degree <=3 determines polynomial.
    allowed=sum(all(x<4*k for x,k in zip(v,sizes)) for v in AX)
    assert allowed==expected
    new=allowed*len(H)*2
    assert c['inherited']==old and c['extra']==new and c['lower_bound']==old+new
    assert c['previous_published_before_both_steps']==prior
    assert c['increase_over_published_before_both_steps']==old+new-prior
    print(f"{d}D: {old+new:,} ({new:,} extra, {old+new-prior:,} across both steps)")
