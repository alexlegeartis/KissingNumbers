#!/usr/bin/env python3
from itertools import combinations,product
from pathlib import Path
import json
import hashlib
from supports import hanani
P=Path(__file__).resolve().parent

def data():
    poly=[7,1,4,4,2,7,3,3,4,7,2,1]
    basis=[sum(((c&1)<<(i+j)) for j,c in enumerate(poly))+(1<<23) for i in range(12)]
    C=[0]
    for u in basis:C += [v^u for v in C]
    assert len(set(C))==4096
    def pinned(folder, file, digest):
        raw=(P.parent/folder/file).read_bytes()
        assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==digest
        return json.loads(raw)
    old='dim88-95-golay-shortening'
    extra='dim88-89-golay-variable-punctures'
    q=pinned(old,'patterns.json','1876250cf8792342b9679f6c6be667c23a4ef7b2')[0]
    assert [q['mask'],q['pivot']]==[3,3]
    s11=[(C[m%4096]&~8,C[m//4096]&~8) for m in q['ids']]
    s11 += [tuple(x) for x in pinned(extra,'patterns_extra.json','9d6e2b0cf7303b6484a9872975bccbcec4f00fdf')]
    s10 = [tuple(x) for x in pinned(old,'weight10.json','b27bed57213ff27abc2b934436d34f44bccf1c6a')]
    assert [len(set(s11)),len(set(s10))]==[982,11]
    assert all((x|y)&3==0 and (x|y).bit_count()==11 for x,y in s11)
    assert all((x|y)&3==0 and (x|y).bit_count()==10 for x,y in s10)
    return s11+s10

def decode(m,w):
    I=[i for i in range(24) if m>>i&1]
    assert len(I)==5 and min(I)>=2 and 0<=w<243
    T=[0,0,0]
    for i in I:T[w%3] |= 1<<i;w//=3
    a=T[0]|T[2]; b=T[1]|T[2]
    assert a|b==m
    return a,b

def qmatch(x,y):
    a,b=x;c,d=y
    z=(a|b)&(c|d)
    same=z&~((a^c)|(b^d))
    return z.bit_count()+same.bit_count()

def dot(x,y):return sum(a*b for a,b in zip(x,y))

def geometry():
    # Integer roots first, then multiply their entries by sqrt(2).
    # Each matching produces an independent orthogonal frame.
    MM=[[(0,1),(2,3)],[(0,2),(1,3)],[(0,3),(1,2)]]
    base=[ [tuple(1 if k==i else s if k==j else 0 for k in range(4)) for i,j in pairs for s in (1,-1)] for pairs in MM]
    assert all(dot(v,v)==2 for L in base for v in L)
    assert all(dot(x,y)==0 for L in base for x,y in combinations(L,2))
    assert all(abs(dot(x,y))<=1 for A,B in combinations(base,2) for x,y in product(A,B))
    assert len({tuple(abs(x) for x in v) for L in base for v in L})==6
    # Triality incl. old dense frame: +-2 axial, all flat.
    flat=list(product((-1,1),repeat=4))
    axial=[tuple(s*2*int(j==i) for j in range(4)) for i in range(4) for s in (1,-1)]
    assert max(abs(dot(x,y)) for x,y in product(flat+axial,sum(base,[])))==2
    assert max(abs(dot(x,y)) for x,y in product(sum(base,[]),sum(base,[])))==2
    # 20 squared radius vs 4-row triality, ordinary one-row root, old dense.
    assert 4*2*(3*2)**2<=20*16 and 4*2*2**2<=20*2 and 4*2*(5*2)**2<=20*88
    # GF(4) projective quadratics evaluated at x=0,1,2,3,infinity.
    def mul(a,b):
        x=0
        while b:
            if b&1:x^=a
            a<<=1
            if a&4:a^=7
            b>>=1
        return x
    Q=[tuple(a^mul(b,x)^mul(c,mul(x,x)) for x in range(4))+(c,) for a,b,c in product(range(4),repeat=3)]
    assert len(set(Q))==64 and min(sum(x!=y for x,y in zip(a,b)) for a,b in combinations(Q,2))==3
    signs=[s for s in range(32) if s.bit_count()%2==0]
    assert len(signs)==16 and min((s^t).bit_count() for s,t in combinations(signs,2))==2
    assert 4*(5-3)<=10 and 4*(5-2*2)<=10
    return len(Q)*len(signs)

def main():
    old=data()
    sqs=[sum(1<<(i+2) for i in x) for x in hanani(8)]
    rows=[decode(*map(int,l.split())) for l in (P/'patterns.txt').read_text().splitlines()]
    assert len(set(rows))==len(rows)
    for a,b in rows:
        assert all(((a|b)&s)!=s for s in sqs)
        for c,d in old:
            s=(c|d).bit_count()
            assert qmatch((a,b),(c,d))**2<=5*s
    for a,b in combinations(rows,2):assert qmatch(a,b)<=5
    per=geometry()
    # Four ordinary supports, one ordinary root row, 10/11 sparse rows
    # and 22 dense rows (or a pure old head) are distinct directions.
    old_bound=15178573330
    out=dict(dimension=89,previous=old_bound,patterns=len(rows),rays_per_pattern=per,addition=len(rows)*per,lower_bound=old_bound+len(rows)*per)
    return out
if __name__=='__main__':print(json.dumps(main(),indent=2))
