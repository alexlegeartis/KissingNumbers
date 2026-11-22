#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import hashlib,json
from importlib import util
P=Path(__file__).resolve().parent

def pin(rel,digest):
    p=P.parent/rel; raw=p.read_bytes()
    assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==digest,rel
    return p
base='dim94-96-hexary-parity-packing/'
shift='dim92-96-shifted-dodecacode-outer-extension/'
module=pin(base+'verify.py','b19c5682d7279be148f722a87f103dd6e4d631c7')
spec=util.spec_from_file_location('flat_parent',module)
m=util.module_from_spec(spec);spec.loader.exec_module(m)
oldsupports=m.oldsupports
def decode(m,w):
    I=[i for i in range(24) if m>>i&1]
    assert 0<m<1<<24 and len(I) in (4,6,7) and min(I)>=2 and 0<=w<3**len(I)
    T=[0,0,0]
    for i in I:T[w%3]|=1<<i;w//=3
    return m,T[0]|T[2],T[1]|T[2]
def qmatch(a,b):
    M,x,y=a;N,s,t=b
    o=M&N
    return o.bit_count()+(o&~((x^s)|(y^t))).bit_count()
def phase(x,C):
    return sum((int((x&c).bit_count()%2)<<i) for i,c in enumerate(C))
def geometry():
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    M=[((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2))]
    A=[[tuple(1 if k==i else s if k==j else 0 for k in range(4))
         for i,j in pairs for s in (-1,1)] for pairs in M]
    assert len(set(sum(A,[])))==12
    assert all(dot(v,v)==2 for v in sum(A,[]))
    assert all(dot(a,b)==0 for L in A for a,b in combinations(L,2))
    assert all(abs(dot(a,b))<=1 for L,K in combinations(A,2) for a,b in product(L,K))
    oldrows=list(product((-1,1),repeat=4))+[tuple(s*2*int(i==j) for i in range(4))
                            for j in range(4) for s in (-1,1)]
    assert all(abs(dot(a,b))<=2 for a,b in product(sum(A,[]),oldrows))
    # Matching rows are sqrt(2) times these integer axes. Shifted Golay rows
    # are twice the same axes. The cross row cap is 2sqrt(2) or 4sqrt(2)
    # according to unequal or equal pairing. A dense or flat row contributes <=2sqrt(2).
    Q4=[(a,b,c,a^b^c) for a,b,c in product(range(4),repeat=3)]
    Q7=[]
    for x in range(16):
        a=phase(x>>1,[1,2,3,4,5,6,7])^(127 if x&1 else 0)
        b=0
        for i,f in enumerate((15,0,8,1)):
            if x>>i&1:b^=f
        for y in range(8):
            c=b^phase(y,[7,3,4,6,2,1,5])
            Q7.append(tuple((a>>i&1)+2*(c>>i&1) for i in range(7)))
    Q6=[v[:6] for v in Q7]
    out={}
    for w,Q,ct,md,eps in ((4,Q4,64,2,1),(6,Q6,128,3,2),(7,Q7,128,4,2)):
        assert len(Q)==len(set(Q))==ct
        assert min(sum(a!=b for a,b in zip(y,z)) for y,z in combinations(Q,2))==md
        signs=[s for s in range(1<<w) if w==4 or s.bit_count()%2==0]
        assert min((a^b).bit_count() for a,b in combinations(signs,2))==eps
        assert 4*(w-md)<=2*w and 4*(w-2*eps)<=2*w
        # Norm squared 4w versus dense radius² 96, raw cap 2sqrt(2)*w.
        assert 8*w*w<=96*w
        # At most two root factors each sqrt(48)*sqrt(2).
        assert 4*48*2<=96*w
        out[w]=len(Q)*len(signs)
    return out

def read_prior():
    C=[0]
    for i in range(12):
        v=sum(((c&1)<<(i+j))for j,c in enumerate([7,1,4,4,2,7,3,3,4,7,2,1]))+(1<<23)
        C += [v^x for x in C]
    assert len(set(C))==4096
    shifted=json.loads(pin(shift+'shifted_patterns.json','7ad938a6be9008825e20d1e45ba9fed0ebfd501d').read_text())[1]
    assert shifted['mask']==3
    old=[]
    for v in map(lambda x:int(x,36),shifted['patterns'].split(',')):
        a,b=C[v%4096],C[v//4096]
        assert (a|b)&3==0 and (a|b).bit_count()==12
        old.append((a|b,a,b))
    assert len(old)==len(set(old))==1064
    U,n=oldsupports(94)
    assert n==22
    V=[int(l.split()[0]) for l in pin(base+'best6_94.txt','5cfd7cb4793ff637bca3f5b03ef2c6cd852bb3b5').read_text().splitlines()]
    assert len(V)==len(set(V))==152 and all(x.bit_count()==6 and x<1<<22 for x in V)
    return old,[v<<2 for v in U|set(V)]

def main():
    assert __debug__
    old,flats=read_prior()
    S=[decode(*map(int,l.split()))for l in (P/'witness.txt').read_text().splitlines()]
    assert len(S)==len(set(S))
    cnt=Counter()
    for z in S:
        m=z[0];w=m.bit_count();cnt[w]+=1
        I=[i for i in range(2,24)if m>>i&1]
        # Dimension-94 groups in physical coordinates have sizes 14,16,...,16.
        # Dimension 95 excludes coordinate 88, after the 88 ordinary coordinates.
        # Thus these are the exact same group intersections in the body for both.
        def group(j):return 0 if j<14 else 1+(j-14)//16
        amps={group(4*(i-2)+j) for i in I for j in range(4)}
        assert len(amps)**2<=3*w  # 4sqrt(2) per old RS row; square <=96w.
        # Earlier flat RS construction uses one block from each old four-row group.
        assert 2*len({i//4 for i in I})**2<=6*w
        for f in flats:assert 2*(f&m).bit_count()**2<=6*w
        for o in old:assert qmatch(z,o)**2<=12*w
    for a,b in combinations(S,2):
        assert qmatch(a,b)**2<=a[0].bit_count()*b[0].bit_count()
    inner=geometry()
    inc=json.loads(pin(base+'claims.json','766a432c247d00a349629e5daaea0686e2e55f47').read_text())
    delta=sum(inner[k]*z for k,z in cnt.items())
    assert cnt=={4:1,6:55,7:29} and delta==463872
    result=[dict(dimension=r['dimension'],previous=r['lower_bound'],lower_bound=r['lower_bound']+delta,gain=delta)
            for r in inc if r['dimension'] in (94,95)]
    out=dict(claims=result,patterns=dict(sorted(cnt.items())),per_pattern=inner)
    assert json.loads(json.dumps(out))==json.loads((P/'claims.json').read_text())
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
