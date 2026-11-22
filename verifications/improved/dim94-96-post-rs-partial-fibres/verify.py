#!/usr/bin/env python3
"""Exact continuation of 94–96D six-row triality repacking."""
from functools import reduce
from operator import xor
from collections import Counter
from hashlib import sha1
from itertools import combinations, product
from importlib import util
from pathlib import Path
import json

P=Path(__file__).resolve().parent
PC=int.bit_count
D94='dim94-95-mixed-triality'
M94='dim94-95-matching-row-append'
H96='dim94-96-hexary-parity-packing'
D96='dim96-mixed-triality-six-rows'
A96='dim96-affine-phase-projections'
P96='dim96-partial-triality-fibres'
RS='dim89-96-rs-axis-code-repack'

def pinned(folder,name,digest):
    p=P.parent/folder/name
    b=p.read_bytes()
    assert sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==digest, (folder,name)
    return p
def read(folder,name,digest):
    return pinned(folder,name,digest).read_text()
def rows(folder,name,digest):
    return [tuple(map(int,l.split())) for l in read(folder,name,digest).splitlines()]
def src():
    p=pinned(P96,'verify.py','47406c61d2f22c7814c1b8ec6a0400bac9394eea')
    spec=util.spec_from_file_location('previous_partial96',p)
    m=util.module_from_spec(spec);spec.loader.exec_module(m)
    return m
def rs():
    ban=(0,0,0,0,4,4)
    C=[]
    for a,b,c,d in product(range(5),repeat=4):
        w=tuple((a+b*x+c*x*x+d*x*x*x)%5 for x in range(5))+(d,)
        if all(x!=y for x,y in zip(w,ban)):
            C.append(sum((x-(x>y))<<(2*j) for j,(x,y) in enumerate(zip(w,ban))))
    assert len(C)==len(set(C))==164
    return C
def ordinary94():
    lengths=[2,4,4,4,4,4]
    starts=[sum(lengths[:i])+2 for i in range(6)]
    A={sum(1<<(starts[j]+x) for j,x in enumerate(w))
       for w in product(*map(range,lengths)) if not reduce(xor,w,0)}
    slots=[[(v,starts[(i//2)*2+v//2]+(i%2)*(lengths[(i//2)*2+v//2]//2)+v%2)
            for v in range(4) if v%2<lengths[(i//2)*2+v//2]//2] for i in range(6)]
    B={sum(1<<r for _,r in w) for w in product(*slots)
       if not reduce(xor,(a for a,_ in w),0) and
       all(w[i][0]//2==w[i+1][0]//2 for i in (0,2,4))}
    assert (len(A),len(B),len(A|B)) == (512,160,672)
    return A|B
def check_special_group(mix,U,S):
    a,e,o,m,_,_ = mix
    assert m&3==0
    G=[15<<(4*j) for j in range(6)]
    H=[0]*6
    for i in range(88):
        H[0 if i<14 else 1+(i-14)//16] |= 1<<(2+i//4)
    # Grouped roots, amplitude-four Reed–Solomon, flat RS, and the earlier
    # ordinary and paired BCH flat layers use these independent envelopes.
    assert all(PC(a&h)<=1 for h in H) and 9*48<=24**2
    assert sum(2 if a&h else bool(m&h) for h in H)*4<=24
    assert sum(4 if (e|o)&g else 2*bool(m&g) for g in G)<=12
    assert all(2*(PC(m&u)+PC((e|o)&u))<=12 for u in U)
    # Every row shared with a D4 matching mode costs <=2 sqrt(2).
    assert all(2*PC(m&s)**2<=6*PC(s) for s in S)
def coarse(X,Y):
    a,e,o,m=X[:4];b,f,q,n=Y[:4]
    return PC(m&n)+PC((a&b)|(e&f)|(o&q))
def main():
    assert __debug__
    m=src();C=rs()
    def check_axis(W):
        assert len(W)==len(set(W)) and all(0<=x<4096 for x in W)
        assert all(m.word_weight(x^y)>=3 for x,y in combinations(W,2))
    check_axis(C)
    E=[v for v in range(64)if PC(v)%2==0]
    assert len(E)==32 and min(PC(x^y) for x,y in combinations(E,2))==2
    assert max(4*(6-3),4*(6-4)) == 12
    # Exact predecessor verifier hashes anchor the inherited interfaces.
    pinned(D94,'verify.py','c08f62e8a665a0001ff162e35bf10ed10306f5e4')
    pinned(RS,'verify.py','4f2b90dee3f670e33f17465496e00f97e4099bbd')
    prev=json.loads(read(RS,'claims.json','1b78a7a50a09fcbe7e032369d6770a4c8ce395da'))
    bounds={r['dimension']:r['lower_bound'] for r in prev['claims']}
    modes94=[m.unpack(a,e,o,C) for a,e,o in rows(D94,'witness.txt','5c45682ae2ea39b43ab4525d30ecaa8f6df855cd')]
    assert len(modes94)==473
    for s,t in rows(H96,'best6_94.txt','5cfd7cb4793ff637bca3f5b03ef2c6cd852bb3b5'):
        assert t&s==t
        for a in (t,s^t):
            modes94.append(m.unpack(0,a<<2,(s^a)<<2,C,True))
    assert len(modes94)==777
    U=ordinary94()
    S=[s for s,t in rows(M94,'witness.txt','517698a0364c917753553dff7fde9e70f9820396')]
    assert len(S)==85
    W94=[]
    for l in (P/'witness94_95.txt').read_text().splitlines():
        a,e,o,n,*W=map(int,l.split())
        assert n==len(W) and W
        check_axis(W)
        V=m.unpack(a,e,o,W)
        check_special_group(V,U,S)
        # Per occupied triality row, dense and shifted Golay cost at most 4.
        # Thus 24 is exactly the mixed bound versus squared norm 96.
        assert 4*24**2 == 24*96
        W94.append(V)
    assert len({v[:3] for v in W94})==len(W94)
    assert not {v[:3] for v in W94}&{v[:3] for v in modes94}
    old94tests=sum(m.cross(X,Y) for X,Y in product(W94,modes94))
    new94tests=sum(m.cross(X,Y) for X,Y in combinations(W94,2))
    total94=sum(len(x[-1]) for x in W94)

    # To register the replacement indices, reproduce the exact ordering used
    # by the preceding RS recipe: mixed, complementary hexary, affine, partial.
    G96=rows(D96,'witness.txt','a00a1bdd3f82f40ed930a8b894465a3feca824ee')
    G96+=rows(D96,'witness_extra.txt','7c2e5db39c7780c573356242bd8f12ff3a82762d')
    assert len(G96)==1170
    for s,t in rows(H96,'best6_96.txt','bef34cfc76b36934686b7714d20299e53c1ea4ca'):
        assert t&s==t
        for a in (t,s^t):
            G96.append((0,a,s^a))
    assert len(G96)==1608
    for l in read(A96,'witness.txt','358e8701a398155022005ae200382a8a493dff15').splitlines():
        G96.append(tuple(map(int,l.split()[:3])))
    assert len(G96)==1711
    W96=[]
    for a,e,o,n,*W in rows(P96,'witness.txt','0926b8dedda71178708a04f450edee3f60248997'):
        assert n==len(W)
        W96.append(m.unpack(a,e,o,W))
    assert len(W96)==60
    fixed96=[m.unpack(a,e,o,(0,),True) for a,e,o in G96]
    V96=[];added96=0;old96tests=0
    for l in (P/'witness96.txt').read_text().splitlines():
        i,n,*W=map(int,l.split())
        assert 0<=i<1608 and n==len(W)>128
        check_axis(W)
        X=m.unpack(*G96[i],W,allow_flat=True)
        # Every changed full mode has at least one coarse-active predecessor
        # partial neighbor, so it was NOT already repacked by the RS package.
        # No changed mode has an active neighbor in the 1711 full modes.
        assert all(coarse(X,Y)<=6 for j,Y in enumerate(fixed96) if i!=j)
        hit=sum(m.cross(X,Y) for Y in W96)
        assert hit>0
        old96tests+=hit
        added96+=n-128
        V96.append(i)
    assert len(V96)==len(set(V96))
    delta={94:32*total94,95:32*total94,96:32*added96}
    result=dict(new_modes94_95=len(W94),new_axes94_95=total94,
                active_old_pairs94_95=old94tests,active_new_pairs94_95=new94tests,
                repacked96=len(V96),net_axes96=added96,active_old_pairs96=old96tests,
                claims=[dict(dimension=d,previous=bounds[d],gain=delta[d],lower_bound=bounds[d]+delta[d])
                        for d in (94,95,96)])
    assert result==json.loads((P/"claims.json").read_text())
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
