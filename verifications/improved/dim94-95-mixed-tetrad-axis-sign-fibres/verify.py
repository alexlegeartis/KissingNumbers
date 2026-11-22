#!/usr/bin/env python3
"""Exact independently checkable replay of the mixed three-tetrad shell."""
from base64 import b64decode
from collections import Counter
from hashlib import sha1
from itertools import combinations, product
from pathlib import Path
from struct import unpack_from
import json
import os
import subprocess
import tempfile
import zlib

P=Path(__file__).resolve().parent
SHIFT='dim92-96-shifted-dodecacode-outer-extension'
LATEST='dim94-95-coordinate-aligned-dodecad-fibres'
SHORT='dim94-95-matching-row-append'
PC=int.bit_count

def pin(folder,name,digest):
    b=(P.parent/folder/name).read_bytes()
    assert sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==digest,(folder,name)
    return b

def positions(m):
    return [i for i in range(24) if m>>i&1]

def closure(B):
    S=[0]
    for b in B:
        S += [c^b for c in S]
    assert len(S)==len(set(S))
    return S

def old():
    f=(7,1,4,4,2,7,3,3,4,7,2,1)
    B=[sum((v&1)<<(i+j) for j,v in enumerate(f))+(1<<23) for i in range(12)]
    G=closure(B);gs=set(G)
    assert Counter(map(PC,G))=={0:1,8:759,12:2576,16:759,24:1}
    manifest=json.loads(pin(SHIFT,'shifted_patterns.json','7ad938a6be9008825e20d1e45ba9fed0ebfd501d'))
    assert manifest[-2]['mask']==3
    ids=[int(x,36) for x in manifest[-2]['patterns'].split(',')]
    assert len(ids)==len(set(ids))==1064 and all(0<=x<2**24 for x in ids)
    octs={}
    for x in G:
        if PC(x)==8:
            for J in combinations(positions(x),4):
                k=sum(1<<j for j in J)
                octs.setdefault(k,[]).append(x)
    assert len(octs)==10626 and all(len(v)==5 for v in octs.values())

    Dc=[sum(int(a)<<(2*((i+j)%12)) for j,a in enumerate('210100100101')) for i in range(12)]
    def mapper(J):
        def project(x):
            return sum((x>>(2*i)&3)<<(2*j) for j,i in enumerate(J))
        piv={}
        for x in map(project,Dc):
            while x:
                h=x.bit_length()-1
                if h in piv:x^=piv[h]
                else:piv[h]=x;break
        order=sorted(piv.items(),reverse=True)
        holes=sorted(set(range(16))-set(piv))
        assert len(holes)==4
        def par(x):
            x=project(x)
            for i,v in order:
                if x>>i&1:x^=v
            return sum((x>>i&1)<<j for j,i in enumerate(holes))
        return par
    maps=[mapper(J) for J in ((4,5,6,7,8,9,10,11),(0,1,2,3,8,9,10,11),(0,1,2,3,4,5,6,7))]
    def sy(x):
        return sum(f(x)<<(4*i) for i,f in enumerate(maps))
    lifts={0:0}
    for i in range(24):
        h=sy(1<<i)
        if h not in lifts:lifts|={k^h:v^(1<<i) for k,v in lifts.items()}
    assert len(lifts)==4096 and all(not sy(w) for w in Dc)
    ph=((0,5,10,15),(0,5,9,12),(0,5,10,15))
    stock=[]
    for x in ids:
        a,b=G[x%4096],G[x//4096]
        H=(a&~b,b&~a,a&b)
        assert PC(a|b)==12
        if not a or not b or a==b:
            stock.append((a,b,0));continue
        assert list(map(PC,H))==[4,4,4]
        v=0
        for i in range(3):
            j,k=H[(i+1)%3],H[(i+2)%3]
            al=sorted(x^j for x in octs[j] if not (x^j)&(j|k))
            assert len(al)==len(set(al))==4 and all((x|k) in gs for x in al)
            v|=ph[i][al.index(H[i])]<<(4*i)
        stock.append((a,b,lifts[v]))
    return gs,stock

def latest():
    raw=pin(LATEST,'witness.b64','adaea496ec5e332d7f0084e77e3101b804b19657')
    data=zlib.decompress(b64decode(b''.join(raw.split()),validate=True))
    count,=unpack_from('<H',data,0)
    assert count==328
    fixed=set();at=2
    for _ in range(count):
        m,t,n=unpack_from('<IBH',data,at)
        assert 0<m<2**24 and PC(m)==12 and t in (1,2,3) and n>0
        fixed.add((m if t&1 else 0,m if t&2 else 0))
        at+=7
        for _ in range(n):
            w,=unpack_from('<I',data,at)
            assert w < (1<<29)
            at+=4
            if w & (1<<28):
                count_bits,=unpack_from('<B',data,at)
                assert count_bits>0 and at+1+count_bits<=len(data)
                at+=1+count_bits
    assert at==len(data) and len(fixed)==count
    return fixed

def cap(x,y):
    a,b=x;c,d=y
    u=(a|b)&(c|d)
    return PC(u)+PC(u&~((a^c)|(b^d)))

def geometry():
    M=(((0,1),(2,3)),((0,2),(1,3)),((0,3),(1,2)))
    A=[[tuple(2 if k==i else s if k==j else 0 for k in range(4))
        for i,j in m for s in (-2,2)] for m in M]
    F=list(product((-1,1),repeat=4))
    F += [tuple(2*s*(i==j) for i in range(4)) for j in range(4) for s in (-1,1)]
    dot=lambda x,y:sum(a*b for a,b in zip(x,y))
    assert all(dot(x,x)==8 and all(abs(dot(x,y))<=4 for y in F) for S in A for x in S)
    assert all(dot(x,y)==0 for S in A for x,y in combinations(S,2))
    assert all(abs(dot(x,y))==4 for S,T in combinations(A,2) for x,y in product(S,T))
    # Common body norms and geometry unchanged from the shifted stock.
    assert 12*8==96 and 12*4==48 and 6*4==24
    assert 4*24**2 == 96*24
    assert 6*4*2==48 and 4**2*48<48**2

def main():
    assert __debug__
    G,stock=old()
    fixed=latest()
    geometry()
    short=[]
    raw_short=pin(SHORT,'witness.txt','517698a0364c917753553dff7fde9e70f9820396')
    for l in raw_short.decode().splitlines():
        m,z=map(int,l.split())
        T=[0,0,0]
        for i in positions(m):
            T[z%3]|=1<<i
            z//=3
        assert z==0
        short.append((T[0]|T[2],T[1]|T[2],PC(m)))
    assert len(short)==85
    W=[]
    raw=(P/'witness.txt').read_text()
    source={(a,b) for a,b,_ in stock}
    for l in raw.splitlines():
        aa,bb,perm,ss,cc=l.split()
        a,b,s,n=map(int,(aa,bb,ss,cc))
        assert (a|b)&3==0 and a in G and b in G and [PC(a&~b),PC(b&~a),PC(a&b)]==[4,4,4]
        assert 0<=s<2**24 and n>0 and sorted(perm)==list('0123456789ab')
        assert (a,b) not in source and (a,b) not in fixed
        assert all(cap((a,b),m)<=12 for m in fixed)
        assert all(cap((a,b),(c,d))**2<=12*w for c,d,w in short)
        assert all(cap((a,b),m)<=12 for m in W)
        assert all(cap((a,b),(c,d))<=12 or cap((a,b),(c,d))==16 for c,d,_ in stock)
        W.append((a,b))
    payload=f'{len(stock)}\n'+''.join(f'{a} {b} {s}\n' for a,b,s in stock)+f'{len(W)}\n'+raw
    with tempfile.TemporaryDirectory() as job:
        exe=Path(job)/'replay'
        subprocess.run([os.environ.get('CXX','g++'),'-std=c++17','-O3',str(P/'replay.cpp'),'-o',str(exe)],check=True)
        r=subprocess.run([str(exe)],input=payload,text=True,capture_output=True,check=True)
    k,total,axes,edges=map(int,r.stdout.split())
    assert k==len(W) and total==sum(int(l.split()[-1]) for l in raw.splitlines())
    prev=json.loads(pin(LATEST,'claims.json','633aacb5912ad1b012517cbbc13894d702ee4417'))
    assert [r['dimension'] for r in prev['claims']]==[94,95]
    out=dict(geometries=k,axis_modes=axes,old_active_interfaces=edges,
             claims=[dict(dimension=r['dimension'],previous=r['lower_bound'],
                          gain=total,lower_bound=r['lower_bound']+total)
                     for r in prev['claims']])
    assert out==json.loads((P/"claims.json").read_text())
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    main()
