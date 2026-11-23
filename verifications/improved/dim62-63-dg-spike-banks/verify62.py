"""Independently replay the 62D new shortened banks and spike family."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from contextlib import redirect_stdout
from collections import Counter,defaultdict
from fractions import Fraction as Q
from io import StringIO
import base64,gzip,struct,json,runpy
def run(here,parent,oldpath,DG,A,cube):
    R=parent/'dim62-two-mark-diagonals'
    P={'verify.py':'e19538b6ab77d37377a885f3b0237139d3cf05024fdee522285c0730401a37f0',
       'merged62.txt.gz.b64':'46e3d3d7042064fd8b8d109eed5a5212d31a30ade783709d60c96997dc152fb0',
       'claims.json':'9370b3d60375ccd16ff0c451469341ccddd5273874ac028167f7fbe8da60b1f6'}
    assert all(sha256((R/p).read_bytes()).hexdigest()==d for p,d in P.items())
    with redirect_stdout(StringIO()):
        D=runpy.run_path(str(R/'verify.py'))
    prior=D['claims']['lower_bounds']['62'];assert prior==172134412
    M=(0,9,20,41);mask=sum(1<<i for i in M)
    heads=((0,1),(1,0),(0,-1),(-1,0))
    rect=((2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1))
    labels=(0,1,3,7,6,4)
    def gt(a,b,r):
        if b==0:return a>r
        if b>0:return a>=r or 3*b*b>(r-a)**2
        return a>r and (a-r)**2>3*b*b
    def fmt(s):
        t=[i for i,k in enumerate(M)if (s>>k)&1]
        u,v=(sum(heads[i][j]for i in t)for j in range(2))
        assert u*u+v*v==len(t) and len(t)<=2 and s.bit_count()==16
        q={sum((w>>j)&1 for j in t)%2 for w,(a,b) in zip(labels,rect)if gt(u*a,v*b,len(t))}
        assert len(q)==1 if t else not q
        assert all(not gt(u*a,v*b,len(t)+2) for a,b in rect)
        return s&~mask,u,v,(1-next(iter(q))) if t else 1,len(t)
    C=sorted(s for s in D['old'] if (s&mask).bit_count()<=1)
    old=set(C);E=[]
    src=gzip.decompress(base64.b64decode((oldpath/'witness62.bin.gz.b64').read_bytes()))
    k=0
    while k<len(src):
        s,d=struct.unpack_from('<QB',src,k);k+=9
        assert s not in old;old.add(s)
        pos=[1<<j for j in range(64)if (s&~mask)>>j&1]
        for _ in range(d):
            i,w=struct.unpack_from('<HH',src,k);k+=4
            g=sum(pos[j]for j in range(len(pos)) if w>>j&1)
            E.append((s,C[i],g,g))
    src=gzip.decompress(base64.b64decode((R/'merged62.txt.gz.b64').read_bytes())).splitlines()
    a,b=map(int,src[0].split())
    old.update(map(int,src[1:1+a]))
    for l in src[1+a:]:
        x,y,g,t=map(int,l.split());E.append((x,y,g,g))
    src=gzip.decompress(base64.b64decode((here/'short62.txt.gz.b64').read_bytes())).splitlines()
    a,b=map(int,src[0].split())
    extra={}
    for l in src[1:1+a]:
        s,u,v=map(int,l.split());f=fmt(s)
        assert s not in extra and s not in old and s in D['old'] and f[1:3]==(u,v) and f[4]==2
        extra[s]=f
    check=set()
    for l in src[1+a:]:
        x,y,g,f,t=map(int,l.split())
        r=(min(x,y),max(x,y))
        assert x in extra and r not in check
        check.add(r);E.append((x,y,g,f))
    assert len(extra)==a==213 and len(check)==b==223
    U={s:fmt(s) for s in old};U.update(extra)
    J=json.loads((here/'sign62.json').read_text());bits=int(J['bits'],16)
    assert J['num']==len(E)==12293 and J['short_banks']==a and bits>>len(E)==0
    eq=defaultdict(list)
    for i,(x,y,g,f) in enumerate(E):
        assert x!=y and g&x==g and f&y==f and g.bit_count()==f.bit_count()==8
        assert g&~mask==f&~mask and (g&~mask).bit_count()>=6
        t=bits>>i&1
        eq[x].append((g,t));eq[y].append((f,t^1))
    count=0
    for x, y in list(combinations(extra,2)):
        a,f,g,_,_=U[x];b,h,j,_,_=U[y]
        z=(a&b).bit_count()+f*h+g*j
        if z>8:
            assert z<=10 and (min(x,y),max(x,y)) in check
            count+=1
    for x in extra:
        a,f,g,_,_=U[x]
        for y in old:
            b,h,j,_,_=U[y]
            z=(a&b).bit_count()+f*h+g*j
            if z>8:
                assert z<=10 and (min(x,y),max(x,y)) in check
                count+=1
    assert count==len(check)
    rows=D['rows']
    for s,(ordinary,_,_,p,m) in U.items():
        assert all(not ((s&w).bit_count()&1) for w in rows)
        if m<2: assert p!=((s&mask&((1<<M[0])|(1<<M[2]))).bit_count()&1)
        H={0:0}
        for g,t in [(s,p)]+eq[s]:
            if g in H:assert H[g]==t
            else:H.update({a^g:b^t for a,b in list(H.items())})
        assert len(H) in (2,4,8,16,32)
        assert all(g in (0,s) or g.bit_count()==8 for g in H)
    # All source mark parities remain fixed. Protect every new support
    # and all two-mark banks by exact magnitude against the spike layer.
    L=5
    chosen=sum(1<<j for i,j in enumerate(M)if (L>>i)&1)
    F=Counter(x&mask for x in A)
    assert len(F)==16 and set(F.values())=={8}
    for c in cube[:60]:
        assert all(F[(x^c^chosen)&mask]==8 for x in DG)
    assert len(set(cube[:60]))==60
    assert max((L^label).bit_count() for label in labels)==2
    ampl=Q(22,5);N=59+ampl**2
    gate=lambda dot,u,v:dot<=0 or 4*dot*dot<=u*v
    assert gate(60-2*24+ampl**2-1,N,N)
    assert gate(60-2*16+2*(ampl-1),N,N)
    assert gate(60-2*(16-2)+ampl-1,N,64)
    assert gate(15+ampl-2,N,16) and gate(14+ampl-2,N,16) and gate(13+ampl,N,16)
    assert gate(3+ampl,N,4) and gate(ampl,N,1)
    new_dg=60*1024*8;new_sparse=len(extra)*512
    assert new_dg==491520 and new_sparse==109056
    r=prior+new_dg+new_sparse
    print(json.dumps(dict(status='PASS',dimension=62,lower_bound=r,spikes=new_dg,added_banks=len(extra),sparse_gain=new_sparse,old_and_new_guards=len(E))))
    return r
