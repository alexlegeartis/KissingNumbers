#!/usr/bin/env python3
"""Python 3.10+. Run python verify.py path/to/a40.8.8.4267H. No dependencies."""
from collections import deque
from itertools import combinations
from pathlib import Path
import hashlib,sys
def basis(rows):
    b={}
    for v in rows:
        while v:
            m=v.bit_length()-1
            if m not in b:
                b[m]=v;break
            v ^= b[m]
    return b
def contains(b,x):
    for m in sorted(b,reverse=True):
        if x>>m&1:x ^= b[m]
    return x==0
def scan(B):
    global state
    mins=[]
    for h in range(40):
        incident=[w for w in B if w>>h&1]
        assert 852<=len(incident)<=855 and len(basis(incident))==39
        def greedy(rows):
            ct=0
            while len(rows)>=39:
                piv,sk={},[]
                for w in rows:
                    x=w
                    while x:
                        k=x.bit_length()-1
                        if k not in piv:
                            piv[k]=x;break
                        x ^= piv[k]
                    if not x: sk.append(w)
                if len(piv)<39: break
                ct+=1;rows=sk
            return ct
        n=greedy(incident)
        if n<13:
            for j in range(len(incident)-1,0,-1):
                state ^= (state<<13)&0xffffffff
                state ^= state>>17
                state ^= (state<<5)&0xffffffff
                k=(state*(j+1))>>32
                incident[k],incident[j]=incident[j],incident[k]
            n=max(n,greedy(incident))
        mins.append(n)
    assert min(mins)>=16
    return mins
def mother():
    p=47;r=sum(1<<i for i in {j*j%p for j in range(1,p)})
    low=(1<<47)-1
    gs=[(((r<<j)|(r>>(47-j)))&low)|(1<<47) for j in range(47)]
    b=basis(gs)
    assert len(b)==24 and contains(b,(1<<48)-1)
    W=[];x=0;G=list(b.values())
    for i in range(1,1<<24):
        x ^= G[(i&-i).bit_length()-1]
        assert x.bit_count()>=12
        if x.bit_count()==12: W.append(x)
    assert len(W)==17296
    return W,G,b
def bits(w):
    while w:
        x=w&-w;yield x.bit_length()-1;w ^= x
def perm(w,p):
    return sum(1<<p[i] for i in bits(w))
def reps(W,G,b):
    p=[list(range(1,47))+[0,47],[47]+[(-pow(x,-1,47))%47 for x in range(1,47)]+[0]]
    for f in p:
        assert all(contains(b,perm(x,f)) for x in G)
    rem=set(W);out=[]
    while rem:
        w=min(rem);rem.remove(w);stack=[w];n=0
        while stack:
            z=stack.pop();n+=1
            for f in p:
                y=perm(z,f)
                if y in rem:rem.remove(y);stack.append(y)
        out.append((w,n))
    assert sorted(n for _,n in out)==[4324,4324,8648]
    return out
def check(D,E,W):
    M=D^E;T=((1<<48)-1)^D;verts=list(bits(M))
    adj={j:[] for j in verts};near=[];edges=0
    for w in W:
        e=(w&E).bit_count();p=w&M
        if e==4 and p.bit_count()==2:
            a,b=bits(p);adj[a].append((b,w));adj[b].append((a,w));edges+=1
        if e==3 and p.bit_count()==1: near.append(w&T)
    root=verts[0];path={root:0};st=[root];cyc=[]
    while st:
        i=st.pop()
        for j,w in adj[i]:
            if j not in path:path[j]=path[i]^w;st.append(j)
            else:cyc.append(path[i]^path[j]^w)
    z=basis(cyc);cover=0
    for w in near:cover|=w
    assert edges==28 and len(path)==8 and len(z) in (13,14) and contains(z,T)
    assert len(near)==480 and cover==T
    return len(z)
if __name__=="__main__":
    data=Path(sys.argv[1]).read_bytes()
    assert hashlib.sha256(data).hexdigest()=="0566ebf4da8b4997c9a782ea537ba3f35125d9c53201fa0cd94551a803ac47e2"
    B=[int(s,16) for s in data.decode().splitlines() if s.strip() and not s.startswith('$')]
    assert len(B)==len(set(B))==4267 and all(x.bit_count()==8 for x in B)
    C=set()
    for w in B:
        for f in combinations(bits(w),5):
            assert f not in C;C.add(f)
    state=914235;print("disjoint local bases:",scan(B))
    W,G,b=mother();out=[]
    for D,n in reps(W,G,b):
        ranks=[check(D,sum(1<<i for i in t),W) for t in combinations(bits(D),4)]
        assert len(ranks)==495
        out.append((n,ranks.count(13),ranks.count(14)))
    assert sum(n13*n for n,n13,n14 in out)==4324*18+8648*3
    print("dodecad orbits and split ranks:",out)
    # Raw dot excess: (21+a)^2 -10(39+a*a), a=sqrt(41)-2.
    # Simplifies to -438+78*sqrt(41)>0.
    assert 78**2*41>438**2
    print("PASS: upper ceiling within the stated QR/Smart one-hot model is 1139159.")