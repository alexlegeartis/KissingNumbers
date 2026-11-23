#!/usr/bin/env python3
"""Exact signed-overlap and DG replacement replay in 64 and 65 dimensions."""
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from subprocess import run, check_output
from tempfile import TemporaryDirectory
import gzip
import json
from dg import code

assert __debug__
P = Path(__file__).resolve().parent
R = P.parent

def checked(path, digest, gz=False):
    data = path.read_bytes()
    if gz:
        data = gzip.decompress(data)
    assert sha256(data).hexdigest() == digest, path
    return data

signed = R/'dim63-67-signed-overlap'
atlas = R/'dim63-68-affine-hyperplane-heads'
quad65 = R/'dim65-optimal-quadruple-layer'
cubic = checked(signed/'cubic80.txt', '5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f')
w0 = checked(signed/'signed_witness.txt.gz', '5173fabd124bd23e077d00c956534e7f56543101b28903cbac017ed252c34be7', True)
w1 = checked(atlas/'joint_witness.txt', 'ff1375eff3010defcdd57145c402d17664adac39eea8d1fcddc084e08d67b87e')
p64 = json.loads(checked(atlas/'claims.json', '1d6e67a0b321db8ec27db97efb64a62fcbc20dee1674bf2836c3defd45ed7e9a'))
p65 = json.loads(checked(quad65/'claims.json', 'a6ab64b9f793c54568be50a995ced72e58754698d41e2375728c7c8b9733c160'))
J = json.loads(checked(P/'phase64screen.json', '6f414a7994352ed0beafad182d7e664f293fed94dfb87ec468ab00464279094e'))
assert J['nondual'] == []
X = int(J['bits'],16)

def mon(u): return sum(1<<i for i in range(64) if (i&u)==u)
r2 = [mon(u) for u in range(64) if u.bit_count()<=2]
r3 = [mon(u) for u in range(64) if u.bit_count()==3]
cube = []
for c in map(int, cubic.split()):
    u = 0
    for i,r in enumerate(r3):
        if c>>i&1: u ^= r
    cube.append(u)
assert len(set(cube)) == 80
def join(a,b):
    assert -(1<<31)<=a<1<<31 and -(1<<31)<=b<1<<31
    return (a&0xffffffff) | ((b&0xffffffff)<<32)
def dep(i,m):
    p=0
    while m:
        b=m&-m
        if i&1: p|=b
        i>>=1
        m^=b
    assert not i
    return p

def run_helper(t, name, input_file, capture=False):
    z=t/name
    run(['g++','-O3','-std=c++20','-march=native',str(P/(name+'.cpp')),'-o',str(z)],check=True)
    if capture: return check_output([str(z),str(input_file)],text=True)
    run([str(z),str(input_file)],check=True)

with TemporaryDirectory() as tmp:
    T=Path(tmp)
    (T/'cubic.txt').write_bytes(cubic)
    B=dict(enumerate(map(int,run_helper(T,'verify_body',T/'cubic.txt',True).split())))
    E=[]
    nxt=60524
    def guard(i,j,a,b):
        assert i>j and i in B and j in B
        assert 0<=a<1<<16 and 0<=b<1<<17
        g=dep(a,B[i])
        assert g.bit_count()==8 and g==dep(b&65535,B[j])
        E.append((i,j,g))
    for l in w0.splitlines():
        i,j,a,b=map(int,l.split())
        if j==-1:
            assert i==nxt and i not in B
            nxt+=1
            B[i]=join(a,b)
        else:
            assert i==nxt-1
            if j>=0: guard(i,j,a,b)
            else: assert j==-2
    assert nxt==62878
    n0=len(B)
    nxt=67859
    for l in w1.splitlines():
        i,j,t,a,b=map(int,l.split())
        if j==-1:
            assert i==nxt and i not in B and 0<=t<16
            nxt+=1
            if t==0: B[i]=join(a,b)
        else:
            assert i==nxt-1
            if j>=0 and i in B and j in B: guard(i,j,a,b)
    assert nxt==68051 and n0==40350
    assert (len(B),len(E))==(J['masks'],J['guards'])==(40414,14723) and X>>len(E)==0
    H=defaultdict(list)
    for k,(i,j,g) in enumerate(E):
        b=X>>k&1
        H[i].append((g,b));H[j].append((g,b^1))
    for i,s in B.items():
        assert s.bit_count()==16 and all(not ((r&s).bit_count()&1) for r in r2+cube)
        F={0:0}
        for g,p in [(s,0)]+H[i]:
            if g in F: assert F[g]==p
            else: F.update({x^g:y^p for x,y in list(F.items())})
        assert len(F) in (2,4,8,16,32)
        assert all(g in (0,s) or g.bit_count()==8 for g in F)
    inds={k:i for i,k in enumerate(B)}
    with (T/'pairs.txt').open('w') as f:
        print(len(B),len(E),file=f)
        print(*B.values(),sep='\n',file=f)
        for i,j,_ in E: print(inds[i],inds[j],file=f)
    run_helper(T,'check_pairs',T/'pairs.txt')

C,A=code()
assert len(C)==1024 and len(A)==128
coords=[(i,i+1) for i in range(0,64,2)]
assert sorted(i for s in coords for i in s)==list(range(64))
assert len(cube[:len(coords)])==len(set(cube[:len(coords)]))==32
# Each DG word and cubic translate is an actual RM80 parent word.
qparts={0}
for u in range(64):
    if u.bit_count()==2:
        v=mon(u)
        qparts|={w^v for w in qparts}
assert len(qparts)==2**15 and set(C)<=qparts
a1={0}
for u in range(64):
    if u.bit_count()<=1:
        v=mon(u)
        a1|={w^v for w in a1}
assert set(A)==a1 and len(a1)==2**7
assert len(set(C))==len(C) and all(w==0 or w.bit_count()>=32 for w in A)

a,b=Q(28,5),Q(29,20)
assert a>b>1
N=62+a*a+b*b
assert N==Q(38185,400)
def gate(d,m=N,k=N):
    assert d<=0 or 4*d*d<=m*k,(d,m,k)
# For a single owner's two swapped profiles, two opposite weighted products.
gate(62-2*a*b)
# Same pair and distinct owners: same profile or swapped profiles.
gate(16+a*a+b*b-2)
gate(16+2*(a*b+1))
# Different pairs, disjoint modified sites and distinct parent cubic cosets.
gate(32+2*(a+b))
# The literal owner's old dense word is removed.
gate(32+a+b,N,64)
# Sparse intersects both spikes, large only, small only, or neither.
# When the large site is present, odd donor parity versus even old
# sparse signs ensures a mismatch at a site of magnitude at least 1.
for dot in (12+a+b,13+a,15+b,16): gate(dot,N,16)
# Any support of four distinct sites, including the optimal 65D stock.
gate(2+a+b,N,4)
gate(0,N,1)  # the two retained 65D tail poles

D=len(cube)*2**22
S=len(B)*2048
F=len({frozenset((u,v,w,u^v^w)) for u,v,w in combinations(range(64),3)})
assert F==10416
old64=D+S+16*F+2*64
old65_atlas=126*2*2048
F65=(65*(64*(63//2)//3))//4
assert F65==10741
old65=D+S+16*F65+2*65+old65_atlas
assert old64==p64['lower_bounds']['64']==418478976
assert old65==p65['lower_bounds']['65']==419000274
owners=len(coords)*len(C)*len(A)
assert owners==4194304
n64=old64-128+owners
n65=old65-128-old65_atlas+owners
out={'status':'PASS','lower_bounds':{'64':n64,'65':n65}}
assert out==json.loads((P/'claims.json').read_text())
print(json.dumps(out,sort_keys=True))
