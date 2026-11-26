from pathlib import Path
import os,sys,time,itertools,math
import numpy as np
from fractions import Fraction as F
import hashlib
src=Path(sys.argv[1]).resolve()
assert hashlib.sha256((src/'common.py').read_bytes()).hexdigest()=='2253a0f70c1e857c993d1f0ad4ff49b3e33e0598f2d6dfb033870a14aa402dc7'
sys.path.insert(0,str(src));from common import leech_shell
fn=src/'data/d29_certificate.npz'
if not fn.exists():fn=src/'v2.npz'
assert hashlib.sha256(fn.read_bytes()).hexdigest()=='b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824'
Z=np.load(fn);A=Z['added_integer'].astype(np.int32)
I=np.array([0,2,4,5,6,10,11,23]); rr=np.flatnonzero(np.sum(A[:,24:]**2,1)==48);assert len(rr)==1;ref=A[rr[0],24:];h0=A[rr[0],I]
shell=leech_shell().astype(np.int16);k=np.ones(len(shell),bool);k[Z['deleted_leech_indices']]=False;S=np.unique(shell[k][:,I],axis=0);S=S[(S*S).sum(1)>=52];del shell,k
A=np.delete(A,rr,axis=0); A=np.unique(np.c_[A[:,I],A[:,24:]],axis=0);T,inv=np.unique(A[:,8:],axis=0,return_inverse=True);inv=inv.reshape(-1);re=np.argsort(inv);B=A[re,:8].astype(np.int32); _,st=np.unique(inv[re],return_index=True);Ti=T.copy();T=T.astype(np.int32)
print('rare orbit',h0,ref,'shell projection',len(S),'aux',len(A),'tails',len(T),flush=True)
oldh=h0.copy(); sca=4; h0=h0*sca; Ti=Ti*sca;T=Ti.astype(np.int32)
# Weyl-reflected shapes on roots of the original E8 shell; each retained projection is checked, so signed-permutation enlargement is legitimate even outside the E8 orbit.
from collections import deque
def canon(a):
    v=tuple(sorted(map(abs,a),reverse=True)); return (v,0 if v[-1]==0 else sum(x<0 for x in a)%2)
q=deque([canon(h0), (canon(h0)[0],1-canon(h0)[1])]); seen=set(q);signs=[s for s in itertools.product([-1,1],repeat=8) if s.count(-1)%2==0]
while q:
    v,p=q.popleft();v=list(v);v[-1]*=1-2*p
    for s in signs:
        val=sum(x*y for x,y in zip(v,s));assert val%4==0
        a=canon([x-val//4*y for x,y in zip(v,s)])
        if a not in seen:seen.add(a);q.append(a)
sh=set(v for v,p in seen);per=set(t for v in sh for t in itertools.permutations(v));print('W classes',len(seen),'magnitudes',len(sh),flush=True)
print('unsigned perms',len(per),flush=True)
heads=[]
for v in per:
    v=np.array(v);nz=np.flatnonzero(v); X=np.tile(v,(1<<len(nz),1));X[:,nz]*=(1-2*((np.arange(len(X))[:,None]>>np.arange(len(nz))[None,:])&1)).astype(np.int16)
    ok=np.all(X.astype(np.int32)@S.T<=256,axis=1);heads.extend(X[ok])
print('functional candidate heads',len(heads),'vs',sum(2**sum(x!=0 for x in v) for v in per),flush=True); checks=0;maxraw=-1;possible=0;branches=0;empty=0
# Build signed slab directions and all actual orthogonal frames; each bound below is a literal weighted sum of squared slab intervals.
import math
dirs={}; mt=1
for i,t in enumerate(Ti):
    if not t.any():continue
    q=int(np.gcd.reduce(abs(t)));mt=math.lcm(mt,q);u=tuple(map(int,t//q));dirs.setdefault(u,[]).append((i,q))
base=[u for u in dirs if tuple(-x for x in u) in dirs and next(x for x in u if x)>0]
orth=np.array(base,int);lens=np.sum(orth*orth,1);unitden=math.lcm(*map(int,lens))
frames=[]
def complete(row,cand):
    if len(row)==5:frames.append(row);return
    for k,j in enumerate(cand):
        nxt=[a for a in cand[k+1:] if orth[a]@orth[j]==0]
        if len(nxt)>=4-len(row):complete(row+[j],nxt)
complete([],list(range(len(base))))
frames=np.array(frames,dtype=np.int32).reshape(-1,5); wgt=unitden//lens; sc=mt*mt*unitden
print('integer slab frames',len(base),len(frames), 'scale',sc,flush=True)
tabpos=[(np.array([i for i,q in dirs[u]]),np.array([mt//q for i,q in dirs[u]])) for u in base]
tabneg=[(np.array([i for i,q in dirs[tuple(-x for x in u)]]),np.array([mt//q for i,q in dirs[tuple(-x for x in u)]])) for u in base]
fcounts=[0,0,0,0]; openheads=[]
for h in heads:
    if np.array_equal(h,h0):continue
    cap=256-np.maximum.reduceat(B@h,st);checks+=1
    up=np.array([min(cap[i]*q) for i,q in tabpos]);lo=-np.array([min(cap[i]*q) for i,q in tabneg])
    if np.any(lo>up): fcounts[0]+=1;continue
    if len(frames):
        U=wgt*np.maximum(lo*lo,up*up);L=wgt*np.maximum(np.maximum(lo,-up),0)**2
        if U[frames].sum(1).min()<48*sc:fcounts[1]+=1;continue
        if L[frames].sum(1).max()>48*sc:fcounts[2]+=1;continue
    fcounts[3]+=1;openheads.append(h.copy())
    raise AssertionError('new head outside the claimed certificate')
assert sum(2**sum(x!=0 for x in v) for v in per)==401184
assert len(heads)==3528 and len(per)==3570 and fcounts==[3435,92,0,0]
# On the original head use the SAME bank (rare removed) to check the triangular prism.
h=h0;cap=256-np.maximum.reduceat(B@h,st);from fractions import Fraction as F
D={}
for t,cut in zip(Ti,cap):
    if not t.any():assert cut>=0;continue
    g=int(np.gcd.reduce(abs(t))); key=tuple(map(int,t//g));D[key]=min(D.get(key,F(10000)),F(int(cut),g))
for i in (1,4):
    for sign in (-1,1):
        k=[0]*5;k[i]=sign;assert D[tuple(k)]==0
p={}
for t,cut in D.items():
    a=np.array([t[j] for j in (0,2,3)]);g=int(np.gcd.reduce(abs(a)))
    if not g:assert cut>=0;continue
    k=tuple(map(int,a//g));p[k]=min(p.get(k,F(10000)),cut/g)
for k,v in {(-1,0,0):4,(0,1,0):4,(1,-1,0):-4,(0,0,1):4,(0,0,-1):4}.items():assert p[k]==v
points=np.array([[x,0,y,z,0] for x,y in [(-4,0),(-4,4),(0,4)] for z in [-4,4]])
assert np.all(points@Ti.T<=cap)
nrm=np.sum(points*points,axis=1);assert max(nrm)==48 and sum(nrm==48)==2
diam=max(int((u-v)@(u-v)) for u in points for v in points);assert diam==96
print('PASS: all new heads excluded; prism has precisely two alternative norm-48 tails and permits at most one donor.')
