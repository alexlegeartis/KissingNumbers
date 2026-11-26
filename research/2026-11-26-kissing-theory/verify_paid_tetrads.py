#!/usr/bin/env python3
"""Replay the exact integer donor prices against Takhanov v0.2 d26_blocks.json.
Usage: python3 verify_paid_tetrads.py path/to/d26_blocks.json
"""
import hashlib,itertools,json,sys
from pathlib import Path
raw=Path(sys.argv[1]).read_bytes()
assert hashlib.sha256(raw).hexdigest()=="5018425664146ee2ace418021537b8b13863fc57c7dc64376ab6f9e811aa5162"
w=json.loads(raw)
B=w["five_coordinate_supports"];Q=w["four_coordinate_supports"]
T=w["template_numerators_and_denominators"]
def sg(s,i): return -1 if (s>>i)&1 else 1
def blocks(U,N):
    return [[s for s in range(16) if any(
        (v:=sum(u[i]*sg(s,k) for k,i in enumerate(q)))>0 and v*v>N
        for u in U)] for q in Q]
def price(U,N,b):
    # Precompute signed projections. For each orientation compute the union
    # of template rows hit by at least one donor, with integer comparisons.
    short=[tuple(u[k] for k in b) for u in U]
    best=27;hist={};receipt=None
    for p in itertools.permutations(range(5)):
        for s in range(32):
            lost=[]
            for j,t in enumerate(T):
                a=[t[k]*sg(s,i) for i,k in enumerate(p)]
                if any((v:=sum(x*y for x,y in zip(a,u)))>0 and 4*v*v>N*t[5]**2 for u in short):
                    lost.append(j)
            k=len(lost); hist[k]=hist.get(k,0)+1
            if k<best: best=k; receipt=(p,s,lost)
    return best,hist,receipt
def check(U,N):
    assert all(sum(k*k for k in u)==N and 4*max(map(abs,u))**2<=N for u in U)
    assert all((sum(abs(u[i]) for i in range(4))+2*abs(u[8]))**2<=4*N for u in U)
    for a,b in itertools.combinations(U,2):
        assert 2*sum(x*y for x,y in zip(a,b))<=N
def name(u):
    a=[0]*10
    for k,x in u.items(): a[k]=x
    return a
U=[name({0:25,1:21*sg(s,0),2:21*sg(s,1),3:21*sg(s,2),8:6,9:23}) for s in (0,3,5,6)]
check(U,2513)
assert [(i,z) for i,z in enumerate(blocks(U,2513)) if z]==[(5,[0,2])]
for k in range(1,5):
    P=[price(U[:k],2513,b) for b in B]
    assert [r[0] for r in P]==[[1,3,5,7][k-1],0,0,0]
    print("quartet",k,P[0])
A=[name({0:3*sg(s,0),1:3*sg(s,1),2:2*e,3:2*e,8:e,9:3*sg(s,2)}) for s in range(8) for e in (1,-1)]
check(A,36)
assert [(i,len(z)) for i,z in enumerate(blocks(A,36)) if z]==[(5,16),(9,16)]
for U in (A[:1],A):
    P=[price(U,36,b) for b in B]
    assert [r[0] for r in P]==[1 if len(U)==1 else 26,0,0,0]
    print("two-tube",len(U),P[0])
print("PASS: incumbent orientation prices only; no new kissing bound.")

# Conditional two-chamber +1 reduction, last existence step not supplied.
U1=[name({0:25,1:21*sg(s,0),2:21*sg(s,1),3:21*sg(s,2),8:6,9:23}) for s in (0,3,5,6)]
U4=[name({4:25,5:21*sg(s,0),6:21*sg(s,1),7:21*sg(s,2),8:6,9:23}) for s in (0,3,5,6)]
for X,head in ((U1,range(4)),(U4,range(4,8))):
    assert all(sum(v*v for v in x)==2513 and 4*max(map(abs,x))**2<=2513 for x in X)
    assert all((sum(abs(x[i]) for i in head)+2*abs(x[8]))**2<=4*2513 for x in X)
    assert all(2*sum(a*b for a,b in zip(x,y))<=2513 for x,y in itertools.combinations(X,2))
assert all(sum(a*b for a,b in zip(x,y))==565 and 2*565<2513 for x in U1 for y in U4)
assert [(i,z) for i,z in enumerate(blocks(U1+U4,2513)) if z]==[(5,[0,1,2])]
assert all(len(set(x)&set(y))<=2 for x,y in itertools.combinations(B,2))
assert all(len(set(q)&set(B[j]))<=2 for q in Q for j in (0,3))
for b,s in ((B[1],16),(B[2],4)):
    for x in U1+U4:
        for t in T:
            v=sum(x[j]*t[i]*sg(s,i) for i,j in enumerate(b))
            assert v<=0 or 4*v*v<=2513*t[5]**2
assert 200864+20+224-3+52+8+2*24==201213
print("PASS: two chambers reduce +1 to a common 24-row clipped code, not supplied.")

N=1000911
def balanced(k,s):
    return name({k:500,k+1:405*sg(s,0),k+2:405*sg(s,1),
                 k+3:405*sg(s,2),8:94,9:500})
W=[[balanced(k,s) for s in (0,3,5,6)] for k in (0,4)]
F=sum(W,[])
for X,head in zip(W,(range(4),range(4,8))):
    assert all(sum(v*v for v in x)==N and 4*max(map(abs,x))**2<=N for x in X)
    assert all((sum(abs(x[i]) for i in head)+2*abs(x[8]))**2<=4*N for x in X)
    assert all(2*sum(a*b for a,b in zip(x,y))<=N for x,y in itertools.combinations(X,2))
assert all(2*sum(a*b for a,b in zip(x,y))<=N for x in W[0] for y in W[1])
assert [(i,z) for i,z in enumerate(blocks(F,N)) if z]==[(5,[0,1,2])]
assert all(sum(abs(x[i]) for i in b)**2<=N for x in F for b in B[1:3])
print("PASS: balanced quartet has the same three paid rows; the 24-row witness is open.")
