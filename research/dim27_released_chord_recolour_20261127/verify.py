#!/usr/bin/env python3
"""Verify the 27D fixed-chord and discrete-recolouring obstructions.
Python 3, NumPy; give upstream d27_certificate.npz as sole argument.
Binary32 dot-product accumulation is exact here (< 2**24 in absolute sum).
"""
from collections import Counter
import hashlib
import sys
import numpy as np

PIN = '8b704bb49503e3a4f6ad6b567690a7a2a9abfb30600f453313eb0510e6d29a94'
OCTAD = [0,2,4,5,6,10,11,23]
J = [j for j in range(24) if j not in OCTAD]

def leech():
    g=np.zeros((12,24),dtype=np.int16)
    for i in range(12):
        for j in (0,2,4,5,6,10,11):
            g[i,(i+j)%23]=1
    g[:,23]=g[:,:23].sum(1)%2
    c=((np.arange(4096)[:,None]>>np.arange(12)&1)@g)%2
    assert len(np.unique(c,axis=0))==4096
    bank=[]
    for i in range(24):
        for j in range(i):
            for a,b in ((8,8),(8,-8),(-8,8),(-8,-8)):
                r=np.zeros(24,dtype=np.int16);r[[i,j]]=[a,b];bank.append(r)
    signs=((np.arange(256)[:,None]>>np.arange(8)&1).sum(1)%2==0)
    sign=4*(1-2*((np.arange(256)[:,None]>>np.arange(8)&1)[signs]))
    for w in c[c.sum(1)==8]:
        x=np.zeros((128,24),np.int16);x[:,w==1]=sign;bank.extend(x)
    for w in c:
        x=np.tile(4*w-2,(24,1)).astype(np.int16);x[np.arange(24),np.arange(24)]*=-3
        bank.extend(x)
    x=np.unique(np.array(bank,dtype=np.int32),axis=0)
    assert x.shape==(196560,24) and (np.sum(x*x,axis=1)==128).all()
    return x

path=sys.argv[1]
assert hashlib.sha256(open(path,'rb').read()).hexdigest()==PIN
D=np.load(path,allow_pickle=False); A=D['added_integer'].astype(np.int32)
L=leech(); d=D['deleted_leech_indices']; X=L[d]
assert len(np.unique(d))==13248 and np.max(abs(L))<=8
index=np.full(len(L),-1,np.int32);index[d]=np.arange(len(d))
E,first,multiplicity=np.unique(np.sort(D['midpoint_endpoints'],axis=1),axis=0,return_index=True,return_counts=True)
assert (multiplicity==2).all() and len(E)==len(X)
assert (index[E]>=0).all(); E=index[E]
C=(X[E[:,0]]+X[E[:,1]])//2
M=D['midpoint_rows']; TM=A[M[first],24:]
assert (np.sum(X[E[:,0]]*X[E[:,1]],axis=1)==64).all()
assert len(np.unique(C,axis=0))==len(E)
assert set(map(tuple,A[M]))==set(map(tuple,np.vstack((np.c_[C,TM],np.c_[C,-TM]))))
assert len(np.unique(A[M],axis=0))==2*len(E)
assert (np.sort(abs(TM),axis=1)==[0,4,4]).all()
assert np.max(abs(A))<=8
adj=[[] for _ in X]; ae=[[] for _ in X]
for e,(i,j) in enumerate(E):
    adj[i].append(j); adj[j].append(i); ae[i].append(e); ae[j].append(e)
assert all(len(t)==2 for t in adj)
adj=np.array(adj);ae=np.array(ae)
x=X.astype(np.float32); hist=Counter(); mixed=0; survivors=[]
for i in range(0,len(X),256):
    r,s=np.where(x[i:i+256]@x[i:].T==64)
    r+=i;s+=i
    k=r<s;r,s=r[k],s[k]
    k=~(adj[r]==s[:,None]).any(1);r,s=r[k],s[k]
    p=np.einsum('ijk,ik->ij',X[adj[r]],X[s]);q=np.einsum('ijk,ik->ij',X[adj[s]],X[r])
    u,v=p.max(1),q.max(1)
    is_mixed=((X[r][:,OCTAD]!=X[s][:,OCTAD]).any(1)&(X[r][:,J]!=X[s][:,J]).any(1))
    mixed+=int(is_mixed.sum())
    hist.update(zip(np.minimum(u,v).tolist(),np.maximum(u,v).tolist()))
    k=(64-u)*(64-v)>=4096
    assert (u[k]==0).all() and (v[k]==0).all() and not is_mixed[k].any()
    if k.any():survivors.append(np.stack((r[k],s[k]),axis=1))
assert hist=={(0,0):186726,(0,32):41472,(0,64):37024,(32,32):2206720,
              (32,64):468992,(64,64):12384}
assert mixed==2756992
edges=np.vstack(survivors);B=(X[edges[:,0]]+X[edges[:,1]])//2
B=B.astype(np.float32); cf=np.ascontiguousarray(C.T,dtype=np.float32)
bound=[]
for k in range(3):
    z=np.eye(3,dtype=np.int32)[k]*8
    f=A[(A[:,24:]==z).all(1),:24];g=A[(A[:,24:]==-z).all(1),:24]
    assert len(f)==48 and set(map(tuple,f))==set(map(tuple,g))
    assert np.all(f@f.T<=64) and set(map(tuple,f))==set(map(tuple,-f))
    bound.append(64-(B@f.astype(np.float32).T).max(1).astype(np.int32))
b=np.stack(bound,axis=1)
assert (b>=0).all() and (b<=32).all()
clear=np.sum(b*b,axis=1)>=2048 # squared radius 32 with units t*8
assert int(clear.sum())==2406
# Two opposite parity vertices in a saturated-contact component force any
# six-line recolouring of them into the uniquely orthogonal pair of lines.
par=np.arange(len(C)); flip=np.zeros(len(C),np.int8)
def root(i):
    p=0;j=i
    while par[j]!=j:
        p^=int(flip[j]);j=par[j]
    return j,p
cfrows=C.astype(np.float32)
for i in range(0,len(C),192):
    r,s=np.where(cfrows[i:i+192]@cf==64);r+=i
    for e,f in zip(r[r<s],s[r<s]):
        assert int(TM[e]@TM[f])==0
        re,pe=root(e);rf,pf=root(f)
        if re!=rf:
            par[max(re,rf)]=min(re,rf);flip[max(re,rf)]=pe^pf^1
        else:
            assert pe^pf==1
RP=np.array([root(i) for i in range(len(C))])
for z in np.flatnonzero(clear):
    ix=np.flatnonzero(B[z]@cf==64)
    cb=RP[ix]
    assert any(0 in cb[cb[:,0]==r,1] and 1 in cb[cb[:,0]==r,1]
               for r in np.unique(cb[:,0]))
print('PASS: %d new released chords, %d mixed; all closed under arbitrary new real tails and every six-line recolouring.'
      %(sum(hist.values()),mixed))
