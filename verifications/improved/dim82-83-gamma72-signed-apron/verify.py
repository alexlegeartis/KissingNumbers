#!/usr/bin/env python3
"""Exact replacement of the norm-4+8 packet component in signed Γ72 tails."""
from pathlib import Path
import numpy as np
import json,runpy,sys
from hashlib import sha256
from functools import lru_cache
assert __debug__
P=Path(__file__).resolve().parent
BASE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent
old=BASE/'dim82-83-gamma72-signed-packets'
PIN={'verify.py': 'ce315771a9abdffe250b30beae646d17b0353f9a0ca77e17bafe489f6f4c6ce5', 'claims.json': '8ef197d96fb09b4f5e0dafa1c141e04ba7746bb79331739ba77dafe4ae5a991f', 'refill.py': '1ae9ee8adff4f4fb23ae9a03e8b35af3cf1be9b85dab1c101274c60f9641cb2a', 'tail82.npz': 'bc4d17399173208c664700e5a97ff6e8b04ffe54b72b5ce20a54c951e9e8015c', 'tail83.txt': '8409b9a22f93422bde6a245d7fd628d4548933e21c0ba375dfa9013e6a5d0813', 'witness82.npz': '8a10964cf804477b979338a8579654d2a88b4f539fac8292ea0be00558895961', 'witness83.npz': '5997cf14725c85a72e084cfab500304e350d67afce364232c62245fb172dde8b'}
assert all(sha256((old/name).read_bytes()).hexdigest()==digest for name,digest in PIN.items())
sys.path.insert(0,str(old))
sys.argv=[str(old/'verify.py'),str(BASE),'82','83']
R=runpy.run_path(sys.argv[0])
M,tail=R['M'],R['tail']
from refill import pack24
@lru_cache(9000)
def signs(i,j):
 D=M.R8[i]@M.V8[j].T
 return pack24(D==4),pack24(D==-4)
@lru_cache(None)
def rel(a,b):
 d=M.delta[a,b];assert d
 a0,b0=M.match(1,a,b,d)
 a1,b1=M.match(1,a,b,np.uint64(0))
 assert np.all(M.U[a,a0,2]==M.U[b,b0,2])
 assert np.all(M.U[a,a1,:2]==M.U[b,b1,:2])
 tab=np.array([signs(int(i),int(j))for i,j in zip(M.U[a,a1,2],M.U[b,b1,2])],dtype='u4').reshape(-1,2,24)
 return a0,b0,a1,b1,tab
base=json.loads((old/'claims.json').read_text())
new=json.loads((P/'claims.json').read_text())
for n,c in new['lower_bounds'].items():
 nn,E=tail(int(n));Z=np.load(old/f'witness{n}.npz',allow_pickle=False)
 S=np.load(P/f'apron{n}.npz',allow_pickle=False)
 A=S['full48'];assert A.dtype==np.dtype('u4') and A.shape==(len(Z['row']),M.N[1],2)
 assert not (A>>24).any() and np.array_equal(S['row'],Z['row'])
 rows=np.array(Z['row'],dtype='i8')
 for j,(i,v,a) in enumerate(rows):
  m=A[j]
  for k,(h,w,b) in enumerate(rows[:j]):
   if w!=v or abs(i-h)==nn:continue
   if a==b:assert not (m&A[k]).any()
   elif E[i,h]:
    a0,b0,a1,b1,Sign=rel(a,b)
    assert not (m[a0]&A[k,b0]).any()
    if not len(a1):continue
    X0=Sign[:,0];X1=Sign[:,1]
    P0=A[k,b1,0,None];P1=A[k,b1,1,None]
    V0=pack24(((X0&P0)|(X1&P1))!=0)
    V1=pack24(((X1&P0)|(X0&P1))!=0)
    assert not (V0&m[a1,0]).any() and not (V1&m[a1,1]).any()
 before=sum(int(x).bit_count()for x in Z['full48'].ravel())
 after=sum(int(x).bit_count()for x in A.ravel())
 gain=2*after-4*before
 assert gain>0 and c==base['lower_bounds'][n]+gain and c>base['previous'][n]
 print('PASS signed 4+8 replacement',n,c,'gain over fixed signed',gain,flush=True)
