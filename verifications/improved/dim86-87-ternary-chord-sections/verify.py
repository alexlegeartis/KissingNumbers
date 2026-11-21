#!/usr/bin/env python3
"""Exact local replay of the 86D chord and sparse union and its 87D extension."""
from fractions import Fraction as Q
from itertools import combinations,product
from pathlib import Path
import hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent
V=HERE.parent
QR=V/'dim68-91-lattice-code-continuation/ternary_qr48_20261117/certificate.json'
S=V/'dim90-92-hexagonal-phase-layers/hexagon_sparse15_phase7812_20261117'
OUTER=S/'outer_42_12_14_25_41_42_44.npz'
PHASE=S/'phase7812.txt'
SOURCES={
 QR:'1505ab829025e938ad09de3575656ccf5f547537d8fd5432533596fa963b22e8',
 OUTER:'e8517c99f17d582d9753de332a9245f75edf2ded261753a2d9237e1550dc98d7',
 PHASE:'d5ba3940239da9750271082790a2ac3f62aa52a14c4ee66727f9a66121f593e2',
}
def rank3(A):
 A=np.array(A,dtype='i8')%3;rank=0
 for i in range(A.shape[1]):
  v=np.flatnonzero(A[rank:,i])
  if not len(v):continue
  a=rank+v[0]; A[[rank,a]]=A[[a,rank]]
  A[rank]=A[rank]*A[rank,i]%3
  for j in range(rank+1,len(A)):
   A[j]=(A[j]-A[j,i]*A[rank])%3
  rank+=1
 return rank

def min_dist(W,d):
 W=np.asarray(W)
 assert W.ndim==2 and 0<W.shape[1]<=64 and W.min()>=0 and W.max()<=2
 rot=np.ones(W.shape[1],dtype='u8')<<np.arange(W.shape[1],dtype='u8')
 A=(W==1)@rot;B=(W==2)@rot
 assert len(set(zip(A,B)))==len(W)
 best=W.shape[1]
 pop=getattr(np,'bitwise_count',None)
 def fallback(v):
  v=v-((v>>np.uint64(1))&np.uint64(0x5555555555555555))
  v=(v&np.uint64(0x3333333333333333))+((v>>np.uint64(2))&np.uint64(0x3333333333333333))
  v=(v+(v>>np.uint64(4)))&np.uint64(0x0f0f0f0f0f0f0f0f)
  return (v*np.uint64(0x0101010101010101))>>np.uint64(56)
 if pop is None:pop=fallback
 for i in range(0,len(W),128):
  j=min(i+128,len(W))
  x=pop((A[i:j,None]^A[None,:j])|(B[i:j,None]^B[None,:j]))
  x=x[np.arange(i,j)[:,None]>np.arange(j)[None,:]]
  if x.size:best=min(best,int(x.min()))
 assert best==d
 return len(W)*(len(W)-1)//2

def results():
 assert __debug__,'Run without -O'
 for f,s in SOURCES.items():assert hashlib.sha256(f.read_bytes()).hexdigest()==s
 G=np.array(json.loads(QR.read_text())['generator_Q48'],dtype='i8')
 assert G.shape==(24,48) and G.min()==0 and G.max()==2
 assert rank3(G)==24 and (G@G.T%3==0).all()
 # First symbol fixed at 0; next eight restricted to {1,2}.
 assert rank3(G[:,:9])==9
 nbody=39; heads=8; norm=Q(nbody)+heads*Q(3,4)
 assert 2*nbody+heads==86 and norm==45
 dense=(2**8)*3**(24-9)
 # Parent minimum 15 is the published and previously pinned QR48 theorem.
 # Self-orthogonality gives wt divisible by 3. The fixed coordinate is zero in a difference.
 assert all(-norm/2<=norm-Q(3,2)*d<=norm/2 for d in range(15,48,3))
 # Full triangle entries do not contain their own negatives.
 tri=[(2,0),(-1,1),(-1,-1)] # exact metric diag(1,3)/4
 assert not set(tri)&{(-a,-b)for a,b in tri}
 assert {Q((a*c+3*b*d),4)for a,b in tri for c,d in tri}=={Q(-1,2),Q(1)}
 center=tuple(Q(a+b,2)for a,b in zip(tri[1],tri[2]))
 chord=[tuple(Q(a)-b for a,b in zip(p,center))for p in tri[1:]]
 assert chord==[(0,1),(0,-1)] and {Q(a*a+3*b*b,4)for a,b in chord}=={Q(3,4)}
 # Center subtraction is constant in each shortened plane and leaves all pair distances unchanged.
 W=np.load(OUTER,allow_pickle=False)['words']
 zeros=[7,9,35]
 assert W.shape==(37468,42) and W.dtype==np.uint8
 W=W[np.all(W[:,zeros]==0,axis=1)]; W=np.delete(W,zeros,axis=1)
 assert W.shape==(9662,39) and np.all(np.count_nonzero(W,axis=1)==15)
 pairs_W=min_dist(W,15)
 L=PHASE.read_bytes().splitlines()
 assert len(L)==7812 and all(len(w)==15 for w in L)
 Z=np.array([[b-48 for b in w]for w in L],dtype='i4')
 pairs_Z=min_dist(Z,5)
 roots=[tuple(a-b for a,b in zip(tri[i],tri[j]))for i,j in [(0,2),(1,0),(2,1)]]
 def bdot(p,q):return Q(p[0]*q[0]+3*p[1]*q[1],4)
 assert [bdot(h,h)for h in roots]==[3,3,3]
 assert all(bdot(x,y)==Q(-3,2)for x,y in combinations(roots,2))
 assert max(abs(bdot(h,v))for h in roots for v in tri)==Q(3,2)
 assert 15*3==norm and 15*Q(3,2)==norm/2
 assert min(2*3+2*bdot(x,y)for x in roots for y in roots)==3
 assert min(2*3-2*bdot(x,y)for x,y in combinations(roots,2))==9
 sparse=len(W)*len(Z)
 # Independent doubled-coordinate expansion of all 240 pure-head E8 roots.
 A=[]
 for a,b in combinations(range(8),2):
  for p,q in product((-2,2),repeat=2):
   x=[0]*8;x[a]=p;x[b]=q;A.append(tuple(x))
 A += [s for s in product((-1,1),repeat=8)if sum(v<0 for v in s)%2==0]
 dot=lambda x,y:sum(a*b for a,b in zip(x,y))
 assert len(A)==len(set(A))==240 and all(dot(x,x)==8 for x in A)
 assert max(dot(x,y)for x,y in combinations(A,2))==4
 assert norm*heads*Q(3,4)<(norm/2)**2 # all dense / pure-head comparisons by Cauchy-Schwarz
 n=2*dense+sparse+len(A)
 assert n==7422120168 and n+2==7422120170
 prior={86:6228511388,87:7266268830}
 return {str(d):dict(total=n+(2 if d==87 else 0),
         previous=prior[d],gain=n+(2 if d==87 else 0)-prior[d],dense=2*dense,sparse=sparse,root_tail=240,
         orthogonal_poles=2 if d==87 else 0,
         outer_words=len(W),phase_words=len(Z),outer_pairs=pairs_W,phase_pairs=pairs_Z) for d in (86,87)}
if __name__=='__main__':
 print(json.dumps(results(),indent=2,sort_keys=True))
