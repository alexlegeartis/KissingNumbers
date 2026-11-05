"""Rebuild a 14D E7-frame code and exactly check a restricted exchange bound.

Run from any directory. The search domain is explicitly restricted in README.
"""
from pathlib import Path
import collections,itertools,time,sys
import numpy as np
from scipy.sparse.csgraph import connected_components
p=Path(__file__).resolve().parent
bases=np.fromfile(p/'seed_e7_bases.bin',np.int8).reshape(9,7,7).astype(np.int16)
R=bases.reshape(-1,7);assert np.all(np.sum(R*R,axis=1)==4)
dot=R@R.T;orth=(dot==0);cliques=[]
def rec(chosen,allowed):
 if len(chosen)==7:cliques.append(tuple(chosen));return
 while allowed:
  if len(chosen)+len(allowed)<7:return
  u=allowed.pop(0);rec(chosen+[u],[v for v in allowed if orth[u,v]])
rec([],list(range(63)));assert len(cliques)==135
assignment=np.fromfile(p/'seed_assignment.bin',np.uint8);assert len(assignment)==63
rows=[]
for i,u in enumerate(R):
 for v in R[list(cliques[int(assignment[i])])]:
  for s,t in itertools.product((-1,1),repeat=2):rows.append(np.r_[s*u,t*v])
for half in range(2):
 for i,j in itertools.combinations(range(7),2):
  for s,t in itertools.product((-2,2),repeat=2):
   v=np.zeros(14,np.int16);v[7*half+i]=s;v[7*half+j]=t;rows.append(v)
X=np.asarray(rows,np.int16);assert X.shape==(1932,14) and len(set(map(tuple,X)))==1932 and np.all((X*X).sum(1)==8)
G=X@X.T;np.fill_diagonal(G,-9);assert G.max()<=4
old=X[np.array([v[np.flatnonzero(v)[0]]>0 for v in X])];assert len(old)==966
existing={tuple(v) for v in old};Q=[];B=[];hist=collections.Counter();start=time.time()
def consume(C):
 D=C@old.T;bad=abs(D)>4;counts=bad.sum(axis=1)
 for i,v in enumerate(C):
  if tuple(v) in existing:continue
  n=int(counts[i]);hist[n]+=1
  if n<=14:Q.append(v.copy());B.append(np.flatnonzero(bad[i]))
for i,j in itertools.combinations(range(14),2):
 C=np.zeros((2,14),np.int16);C[:,i]=2;C[:,j]=[2,-2];consume(C)
for i in range(14):
 for S in itertools.combinations([j for j in range(14) if j!=i],4):
  pos=sorted((i,*S));C=np.zeros((16,14),np.int16)
  for k,bits in enumerate(itertools.product((-1,1),repeat=4)):
   for j,s in zip(pos,(1,*bits)):C[k,j]=s*(2 if j==i else 1)
  consume(C)
patterns=np.array([[1,*s] for s in itertools.product((-1,1),repeat=7)],np.int16)
for S in itertools.combinations(range(14),8):
 C=np.zeros((128,14),np.int16);C[:,S]=patterns;consume(C)
assert not any(hist[k] for k in range(4)) and hist[4]==672 and hist[8]==7056
assert not any(hist[k] for k in (5,6,7,9,10,11,12,13))
Q=np.asarray(Q,np.int16);cost=np.array(list(map(len,B)));assert len(Q)==9392
sources={c:set(int(j) for i,b in enumerate(B) if cost[i]==c for j in b) for c in (4,8)}
assert len(sources[4])==84 and len(sources[8])==882 and not sources[4]&sources[8]
def has5(N,P,depth=0):
 if depth==5:return True
 if P.bit_count()<5-depth:return False
 while P:
  bit=P&-P;u=bit.bit_length()-1;P^=bit
  if has5(N,P&N[u],depth+1):return True
 return False
def no_five(G):
 D=~G.copy();np.fill_diagonal(D,False);n=len(D)
 masks=[sum(1<<int(j) for j in np.flatnonzero(D[i])) for i in range(n)]
 return not has5(masks,(1<<n)-1)
for c in (4,8):
 for oldline in sources[c]:
  ids=[i for i,b in enumerate(B) if cost[i]==c and oldline in b];assert len(ids)==(32 if c==4 else 64)
  P=Q[ids].copy();dots=P@old[oldline];assert np.all(abs(dots)>4);P[dots<0]*=-1
  conflict=(P@P.T)>4;np.fill_diagonal(conflict,False)
  if c==4:assert no_five(conflict)
  else:
   k,lab=connected_components(conflict);assert k==2 and sorted(np.bincount(lab))==[32,32]
   for component in (0,1):assert no_five(conflict[np.ix_(lab==component,lab==component)])
if '--export' in sys.argv:
 cache=p/'cache';cache.mkdir(exist_ok=True)
 blocks=np.full((len(Q),14),65535,np.uint16)
 for i,b in enumerate(B):blocks[i,:len(b)]=b
 np.savez_compressed(cache/'candidates14.npz',Q=Q.astype(np.int8),B=blocks,old=old.astype(np.int8))
 print('Exported',len(Q),'candidate lines with at most 14 seed-line conflicts')
print('PASS 1932-point E7-frame code, complete norm-8 shell, and pointwise local capacity bound; seconds',round(time.time()-start,1))
