"""Exact replay of the 82D and 83D extra norm-four layers."""
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import json,runpy,sys
import numpy as np
assert __debug__
P=Path(__file__).resolve().parent
base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent
old=base/'dim82-83-gamma72-signed-apron'
for n,h in json.loads((P/'pin.json').read_text()).items():
 assert sha256((old/n).read_bytes()).hexdigest()==h
original=base/'dim82-83-gamma72-signed-packets'
sys.argv=[str(old/'verify.py'),str(base)]
sys.path.insert(0,str(old))
K=runpy.run_path(sys.argv[0])
R,H=K['R']['R'],list(map(int,K['R']['H']))
D=dict(zip(map(tuple,R['A8']),R['V8']))
def positive(a,b,k):
 return ((a>=0)&(b>=0)&((a>0)|(b>0)))|((a>0)&(b<0)&(a*a>k*b*b))|((a<0)&(b>0)&(k*b*b>a*a))
def roots(i,d,n):
 g=[]
 for j in range(1,(d+1)//2):
  for a,b in product((-1,1),repeat=2):
   z=[0]*n;z[(i+j)%d]=a;z[(i-j)%d]=b;g.append(z)
 return g
assert len({tuple(z) for i in range(5)for z in roots(i,5,10)})==40
W=json.loads((P/'witness.json').read_text())
A=np.load(original/'tail82.npz',allow_pickle=False)['A'].astype('i8')
G=np.loadtxt(original/'tail83.txt',dtype='i8').reshape(-1,11,2)
def accept83(x):
 a=x@G[:,:,0].T; b=x@G[:,:,1].T
 q=(x*x).sum(1)
 return ~positive(a*a+2*b*b-18*q[:,None],2*a*b,2).any(1)
g83=[]
for i in range(11):
 X=np.asarray(roots(i,11,11),dtype='i8')
 g83.append(X[accept83(X)].tolist())
assert sum(map(len,g83))==208
c=W['83'];p=c['pairs'];omit=c['omitted']
assert len(p)==5 and all(len(r)==2 for r in p)
assert sorted([omit]+[i for o in p for i in o])==list(range(11))
labels=[tuple(t)for t in c['labels']]
assert len(labels)==16 and all(0<=j<5 and 0<=a<32 and not (a>>j&1)for j,a in labels)
assert {x for j,a in labels for x in (a,a^(1<<j))}==set(range(32))
assert all(i!=j or (a^b).bit_count()>=2 for (i,a),(j,b)in combinations(labels,2))
for j,a in labels:
 active=[p[k]for k in range(5)if k!=j]
 parity=[a>>k&1 for k in range(5)if k!=j]
 for eps in (0,1):
  packet=[]
  for z in product((0,1),repeat=4):
   if sum(z)%2!=eps:continue
   x=[0]*11
   for (u,v),b,t in zip(active,parity,z):x[u],x[v]=1-2*t,1-2*(t^b)
   packet.append(x)
  g83.append(packet)
assert len(g83)==43 and sum(map(len,g83))==464
g82=[roots(i,5,10)for i in range(5)]+W['82']
out={}
for n,groups in [('82',g82),('83',g83)]:
 d=int(n)
 assert all(groups) and len(groups)<len(H)
 x=np.asarray([v for g in groups for v in g],dtype='i8')
 I=np.array([i for i,g in enumerate(groups,start=1)for _ in g])
 assert x.shape[1]==d-72 and np.max(abs(x))<=6 and len(np.unique(x,axis=0))==len(x)
 assert all(H[i]!=0 for i in I)
 if d==82:
  wgt=np.array([3]*5+[1]*5,dtype='i8')
  q=(x*x*wgt).sum(1)
  # Physical tail is 2*(sqrt(3)*x_front,x_back)/sqrt(q).
  assert ((x@A.T)**2 <=24*q[:,None]).all()
  raw=(x*wgt)@x.T
 else:
  q=(x*x).sum(1)
  assert accept83(x).all()
  raw=x@x.T
 assert (q>0).all()
 # Same phase: tails have dot <=0. Different phases: tails have dot <=2.
 S=I[:,None]==I[None,:];eye=np.eye(len(x),dtype=bool)
 assert (raw[S&~eye]<=0).all()
 assert ((raw<=0)|(4*raw**2<=q[:,None]*q[None,:]))[~S].all()
 for i in set(I):
  V=D[0,H[i]]
  assert V.shape==(24,24) and np.array_equal(V@R['F']@V.T,8*np.eye(24,dtype='i8'))
 previous=json.loads((old/'claims.json').read_text())['lower_bounds'][n]
 out[n]=previous+144*len(x)
 print('PASS',n,out[n],'extra norm-four tails',len(x),'phase count',len(set(I)),flush=True)
assert out==json.loads((P/'claims.json').read_text())['lower_bounds']
