"""Exact replay of the D9 and even-octad Γ72 tail addition in 81D."""
from hashlib import sha256
from itertools import product
from pathlib import Path
import argparse,json,runpy,sys
import numpy as np
assert __debug__
P=Path(__file__).resolve().parent
par=argparse.ArgumentParser()
par.add_argument('--siblings',type=Path,default=P.parent)
par.add_argument('--incremental',action='store_true',help='Check the new geometry against the pinned predecessor; its replay is separate.')
a=par.parse_args()
base=a.siblings.resolve(); old=base/'dim81-84-antipodal-tail-codes'
for f,h in json.loads((P/'pins.json').read_text()).items():assert sha256((old/f).read_bytes()).hexdigest()==h
if not a.incremental:
 sys.argv=[str(old/'verify.py')]
 K=runpy.run_path(sys.argv[0])
 R=K['R']; H=list(map(int,R['H']))
 D=dict(zip(map(tuple,R['A8']),R['V8']))
 assert K['out']['81']==6251523810
A=np.asarray(json.loads((old/'P9a_tails.json').read_text()),dtype='i8')
assert A.shape==(153,9) and abs(A).max()<=2 and ((A*A).sum(1)==4).all()
groups=[]
for i in range(9):
 g=[]
 for j in range(1,5):
  for u,v in product((-1,1),repeat=2):
   x=[0]*9;x[(i+j)%9]=u;x[(i-j)%9]=v;g.append(x)
 groups.append(g)
assert sum(map(len,groups))==144
C={sum((((i&s).bit_count()^b)&1)<<i for i in range(8))for s in range(8)for b in range(2)}
assert len(C)==16 and {u.bit_count()for u in C}=={0,4,8}
U={u for u in range(256)if u.bit_count()%2==0}
while U:
 u=min(U); T={u^c for c in C}
 assert T<=U;U-=T
 groups.append([[1-2*(x>>i&1)for i in range(8)]+[0]for x in sorted(T)])
assert len(groups)==17
X=np.asarray([x for g in groups for x in g],dtype='i8')
J=np.array([i for i,g in enumerate(groups,start=1)for _ in g])
assert X.shape==(272,9) and len(np.unique(X,axis=0))==len(X)
q=(X*X).sum(1)
# Old mixed tails are A/sqrt(2), the new tails are 2*X/sqrt(q).
assert ((X@A.T)**2<=2*q[:,None]).all()
M=X@X.T;S=J[:,None]==J[None,:];eye=np.eye(len(X),dtype=bool)
assert (M[S&~eye]<=0).all()
assert ((M<=0)|(4*M**2<=q[:,None]*q[None,:]))[~S].all()
if not a.incremental:
 for i in set(J):
  assert H[i]!=0
  V=D[0,H[i]]
  assert V.shape==(24,24) and np.array_equal(V@R['F']@V.T,8*np.eye(24,dtype='i8'))
p=json.loads((old/'claims.json').read_text())['lower_bounds']['81']
assert p==6251523810 and p+144*len(X)==6251562978
assert json.loads((P/'claims.json').read_text())['lower_bounds']['81']==p+144*len(X)
print('PASS 81D 6,251,562,978:',len(X),'norm-four tails in',len(groups),'phases')
