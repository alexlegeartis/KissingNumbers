"""Exact pure-layer replacement on the frozen 59D unequal-coordinate cube code."""
import argparse,json,itertools,io
from pathlib import Path
from hashlib import sha256
import numpy as np
P=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--parent',type=Path,default=P.parent/'dim59-p48-frame-cube',help='weighted-cube certificate directory')
p.add_argument('--owners',type=Path,default=P.parent/'dim49-60-p48-pointwise-refill',help='pointwise owner certificate directory')
p.add_argument('--body',type=Path,default=P.parent/'dim50-60-p48-global-owner-refill',help='global body certificate directory')
p.add_argument('--witness',type=Path,default=P/'witness.json')
o=p.parse_args()
if not __debug__:raise SystemExit('Use assertions.')
def pin(a,h):
 z=a.read_bytes();assert sha256(z).hexdigest()==h;return z
caps=np.load(io.BytesIO(pin(o.owners/'source/capdirs_11.npz','62e21350d0f62718a16a012b8932caf9c5cb16db4f2bd00e99c8f5934e18ba62')),allow_pickle=False)
tet=json.loads(pin(o.body/'tail55/witness.json','cc4151354d9b6eaf6551275a2d9424eb5b8eba11e0a4c3bf5db0bec31e9aa053'))['59']
old=json.loads(pin(o.parent/'witness.json','68f8687a602a51559bdbfecac75cacf90097a47e49b70bfdf711091c22527cf7'))
Y=np.array(json.loads(o.witness.read_text())['poles'],dtype='i8')
n,d=Y.shape;assert d==11 and n>0 and abs(Y).max()<1000
N=(Y*Y).sum(1);assert (N>0).all()
J=Y@Y.T;k,l=np.triu_indices(n,1);assert (4*J[k,l]**2<=N[k]*N[l]).all()
cube=np.array([u for u in itertools.product((-1,1),repeat=11)if u.count(-1)%2==0],dtype='i8')*np.r_[np.full(8,50),np.full(3,49)]
assert len(cube)==1024
L=np.abs(Y)@np.r_[np.full(8,50),np.full(3,49)]
assert (3*L*L<=27203*N).all()
supp=np.array(tet['supports'],dtype='i8');assert supp.shape==(8,3)
K=np.abs(Y[:,supp]).sum(2);assert (5*K*K<=6*N[:,None]).all()
Z=np.array(old['weighted_tails'],dtype='i8')[:,1:];scale=int(old['scale'])
K=Y@Z.T;assert (4*K*K<=N[:,None]*scale**2).all()
A=caps['A'];B=caps['B']
def nn(a,b):
 # exact nonnegative test a+b sqrt(2)>=0
 if a>=0 and b>=0:return True
 if a<0 and b<0:return False
 if a>=0:return a*a>=2*b*b
 return a*a<=2*b*b
D=cube@A.T;E=cube@B.T;poles=set(range(604))
for j,(a,b)in enumerate(zip(D.T,E.T)):
 for x,y in zip(a.tolist(),b.tolist()):
  if nn(x,y) and not nn(12*27203-x*x-2*y*y,-2*x*y):
   poles.remove(j);break
assert len(poles)==28
T=Y@A.T;S=Y@B.T
remaining=[]
for j in sorted(poles):
 if all(nn(9*int(q)-a*a-2*b*b,-2*a*b) for q,a,b in zip(N,T[:,j].tolist(),S[:,j].tolist())):remaining.append(j)
delta=2*n+len(remaining)-28
assert delta>0
print(json.dumps(dict(status='PASS',prior=56969091,old_pure_retained=len(remaining), new_pure=2*n,improvement=delta,lower_bound=56969091+delta)))

assert 56969091+delta == json.loads((P/'claims.json').read_text())['lower_bounds']['59']
