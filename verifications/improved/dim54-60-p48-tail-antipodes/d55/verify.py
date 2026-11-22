"""Exact E7 minuscule-tail append over the signed P48 construction."""
import argparse,gzip,json
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--siblings',type=Path,default=P.parent);a=ap.parse_args().siblings
D=a/'dim49-60-p48-signed-refill'; B=a/'dim49-60-p48-pointwise-refill'
def pin(path,sha):
 assert sha256(path.read_bytes()).hexdigest()==sha,path
 return path
u=np.load(pin(B/'heads55.npz','23d891f182411220d4e85afa39d457c34bd20a5df62e825f02756485830d52dd'),allow_pickle=False)
g=np.load(pin(B/'source/p48p_gram.npy','01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'),allow_pickle=False)
d=json.loads(gzip.decompress(pin(D/'p48_signed_delta55.json.gz','c74653824c84c961a5408207324e07fc5e132c4c4ac8828cb0e05ea0866fa87e').read_bytes()))
assert d['dimension']==55 and d['bound']==53393546 and np.array_equal(u['gram'],g)
def prod(x,y):
 assert x.dtype.kind in 'iu' and y.dtype.kind in 'iu'
 def m(a):return max(abs(int(a.min())),abs(int(a.max()))) if a.size else 0
 assert x.shape[-1]==y.shape[0] and x.shape[-1]*m(x)*m(y)<2**62
 return x.astype('i8')@y.astype('i8')
w=json.loads((P/'witness.json').read_text())
h=np.asarray(w['head'],dtype='i8');assert len(h)==48 and prod(prod(h,g),h)==6
V=u['base'];T=u['transforms'];E=u['extra_points'];O=u['extra_offsets']
M=np.unpackbits(u['retained_bits'],axis=1)[:,:len(V)].astype(bool);C={c['class_index']:c['add']for c in d['changes']}
omega=np.random.default_rng(41).integers(-10**7,10**7,48,dtype='i8'); H=abs(prod(h,omega))
for j,A in enumerate(T):
 t=prod(A,omega)
 for x in [V[M[j]],E[O[j]:O[j+1]],np.array(C.get(j,[]),dtype='i8').reshape(-1,48)]:
  assert not np.any(abs(prod(x,t))==H)
I=np.eye(8,dtype='i8')
E7=[]
for i,j in combinations(range(8),2):
 for p,q in product((-2,2),repeat=2):
  v=p*I[i]+q*I[j]
  if v[6]+v[7]==0:E7.append(v)
E7.extend(np.array(v,dtype='i8')for v in product((-1,1),repeat=8)if v[6]+v[7]==0 and sum(x<0 for x in v)%2==0)
E7=np.array(E7);assert len({tuple(v)for v in E7})==126 and np.all((E7**2).sum(1)==8)
Z=[np.array((*v,0,0),dtype='i8')for v in product((-1,1),repeat=6)if sum(x<0 for x in v)%2==0]
Z +=[a*I[i]+b*(I[6]-I[7])for i in range(6)for a,b in product((-2,2),(-1,1))]
Z=np.array(Z);assert len({tuple(v)for v in Z})==56 and np.all((Z**2).sum(1)==6) and np.all(Z[:,6]+Z[:,7]==0)
assert np.max(np.abs(prod(E7,Z.T)))<=4
J=prod(Z,Z.T);assert J[np.triu_indices(len(Z),1)].max()<=2
assert w['previous']==d['bound'] and w['add']==112 and w['lower_bound']==d['bound']+112
print(json.dumps(w))
