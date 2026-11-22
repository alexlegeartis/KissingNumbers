"""Literal tetrahedral pole extension in dimensions 55, 57 and 59."""
import argparse,gzip,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent
ap=argparse.ArgumentParser()
ap.add_argument('--siblings',type=Path,default=P.parent)
ap.add_argument('--dims',type=int,nargs='+',default=[55,57,59])
a=ap.parse_args()
D=a.siblings/'dim49-60-p48-signed-refill'; B=a.siblings/'dim49-60-p48-pointwise-refill'
def pin(path,sha):
 assert sha256(path.read_bytes()).hexdigest()==sha,path
 return path
def m(a):return max(abs(int(a.min())),abs(int(a.max()))) if a.size else 0
def prod(x,y):
 assert x.dtype.kind in 'iu' and y.dtype.kind in 'iu'
 assert x.shape[-1]==y.shape[0] and x.shape[-1]*m(x)*m(y)<2**62
 return x.astype('i8')@y.astype('i8')
TABLE={
 55:('23d891f182411220d4e85afa39d457c34bd20a5df62e825f02756485830d52dd','c74653824c84c961a5408207324e07fc5e132c4c4ac8828cb0e05ea0866fa87e',None,53393546,126,168),
 57:('b081a40b1845dbf3e10448328449c6449916cec7771736defe45598591a517fe','ea9a5e9a92059fffb9b8e4e0ba93a3d616b8de2c11bd7a6b257f869fb5dea64f','c5a711ad2f02482b3af198896dfe74dc2a3ab07bf0bf66e4d80b9758584ef653',54760693,264,352),
 59:('e92d40b8e36bbfdd81b744fa73692ddbf96cb099f0ce91c6475d8986b9024ab1','cd0c5b368828c13117fa845b97b12b8a930333a86893cb8bb5c0ef4d2446b708','62e21350d0f62718a16a012b8932caf9c5cb16db4f2bd00e99c8f5934e18ba62',56952414,168,224)}
g=np.load(pin(B/'source/p48p_gram.npy','01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'),allow_pickle=False)
cfg=json.loads((P/'witness.json').read_text())
def sign(a,b,R):
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if a*b>0:return 1 if a>0 else -1
 q=a*a-R*b*b
 return ((q>0)-(q<0))*(1 if a>0 else -1)
for n in a.dims:
 w=cfg[str(n)];k=n-48
 q0,q1,q2,before,oldgain,more=TABLE[n];nc=14 if n==55 else 8
 u=np.load(pin(B/f'heads{n}.npz',q0),allow_pickle=False)
 d=json.loads(gzip.decompress(pin(D/f'p48_signed_delta{n}.json.gz',q1).read_bytes()))
 assert d['dimension']==n and d['bound']==before and np.array_equal(u['gram'],g)
 h=np.asarray(w['heads'],dtype='i8');assert h.shape==(nc,4,48)
 for Tj in h: assert not Tj.sum(axis=0).any() and np.array_equal(prod(prod(Tj,g),Tj.T),8*np.eye(4,dtype='i8')-2)
 h=h.reshape(-1,48);assert len({tuple(z)for z in h})==4*nc
 V=u['base'];T=u['transforms'];E=u['extra_points'];O=u['extra_offsets']
 M=np.unpackbits(u['retained_bits'],axis=1)[:,:len(V)].astype(bool);C={c['class_index']:c['add']for c in d['changes']}
 omega=np.random.default_rng(41).integers(-10**7,10**7,48,dtype='i8'); H=abs(prod(h,omega))
 for j,A in enumerate(T):
  t=prod(A,omega)
  for x in [V[M[j]],E[O[j]:O[j+1]],np.array(C.get(j,[]),dtype='i8').reshape(-1,48)]:
   assert set(abs(prod(x,t)).tolist()).isdisjoint(H.tolist())
 if n==55:
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
  norm=6
 else:
  U=np.load(pin(B/f'source/capdirs_{k}.npz',q2),allow_pickle=False)
  A=U['A'];Br=U['B'];R=int(U['R'].item()); N=int(U['norm'].item())
  assert not np.any(U['c']!=1)
  S=[tuple(t)for t in w['supports']]
  assert len(S)==len(set(S)) and all(len(t)==3 and len(set(t))==3 and all(0<=v<k for v in t)for t in S)
  Z=[]
  for x in S:
   for sig in product((-1,1),repeat=3):
    v=np.zeros(k,dtype='i8');v[list(x)]=sig;Z.append(v)
  Z=np.array(Z); assert np.all((Z*Z).sum(1)==3)
  norm=3
  for p,q in zip(prod(A,Z.T).flat,prod(Br,Z.T).flat):
   p=int(p);q=int(q);assert sign(p*p+R*q*q-N,2*p*q,R)<=0
 J=prod(Z,Z.T);assert 3*J[np.triu_indices(len(Z),1)].max()<=norm
 C=np.asarray(w['colors'],dtype='i8');assert C.shape==(len(Z),) and set(C)==set(range(nc))
 for i,j in combinations(range(len(Z)),2):
  if C[i]==C[j]:assert 5*int(J[i,j])<=norm
 assert F(3,8)*(-F(1,3))+F(5,8)==F(1,2)
 assert F(3,8)+F(5,8)*F(1,5)==F(1,2)
 assert F(3,8)*F(1,2)+F(5,8)*F(1,3)<F(1,2)
 assert F(5,8)*F(1,3)<F(1,4)
 assert F(31,8)**2>15 and F(7,5)**2<2 and (F(3,8)+F(31,96))/F(7,5)<F(1,2)
 assert 4*len(Z)-len(h)==more
 print(json.dumps(dict(dimension=n,previous=before+oldgain,add=more-oldgain,
                       lower_bound=before+more)),flush=True)
