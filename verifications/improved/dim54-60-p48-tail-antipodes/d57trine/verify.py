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
def prod(x,y):
 assert x.dtype.kind in 'iu' and y.dtype.kind in 'iu'
 def m(a):return max(abs(int(a.min())),abs(int(a.max()))) if a.size else 0
 assert x.shape[-1]==y.shape[0] and x.shape[-1]*m(x)*m(y)<2**62
 return x.astype('i8')@y.astype('i8')
TABLE={
 57:('b081a40b1845dbf3e10448328449c6449916cec7771736defe45598591a517fe','ea9a5e9a92059fffb9b8e4e0ba93a3d616b8de2c11bd7a6b257f869fb5dea64f','c5a711ad2f02482b3af198896dfe74dc2a3ab07bf0bf66e4d80b9758584ef653',54760693,264),
 59:('e92d40b8e36bbfdd81b744fa73692ddbf96cb099f0ce91c6475d8986b9024ab1','cd0c5b368828c13117fa845b97b12b8a930333a86893cb8bb5c0ef4d2446b708','62e21350d0f62718a16a012b8932caf9c5cb16db4f2bd00e99c8f5934e18ba62',56952414,128)}
gram_sha='01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'
g=np.load(pin(B/'source/p48p_gram.npy',gram_sha),allow_pickle=False)
cfg=json.loads((P/'witnesses.json').read_text())
def sign(a,b,R):
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if a*b>0:return 1 if a>0 else -1
 q=a*a-R*b*b
 return ((q>0)-(q<0))*(1 if a>0 else -1)
for w in cfg:
 n=w['dimension'];k=n-48
 q0,q1,q2,before,more=TABLE[n]
 u=np.load(pin(B/f'heads{n}.npz',q0),allow_pickle=False)
 d=json.loads(gzip.decompress(pin(D/f'p48_signed_delta{n}.json.gz',q1).read_bytes()))
 assert d['dimension']==n and d['bound']==before and np.array_equal(u['gram'],g)
 h=np.asarray(w['head_triples'],dtype='i8');assert h.shape==(8,3,48)
 for Tj in h: assert not Tj.sum(axis=0).any() and np.array_equal(prod(prod(Tj,g),Tj.T),9*np.eye(3,dtype='i8')-3*np.ones((3,3),dtype='i8'))
 h=h.reshape(-1,48);assert len({tuple(z)for z in h})==24
 V=u['base'];T=u['transforms'];E=u['extra_points'];O=u['extra_offsets']
 M=np.unpackbits(u['retained_bits'],axis=1)[:,:len(V)].astype(bool);C={c['class_index']:c['add']for c in d['changes']}
 omega=np.random.default_rng(41).integers(-10**7,10**7,48,dtype='i8'); H=abs(prod(h,omega))
 for j,A in enumerate(T):
  t=prod(A,omega)
  for x in [V[M[j]],E[O[j]:O[j+1]],np.array(C.get(j,[]),dtype='i8').reshape(-1,48)]:
   assert set(abs(prod(x,t)).tolist()).isdisjoint(H.tolist())
 U=np.load(pin(B/f'source/capdirs_{k}.npz',q2),allow_pickle=False)
 A=U['A'];Br=U['B'];R=int(U['R'].item()); N=int(U['norm'].item())
 assert not np.any(U['c']!=1)
 S=[tuple(t)for t in w['supports']]
 assert len(S)==len(set(S)) and all(len(t)==3 and len(set(t))==3 and all(0<=v<k for v in t)for t in S)
 assert all(len(set(x)&set(y))<=1 for x,y in combinations(S,2))
 Z=[]
 for x in S:
  for sig in product((-1,1),repeat=3):
   v=np.zeros(k,dtype='i8');v[list(x)]=sig;Z.append(v)
 Z=np.array(Z);assert np.all((Z*Z).sum(1)==3)
 for p,q in zip(prod(A,Z.T).flat,prod(Br,Z.T).flat):
  p=int(p);q=int(q)
  assert sign(p*p+R*q*q-N,2*p*q,R)<=0
 Q=prod(Z,Z.T);assert Q[np.triu_indices(len(Z),1)].max()<=1
 C=np.asarray(w['colors'],dtype='i8');assert C.shape==(len(Z),) and set(C)==set(range(8))
 for i,j in combinations(range(len(Z)),2):
  if C[i]==C[j]:assert 4*int(Q[i,j])<=3
 assert 3*len(Z)-len(h)==more and w['lower_bound']==before+more
 print(json.dumps(dict(dimension=n,previous=before,add=more,lower_bound=before+more,tail_supports=S)))
