#!/usr/bin/env python3
"""Exact P48 pole splitting, no optimization."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from math import lcm
import argparse,json,gzip,numpy as np

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--siblings',type=Path,default=HERE.parent)
args=parser.parse_args()
P=args.siblings/'dim49-60-p48-pointwise-refill'
D=args.siblings/'dim49-60-p48-signed-refill'

def pin(path, digest):assert sha256(path.read_bytes()).hexdigest()==digest,path;return path
head_file=pin(P/'heads60.npz','cc3a2e07fbe3f5f91c33d61ba69102057c5fd8f90dc509d3c5d83d29adac9b25')
gram_file=pin(P/'source/p48p_gram.npy','01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1')
tail_file=pin(P/'source/tail840_exact.npz','39d71d68bd77cfb96fd973a70b48059ef25377d3a73b6dbcb82de59f93779009')
pole_file=pin(P/'source/coordinates841.txt','be4b651ae8167e2c33bda881f77c9fa2cb5e6a07a8fb92f8360acd4691ed48de')
delta_file=pin(D/'p48_signed_delta60.json.gz','6e516dedf722e13c6c660a873285332d73ee0f813ad8cfd39e2d97c5715612ac')
G=np.load(gram_file,allow_pickle=False).astype('i8')
U=np.load(head_file,allow_pickle=False)
A=np.load(tail_file,allow_pickle=False)
config=json.loads(gzip.decompress(delta_file.read_bytes()))
W=json.loads((HERE/'candidate.json').read_text())
assert config['dimension']==60 and config['bound']==58631571
assert G.shape==(48,48) and np.array_equal(U['gram'],G)
poles=[]
for l in pole_file.read_text().splitlines():
 if not l.strip():continue
 a=list(map(F,l.split()));d=lcm(*(z.denominator for z in a));poles.append([int(z*d)for z in a])
assert len(poles)==841 and all(len(p)==12 for p in poles)
def dot(x,y):return sum(int(a)*int(b)for a,b in zip(x,y,strict=True))
def sign(a,b):
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if a*b>0:return 1 if a>0 else -1
 q=a*a-2*b*b
 return ((q>0)-(q<0))*(1 if a>0 else -1)
def check(x,y,t):
 s=dot(x,y)
 return s<=0 or F(s*s)<=t*t*dot(x,x)*dot(y,y)
ids={x['drop_pole']for x in W}
assert len(ids)==len(W)>0 and all(0<=i<841 for i in ids)
assert all(len(w['head'])==48 and len(w['tail'])==12 for w in W)
# All needed split heads are outside even the enlarged old-class union;
# local coefficient fingerprints can only give false negatives, not false positives.
Hs=np.array([w['head']for w in W],dtype='i8')
V=U['base'].astype('i8');Q=U['extra_points'].astype('i8');T=U['transforms'].astype('i8')
O=U['extra_offsets'];bits=np.unpackbits(U['retained_bits'],axis=1)[:,:len(V)].astype(bool)
C={x['class_index']:x['add']for x in config['changes']}
omega=np.random.default_rng(41).integers(-10**7,10**7,48,dtype='i8')
def m(a):return max(map(abs,map(int,a.flat)),default=0)
assert all(v.dtype.kind in 'iu' for v in (G,V,Q,T))
def ip(a,b):
 assert a.shape[-1]==b.shape[0] and a.shape[-1]*m(a)*m(b)<2**62
 return a@b
assert np.array_equal(np.diag(ip(ip(Hs,G),Hs.T)),np.full(len(W),6))
target=set(np.abs(ip(Hs,omega)).tolist())
for j,Uj in enumerate(T):
 Cj=ip(Uj,omega)
 for X in (V[bits[j]],Q[O[j]:O[j+1]],np.array(C.get(j,[]),dtype='i8').reshape(-1,48)):
  assert target.isdisjoint(np.abs(ip(X,Cj)).tolist()), ('ownership',j)
S=F('0.6547')
assert 3*(S+F(1,2))**2<4
A0,A1=A['A'],A['B'];assert A0.shape==A1.shape==(840,12)
for w in W:
 x=w['tail'];xx=dot(x,x);assert xx>0
 for a,b in zip(A0,A1,strict=True):
  s,t=dot(a,x),dot(b,x)
  assert sign(F(s*s+2*t*t)-16*S*S*xx,2*s*t)<=0,('lift direction',w['drop_pole'])
 for j,p in enumerate(poles):
  if j in ids:continue
  s=dot(x,p)
  assert s<=0 or 3*s*s<=xx*dot(p,p),('pole',w['drop_pole'],j)
for i,j in combinations(range(len(W)),2):
 a,b=Hs[i],Hs[j];raw=abs(int(ip(ip(a,G),b)))
 assert check(W[i]['tail'],W[j]['tail'],F(2,3)-F(raw,18)),('split pair',i,j)
assert F(3,4)-F(1,4)==F(1,2)
print(json.dumps(dict(previous=config['bound'],split=len(W),remove=sorted(ids),add=2*len(W),lower_bound=config['bound']+len(W))))
