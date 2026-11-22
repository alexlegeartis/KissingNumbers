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
delta_file=pin(D/'p48_signed_delta60.json.gz','6e516dedf722e13c6c660a873285332d73ee0f813ad8cfd39e2d97c5715612ac')
G=np.load(gram_file,allow_pickle=False).astype('i8')
U=np.load(head_file,allow_pickle=False)
config=json.loads(gzip.decompress(delta_file.read_bytes()))
W=json.loads((HERE/'witness.json').read_text())
assert config['dimension']==60 and config['bound']==58631571
assert G.shape==(48,48) and np.array_equal(U['gram'],G)
Hs=np.array(W,dtype='i8')
V=U['base'].astype('i8');Q=U['extra_points'].astype('i8');T=U['transforms'].astype('i8')
O=U['extra_offsets'];bits=np.unpackbits(U['retained_bits'],axis=1)[:,:len(V)].astype(bool)
C={x['class_index']:x['add']for x in config['changes']}
omega=np.random.default_rng(41).integers(-10**7,10**7,48,dtype='i8')
def m(a):return max(map(abs,map(int,a.flat)),default=0)
assert all(v.dtype.kind in 'iu' for v in (G,V,Q,T))
def ip(a,b):
 assert a.shape[-1]==b.shape[0] and a.shape[-1]*m(a)*m(b)<2**62
 return a@b
assert np.array_equal(ip(ip(Hs,G),Hs.T),6*np.eye(len(W),dtype='i8'))
target=set(np.abs(ip(Hs,omega)).tolist())
for j,Uj in enumerate(T):
 Cj=ip(Uj,omega)
 for X in (V[bits[j]],Q[O[j]:O[j+1]],np.array(C.get(j,[]),dtype='i8').reshape(-1,48)):
  assert target.isdisjoint(np.abs(ip(X,Cj)).tolist()), ('ownership',j)

from itertools import product
from collections import Counter
from functools import reduce
from operator import xor
assert Hs.shape==(len(W),48) and 1<=len(W)<=16
supports=[]
for j in (0,6):
 for k in (0,2,4):
  supports.append(tuple(j+i for i in range(6) if i not in (k,k+1)))
for c in range(5):
 M=[(c,5),((c+1)%5,(c-1)%5),((c+2)%5,(c-2)%5)]
 for e,f in product(M,repeat=2): supports.append(tuple(e)+tuple(j+6 for j in f))
D=[]
for i in range(12):
 for a in (-2,2):
  v=[0]*12;v[i]=a;D.append(v)
assert len(supports)==51 and len(set(tuple(sorted(e)) for e in supports))==51
for e in supports:
 for signs in product((-1,1),repeat=4):
  v=[0]*12
  for i,a in zip(e,signs,strict=True):v[i]=a
  D.append(v)
D=np.array(D,dtype='i8'); M=D@D.T
assert len(D)==840 and np.all(np.diag(M)==4)
assert np.all(M[np.triu_indices(840,1)]<=2)
assert all(tuple(-v) in set(map(tuple,D)) for v in D)
assert np.max(np.abs(D).sum(1))<=4
balls={}
for w in range(1<<12):
 if w.bit_count()%2==0:
  u=reduce(xor,(j for j in range(12) if (w>>j)&1),0)
  balls.setdefault(u,[]).append(w)
assert set(balls)==set(range(16))
for words in balls.values():
 assert len(words)==128 and all((a^b).bit_count()>=4 for a,b in combinations(words,2))
# All new points are (±u/2, sqrt(3) z/2), with u a unit unowned head
# and z a unit sign tail. Old poles have |z·pole|<=1/sqrt(3) by L1.
# Old lifts score <=sqrt(3)/8+1/4 < 1/2 by ownership.
# Same new color scores <=1/4+(3/4)(1/3); different colors:
# orthogonal heads and distance >=2 imply score <=(3/4)(2/3).
assert F(1,4)+F(3,4)*F(1,3)==F(3,4)*F(2,3)==F(1,2)
print(json.dumps(dict(previous=config['bound'],head_colors=len(W),remove_pole_net=1,
                      additions=256*len(W),lower_bound=config['bound']-1+256*len(W))))
