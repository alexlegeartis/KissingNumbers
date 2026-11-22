"""Read-only exact verifier: a matched norm-12/12D tail augmentation in dimension 60."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from itertools import product,combinations
from functools import reduce
from operator import xor
import argparse,gzip,json,numpy as np
if not __debug__:raise SystemExit('Assertions required.')
P=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--siblings',type=Path,default=P.parent)
p.add_argument('--body',type=Path,default=None)
p.add_argument('--refill',type=Path,default=None)
p.add_argument('--witness',type=Path,default=P/'witness.json')
a=p.parse_args()
A=a.siblings/'dim49-60-p48-pointwise-refill';B=a.siblings/'dim49-60-p48-signed-refill'
body=a.body or a.siblings/'dim50-60-p48-global-owner-refill'
refill=a.refill or body/'refill60.json.gz'
def pin(path,d):
 z=path.read_bytes();assert sha256(z).hexdigest()==d,(str(path),sha256(z).hexdigest())
 return z
def readz(path,d):return json.loads(gzip.decompress(pin(path,d)))
def num(a):
 return max(abs(int(a.min())),abs(int(a.max()))) if a.size else 0
def ip(a,b):
 assert a.dtype.kind in 'iu' and b.dtype.kind in 'iu'
 assert a.shape[-1]==b.shape[0] and a.shape[-1]*num(a)*num(b)<2**62
 return a.astype('i8')@b.astype('i8')
G=np.load(A/'source/p48p_gram.npy',allow_pickle=False)
pin(A/'source/p48p_gram.npy','01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1')
pin(A/'heads60.npz','cc3a2e07fbe3f5f91c33d61ba69102057c5fd8f90dc509d3c5d83d29adac9b25')
U=np.load(A/'heads60.npz',allow_pickle=False)
assert np.array_equal(U['gram'],G)
d=readz(B/'p48_signed_delta60.json.gz','6e516dedf722e13c6c660a873285332d73ee0f813ad8cfd39e2d97c5715612ac')
g=readz(body/'global60.json.gz','841ec616c18efe427f7f5bd1665a79b5e8e0c57538606fa14141e183e3276393')
f=readz(refill,'6f34b349b3458a946dc3b90129245542036ec66374fba93c42be777711da95bc')
bank=np.array(readz(body/'bank.json.gz','65d79e04fe0164c344934ed8231b51faf7d316da101093235781cd129d8b3e9f'),dtype='i8')
root=json.loads(pin(body/'tail60/witness.json','d203b2ec80be235d57a7ff37fcfb63c3c2361dac5df60cb3e528499bb76d70e1'))
assert len(g['changes'])==420 and d['bound']==58631571 and bank.shape==(4000,48)
assert g['adds']+f['adds']==62937 and g['owner_gain']+f['gain']==27016
Y={w['class_index']:bank[w['insert']]for w in g['changes']}
for w in f['changes']:
 j=w['class_index'];Y[j]=np.r_[Y[j],np.array(w['extra'],dtype='i8').reshape(-1,48)]
assert sum(map(len,Y.values()))==62937
# Old 840-point antipodal stock, scaled by 2, rederived rather than trusted.
ds=[]
for i in range(12):
 for a0 in (-2,2):
  v=[0]*12;v[i]=a0;ds.append(v)
supports=[]
for j in (0,6):
 for k in (0,2,4):supports.append(tuple(j+i for i in range(6)if i not in (k,k+1)))
for c in range(5):
 Q=[(c,5),((c+1)%5,(c-1)%5),((c+2)%5,(c-2)%5)]
 for e,r in product(Q,repeat=2):supports.append(tuple(e)+tuple(j+6 for j in r))
for e in supports:
 for q in product((-1,1),repeat=4):
  v=[0]*12
  for j,value in zip(e,q,strict=True):v[j]=value
  ds.append(v)
D=np.array([v for v in ds if next(t for t in v if t)>0],dtype='i8')
assert D.shape==(420,12)
H=np.array(root['frame'],dtype='i8');S=np.array(root['sparse'],dtype='i8')
w=json.loads(a.witness.read_text())
W=np.array(w['heads'],dtype='i8')
Z=np.array([q['z']for q in w['tail']],dtype='i8')
I=np.array([q['h']for q in w['tail']],dtype='i8')
perm=np.array(w['permutation'],dtype='i8')
assert H.shape==(16,48) and S.shape==(4,4,48) and W.ndim==2 and W.shape[1]==48
assert Z.shape==(len(I),12) and set(I)==set(range(len(W))) and len(Z)==len(set(map(tuple,Z)))
assert np.all(np.sum(ip(W,G)*W,axis=1)==12)
assert np.all(np.sum(Z*Z,axis=1)==10) and np.max(np.abs(Z).sum(axis=1))<=8
assert not np.any(ip(ip(W,G),H.T))
assert np.max(abs(ip(Z,D.T)))<=4
mat=np.abs(ip(ip(W[I],G),W[I].T))+2*ip(Z,Z.T)
assert np.max(mat[np.triu_indices(len(Z),1)])<=16
assert F(5,8)-F(3,8)<=F(1,2)
# All raw w·x with x norm6 lie in [-6,6] by the mother min6;
# in particular against equator 6/(8 sqrt3) < 1/2.
assert F(6,8)**2<3*F(1,2)**2
# Check against the exact literal sparse mode from the pinned parent tail.
assert F(141,100)**2<2 and F(457,1000)**2>F(5,24)
for r,c in enumerate([6,6,6,3,0]):
 assert F(c,16)/F(141,100)+F(r,4)*F(457,1000)<=F(1,2)
for c in range(4):
 U0=np.max(np.abs(ip(ip(W,G),S[c].T)),axis=1)[I]
 inds=[[j+2*k+((q>>k)&1)for k in range(3)]for j in (0,6)for q in (c,c^7)]
 r=np.max(np.array([abs(Z[:,j]).sum(axis=1)for j in inds]),axis=0)
 assert np.max(r)<=4 and np.all(U0<=np.array([6,6,6,3,0])[r])
# Against any unupgraded cube fibre compatibility follows from exact W·H=0
# and the norm-L1 bound. Check all actual upgraded tetra subcosets.
balls={}
for x in range(4096):
 if x.bit_count()%2:continue
 k=reduce(xor,(i for i in range(12)if x>>i&1),0);balls.setdefault(k,[]).append(x)
L=set(root['subcode']);assert len(L)==16
assert F(547,100)**2<30 and F(223,100)**2<5
for r,c in ((8,1),(6,3),(4,6)):
 assert F(c,4)/F(547,100)+F(r,8)/F(223,100)<=F(1,2)
k=0
for u in root['upgrade']:
 left=set(balls[u['color']])
 for x in u['tetrahedra']:
  p0=min(left);coset={p0^l for l in L};assert coset<=left;left-=coset
  if x is None:continue
  k+=1
  Dn=np.array([[1-2*((q>>i)&1)for i in range(12)]for q in sorted(coset)],dtype='i8')
  r=np.max(np.abs(ip(Z,Dn.T)),axis=1)
  assert np.all((r%2)==0) and np.max(r)<=8
  U0=np.max(np.abs(ip(ip(W,G),np.array(x,dtype='i8').T)),axis=1)[I]
  assert np.all(U0<=np.where(r>=8,1,np.where(r>=6,3,6)))
 assert not left
# Match old lift direction labels at the new exact body-shell contacts.
# Old lift (x/(2 sqrt2), ±D/4), new (±w/(4 sqrt2), Z/4),
# so the cap is exactly |w·x|+|Z·D|<=8.
assert sorted(perm.tolist())==list(range(420))
scores=np.array([np.max(np.abs(ip(Z[I==i],D.T)),axis=0)for i in range(len(W))])
V=U['base'];E=U['extra_points'];T=U['transforms'];O=U['extra_offsets']
M=np.unpackbits(U['retained_bits'],axis=1)[:,:len(V)].astype(bool)
d0={q['class_index']:q['add']for q in d['changes']}
WG=ip(G,W.T)
cost=0; class_costs=[]
for j,Hj in enumerate(T):
 # The pinned pointwise archive stores one representative per antipodal line.
 # Expand both signs here; the signed delta and all later changes are already oriented.
 old=np.vstack((V[M[j]],E[O[j]:O[j+1]]))
 x=np.vstack((old,-old,np.array(d0.get(j,[]),dtype='i8').reshape(-1,48),Y[j]))
 y=ip(Hj,WG)
 assert 48*num(x)*num(y)<2**52
 # Exactly rounded integers even for every partial sum in a binary64 GEMM.
 m=np.abs(x.astype(float)@y.astype(float)); b=np.any(m+scores[:,perm[j]][None,:]>8,axis=1);cost+=int(b.sum());class_costs.append(int(b.sum()))
pre=d['bound']+g['owner_gain']+f['gain']-1+4096+480+28*k
assert k==20 and pre==58663722
claim=dict(status='PASS',dimension=60,previous=pre,carriers=len(W),tails=len(Z),charged_blockers=cost,gain=2*len(Z)-cost,lower_bound=pre+2*len(Z)-cost)
if a.witness.resolve()==(P/'witness.json').resolve():
 assert claim==json.loads((P/'claims.json').read_text())
print(json.dumps(claim),flush=True)
