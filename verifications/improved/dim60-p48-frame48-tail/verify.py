"""Read-only exact verifier: 48-axis norm-12 frame and 12D tail layers in dimension 60."""
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
w=json.loads(a.witness.read_text())
W=np.array(w['frame'],dtype='i8')
controls=np.array(w['controls'],dtype='i8')
active=np.array(w['carriers'],dtype='i8')
assert W.shape==(48,48) and np.array_equal(ip(ip(W,G),W.T),12*np.eye(48,dtype='i8'))
assert np.all((W-W[0])%2==0)
assert controls.shape==(16,) and active.shape==(32,) and sorted(np.r_[controls,active].tolist())==list(range(48))
H=np.vstack(((W[controls[0::2]]+W[controls[1::2]])//2,(W[controls[0::2]]-W[controls[1::2]])//2))
A0=W[active];WG=ip(A0,G)
assert H.shape==(16,48) and np.array_equal(ip(ip(H,G),H.T),6*np.eye(16,dtype='i8'))
assert not np.any(ip(WG,H.T))

groups=[np.array(q,dtype='i8')for q in w['tail_groups']]
assert len(groups)==32 and all(q.ndim==2 and q.shape[1]==12 and len(q) for q in groups)
Z=np.concatenate(groups)
I=np.repeat(np.arange(32),list(map(len,groups)))
assert len(set(map(tuple,Z)))==len(Z)
assert np.all(np.sum(Z*Z,axis=1)==10) and np.max(np.abs(Z).sum(1))<=8
assert np.max(abs(ip(Z,D.T)))<=4
g0=ip(Z,Z.T)*2+12*(I[:,None]==I[None,:])
assert np.max(g0[np.triu_indices(len(Z),1)])<=16
assert F(12,32)+F(10,16)==1 and F(10,16)-F(12,32)<=F(1,2)
# Equator bound 6/(8 sqrt(3)) (lattice minimum six gives |w·x|<=6).
assert F(6,8)**2 < 3*F(1,2)**2

# Dense replacement: all even-word sign cube points ±h/(2sqrt(6))+eps/4.
balls={}
for v in range(4096):
 if v.bit_count()%2:continue
 x=reduce(xor,(k for k in range(12)if (v>>k)&1),0)
 balls.setdefault(x,[]).append(v)
assert set(balls)==set(range(16))
assert all(len(b)==128 and all((i^j).bit_count()>=4 for i,j in combinations(b,2))for b in balls.values())
# Same-colour cube points use t²=.75 and d>=4; other colours are
# orthogonal with even d>=2. The 12-dimensional pole cap is |D·eps|/8<=1/2.
assert F(1,4)+F(3,4)*F(1,3)==F(1,2)==F(3,4)*F(2,3)
assert np.max(abs(D).sum(1))==4
# New cube against legacy lifted owner: raw <=3 because head rays unowned,
# hence maximum 3/(8sqrt(3))+1/4<1/2.
assert F(3,8)**2 < 3*F(1,4)**2
# New cube against q10: zero head product and |Z·eps|<=8.

# Keep any independently compatible old sparse heads, from the pinned
# tetrahedral sparse stock. Every kept h contributes 32 signed tail rays
# in the four support/sign groups and costs two equatorial rays.
keep=np.array(w.get('sparse_keep',[[False]*4]*4))
assert keep.shape==(4,4) and keep.dtype==bool
spar=np.array(root['sparse'],dtype='i8')
assert spar.shape==(4,4,48)
HS=[];ZS=[]
for c in range(4):
 if not keep[c].any():continue
 h=spar[c,keep[c]]
 for b in (0,6):
  for p0 in (c,c^7):
   pos=[b+2*k+((p0>>k)&1) for k in range(3)]
   x=abs(ip(WG,h.T))[I].max(1)
   y=abs(Z[:,pos]).sum(1)
   assert F(141421356,10**8)**2<2 and F(456435465,10**9)**2>F(5,24)
   assert all(F(int(u),16)/F(141421356,10**8)+F(int(v),4)*F(456435465,10**9)<=F(1,2)for u,v in zip(x,y))
   for s in product((-1,1),repeat=3):
    a0=np.zeros(12,dtype='i8');a0[pos]=s
    b0=-1 if sum(t<0 for t in s)%2 else 1
    HS.extend((b0*h).tolist());ZS.extend([a0.tolist()]*len(h))
S0=np.array(HS,dtype='i8').reshape(-1,48);ZS=np.array(ZS,dtype='i8').reshape(-1,12)
assert len(S0)==32*int(keep.sum()) and len(S0)==len(ZS)
if len(S0):
 assert np.all((ip(S0,G)*S0).sum(1)==6)
 assert np.max(abs(ip(ZS,D.T)))<=2
 m=3*ip(ip(S0,G),S0.T)+10*ip(ZS,ZS.T)
 assert np.max(m[np.triu_indices(len(S0),1)])<=24
# Sparse old body cap 3/(8sqrt(2))+sqrt(5/24)/2<1/2
assert F(3,8)/F(141421356,10**8)+F(456435465,10**9)/2<=F(1,2)
# Sparse/cube: tail unit inner <=1/2, and distinct unowned norm6 head
# lines have unit inner <=1/2, so Cauchy–Schwarz gives the mixed bound.
check=np.r_[H,spar[keep]]
assert np.max(abs(ip(ip(H,G),spar[keep].T)),initial=0)<=3

# Literally signed superset of all old owners; restore unretired equator
# head once whenever a charged owner must be deleted. This overcharges
# rather than undercharges old antecedent coordinates already retired.
perm=np.array(w['permutation'],dtype='i8')
assert sorted(perm.tolist())==list(range(420))
scores=np.array([np.max(np.abs(ip(z,D.T)),axis=0)for z in groups])
V=U['base'];E=U['extra_points'];T=U['transforms'];O=U['extra_offsets']
M=np.unpackbits(U['retained_bits'],axis=1)[:,:len(V)].astype(bool)
d0={q['class_index']:q['add']for q in d['changes']}
L=ip(G,check.T); WGt=ip(G,A0.T)
cost=0;candidate_count=0
for j,Hj in enumerate(T):
 old=np.vstack((V[M[j]],E[O[j]:O[j+1]]))
 x=np.vstack((old,-old,np.array(d0.get(j,[]),dtype='i8').reshape(-1,48),Y[j]))
 h=ip(Hj,L);p=ip(Hj,WGt)
 assert np.max(np.abs(ip(x,h)),initial=0)<=3
 q=np.abs(ip(x,p)); assert np.max(q,initial=0)<=6
 cost+=int(np.any(q+scores[:,perm[j]][None,:]>8,axis=1).sum())
 candidate_count+=len(x)
assert d['bound']==58631571
body=d['bound']+g['owner_gain']+f['gain']-1
pre=body+4096+30*int(keep.sum())
assert body==58658586
val=pre+2*len(Z)-cost
print(json.dumps(dict(status='PASS',dimension=60,old_body=body,control_cube=4096,
 retained_sparse_heads=int(keep.sum()),frame_axes=48,active_axes=32,tails=len(Z),
 signed_candidates_tested=candidate_count,charged_blockers=cost,lower_bound=val)))
