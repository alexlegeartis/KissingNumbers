# Independent physical checks for literal new layers and conservative owner charge.
import gzip,json,numpy as np
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
H0=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--siblings',type=Path,default=H0.parent)
ap.add_argument('--body',type=Path)
ap.add_argument('--refill',type=Path)
ap.add_argument('--witness',type=Path,default=H0/'witness.json')
args=ap.parse_args()
P=args.siblings
Q=args.body or P/'dim50-60-p48-global-owner-refill'
refill=args.refill or Q/'refill60.json.gz'
def gz(p):return json.loads(gzip.decompress(p.read_bytes()))
S=np.load(P/'dim49-60-p48-pointwise-refill'/'heads60.npz',allow_pickle=False)
G=S['gram'].astype('i8');T=S['transforms']; B=S['base'];E=S['extra_points']; O=S['extra_offsets']
M=np.unpackbits(S['retained_bits'],axis=1)[:,:len(B)].astype(bool)
old=Q/'tail60/witness.json'
old=json.loads(old.read_text())
w=json.loads(args.witness.read_text())
W=np.array(w['frame'],dtype='i8')
A=W[w['carriers']]
I,J=np.array(w['controls']).reshape(-1,2).T
H=np.r_[(W[I]+W[J])//2,(W[I]-W[J])//2]
assert np.array_equal(np.r_[A,H]@G@np.r_[A,H].T,np.diag([12]*32+[6]*16))
Z=np.array(w['tail_groups'],dtype='i8').reshape(-1,12)
Ac=np.repeat(A,48,axis=0)
assert len(Z)==1536 and len(set(map(tuple,Z)))==1536
assert np.all((Z*Z).sum(1)==10)
for ss in [np.arange(1536)[i:i+64]for i in range(0,1536,64)]:
 c=abs(Ac[ss]@G@Ac.T)+2*Z[ss]@Z.T
 c[np.arange(len(ss)),ss]=0
 assert c.max()<=16
# Cube signs: literal parity, color and all geometry including both head signs.
eps=np.array([[1-2*((x>>i)&1)for i in range(12)]for x in range(4096)if x.bit_count()%2==0],dtype='i8')
colour=np.bitwise_xor.reduce((eps==-1)*np.arange(12),axis=1)
assert len(eps)==2048
for ss in [np.arange(2048)[i:i+64]for i in range(0,2048,64)]:
 c=4*(colour[ss,None]==colour[None,:])+eps[ss]@eps.T
 c[np.arange(len(ss)),ss]=0
 assert c.max()<=8
assert abs(Z@eps.T).max()<=8
# Direct reconstruction of all kept sparse heads/tail signs. Rational
# positive-radical separator: a+b≤1/2 iff gap^2≥4a²b² and gap≥0.
def test(x,y):
 gap=F(1,4)-x-y
 return gap>=0 and gap*gap>=4*x*y
s=np.array(old['sparse'],dtype='i8')
Sx=[];Sy=[]
for a in range(4):
 for b in range(4):
  if not w['sparse_keep'][a][b]:continue
  for j,z in product((0,6),(a,a^7)):
   pos=[j+2*k+((z>>k)&1)for k in range(3)]
   for hi in product((-1,1),repeat=3):
    p=[0]*12
    for k,v in zip(pos,hi):p[k]=v
    Sx.append(np.prod(hi)*s[a,b]);Sy.append(p)
SX=np.array(Sx,dtype='i8');SY=np.array(Sy,dtype='i8')
assert len(SX)==512
QQ=abs(Ac@G@SX.T);RR=abs(Z@SY.T)
for q,r in set(zip(QQ.ravel(),RR.ravel())):
 assert test(F(int(q)**2,512),F(5*int(r)**2,384))
assert np.max((3*SX@G@SX.T+10*SY@SY.T)[np.triu_indices(512,1)])<=24
assert np.max(abs(H@G@SX.T))<=3 and np.max(abs(SY@eps.T))<=3
print('Literal frame, all 1536 q10 tails, signed cube and 512 retained sparse rays PASS',flush=True)
# Independent reconstruct all signed owner candidates and score exact scaled angles.
D=[]
for i in range(12):
 z=[0]*12;z[i]=2;D.append(z)
sp=[]
for j in (0,6):
 for p in (0,2,4):sp.append(tuple(j+i for i in range(6)if i not in (p,p+1)))
for i in range(5):
 P0=[(5,i),((i+1)%5,(i-1)%5),((i+2)%5,(i-2)%5)]
 for a,b in product(P0,repeat=2):sp.append(a+tuple(6+j for j in b))
for a in sp:
 for h in product((-1,1),repeat=4):
  v=np.zeros(12,dtype='i8');v[list(a)]=h
  if v[np.flatnonzero(v)[0]]>0:D.append(v)
D=np.array(D)
assert len(D)==420
Rmax=abs(Z@D.T).reshape(32,48,420).max(1)
assert Rmax.max()<=4
g=gz(Q/'global60.json.gz');f=gz(refill)
bank=np.array(gz(Q/'bank.json.gz'),dtype='i8')
y={v['class_index']:bank[v['insert']]for v in g['changes']}
for v in f['changes']:
 i=v['class_index'];y[i]=np.r_[y[i],v['extra']]
de=gz(P/'dim49-60-p48-signed-refill'/'p48_signed_delta60.json.gz')
p={v['class_index']:v['add']for v in de['changes']}
J=G@np.r_[A,H,s.reshape(-1,48)].T
seen=cost=0
for i,U in enumerate(T):
 k=Rmax[:,w['permutation'][i]]
 a=np.flatnonzero(k>2)
 old=np.r_[B[M[i]],E[O[i]:O[i+1]]].astype('i8')
 S0=np.r_[old,-old,np.array(p.get(i,[]),dtype='i8').reshape(-1,48),y[i]]
 # Work in float64 only for bounded integer products, executed exactly:
 # every intermediate sum is strictly below 2^52 in magnitude.
 J0=U@J
 assert 48*abs(S0).max()*abs(J0).max()<2**52
 z=abs(S0.astype('d')@J0.astype('d'))
 assert z[:,32:].max()<=3 and z[:,:32].max()<=6
 cost+=np.count_nonzero((z[:,a]+k[a]>8).any(1));seen+=len(z)
assert seen==6289330 and cost==236
bound=58658586+4096+30*16+2*len(Z)-cost
assert bound==58665998
print('Independent exact bounded-integer signed replay:',seen,cost,bound,flush=True)
