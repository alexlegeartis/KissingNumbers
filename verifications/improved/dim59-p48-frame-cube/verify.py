"""Exact append of a frame-coloured mixed tail layer to the 59-dimensional P48p lift."""
from pathlib import Path
from hashlib import sha256
import argparse,gzip,json,numpy as np,io
assert __debug__

P=Path(__file__).resolve().parent
a=argparse.ArgumentParser(description=__doc__)
a.add_argument('--siblings',type=Path,default=P.parent)
a.add_argument('--body',type=Path)
a.add_argument('--witness',type=Path,default=P/'witness.json')
o=a.parse_args()
A=o.siblings/'dim49-60-p48-pointwise-refill'
B=o.siblings/'dim49-60-p48-signed-refill'
C=o.body or o.siblings/'dim50-60-p48-global-owner-refill'
def pin(p,s):
 z=p.read_bytes();assert sha256(z).hexdigest()==s,(str(p),sha256(z).hexdigest())
 return z
def read(p,s):return json.loads(gzip.decompress(pin(p,s)))
G=np.load(io.BytesIO(pin(A/'source/p48p_gram.npy','01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1')),allow_pickle=False)
d=np.load(io.BytesIO(pin(A/'heads59.npz','e92d40b8e36bbfdd81b744fa73692ddbf96cb099f0ce91c6475d8986b9024ab1')),allow_pickle=False)
p=np.load(io.BytesIO(pin(A/'source/capdirs_11.npz','62e21350d0f62718a16a012b8932caf9c5cb16db4f2bd00e99c8f5934e18ba62')),allow_pickle=False)
D=read(B/'p48_signed_delta59.json.gz','cd0c5b368828c13117fa845b97b12b8a930333a86893cb8bb5c0ef4d2446b708')
g=read(C/'global59.json.gz','3e4c865907694dea6317fa66386539f6fdd305ac9c01e0bc2d4811cef9a349c3')
f=read(C/'refill59.json.gz','28ee556c901fee014fb7e3a69fbbc2877c4558ad5edf6e0eb0840024d5b7bcfd')
bank=np.array(read(C/'bank.json.gz','65d79e04fe0164c344934ed8231b51faf7d316da101093235781cd129d8b3e9f'),dtype='i8')
tet=json.loads(pin(C/'tail55/witness.json','cc4151354d9b6eaf6551275a2d9424eb5b8eba11e0a4c3bf5db0bec31e9aa053'))['59']
assert D['bound']==56952414 and g['owner_gain']==6493 and f['gain']==8440
base=D['bound']+g['owner_gain']+f['gain']+224
assert base==56967571
assert np.array_equal(d['gram'],G) and bank.shape==(4000,48)
Y={q['class_index']:bank[q['insert']]for q in g['changes']}
for q in f['changes']:
 j=q['class_index'];Y[j]=np.r_[Y[j],np.array(q['extra'],dtype='i8').reshape(-1,48)]
D0={q['class_index']:q['add']for q in D['changes']}

W=json.loads(o.witness.read_text());X=np.array(W['heads'],dtype='i8');Q=np.array(W['weighted_tails'],dtype='i8')
assert len(X)>=16

def num(q):return int(abs(q).max(initial=0))
def mat(a,b):
 assert a.shape[-1]==b.shape[0] and a.shape[-1]*num(a)*num(b)<2**62
 return a.astype('i8')@b.astype('i8')
def fast(a,b):
 # Exact products and all intermediate integer sums represented by binary64.
 assert a.shape[-1]==b.shape[0] and a.shape[-1]*num(a)*num(b)<2**51
 return a.astype('f8')@b.astype('f8')
assert X.ndim==2 and X.shape[1]==48
assert np.array_equal(mat(mat(X,G),X.T),6*np.eye(len(X),dtype='i8'))
assert Q.ndim==2 and Q.shape[1]==12
I=Q[:,0];Z=Q[:,1:];assert I.min()>=0 and I.max()<len(X) and set(I)==set(range(16,len(X)))
scale=int(W['scale']);assert scale==10000 and num(Z)<=scale//2
assert len(set(map(tuple,Q)))==len(Q)
# Retain all tetrahedral lift points from the pinned existing code.
oldH=np.array(tet['heads'],dtype='i8').reshape(-1,48)
assert oldH.shape==(32,48) and num(mat(mat(X,G),oldH.T))<=2
assert len(tet['supports'])==8 and all(len(set(a))==3 for a in tet['supports'])
for a in tet['supports']:
 assert 100*int(np.abs(Z[:,a]).sum(1).max())<=83*scale
# Old head factor ≤1/(6 sqrt(2)) <118/1000 from new head square ≤1/3;
# old triple factor sqrt(5/24)<457/1000, total <= .49731 <1/2.
assert 72*118**2>1000**2 and 24*457**2>5*1000**2
assert 11800+83*457<50000
norm=(Z*Z).sum(1);square=scale**2
assert np.all(3*norm>=2*square) and np.all(4*norm<=3*square)
assert all(v==0 for v in p['c']-1) and int(p['norm'][0])==36 and int(p['R'][0])==2
AA=mat(Z,p['A'].T);BB=mat(Z,p['B'].T)
def nonneg(a,b):
 if a>=0 and b>=0:return True
 if a<0 and b<0:return False
 if a>=0:return a*a>=2*b*b
 return a*a<=2*b*b
for a,b in zip(AA.flat,BB.flat):
 a,b=int(a),int(b)
 assert nonneg(3*scale-a,-b) and nonneg(3*scale+a,b)
# Largest head product for equal axis is when the two head signs agree.
F=mat(Z,Z.T)
u,v=np.triu_indices(len(Z),1);r=square-2*F[u,v]
assert np.min(r,initial=0)>=0 and num(r)**2<2**62 and 4*square**2<2**62
assert np.all(r*r>=(I[u]==I[v])*4*(square-norm[u])*(square-norm[v]))
# A given tail with opposite head signs has scalar product 2*norm/square-1.
assert np.all(4*norm<=3*square)


# Weighted cube, two copies of every even word, at exact tail square 3/4.
from itertools import product
from functools import reduce
from operator import xor
from fractions import Fraction as Fr
signs=np.array([u for u in product((-1,1),repeat=11) if u.count(-1)%2==0],dtype='i8')
assert signs.shape==(1024,11)
weights=np.array([50]*8+[49]*3,dtype='i8')
cube=signs*weights;N=int(weights@weights);assert N==27203
C0=np.array([reduce(xor,np.flatnonzero(z<0),0)for z in signs])
assert set(C0)==set(range(16))
dot=mat(cube,cube.T)
k,l=np.triu_indices(len(cube),1);same=C0[k]==C0[l]
assert np.all(3*dot[k,l]<=(2-same.astype('i8'))*N)
# Each tail coordinate absolute value <= (sqrt(3)/(2 sqrt(N)))*50 < .263.
assert 4*N*263**2>3*50**2*1000**2
# Cube axes orthogonal to weighted six-support axes.
assert 263*int(np.abs(Z).sum(1).max())<500*scale
# Old tetra cap: head <=1/sqrt(96)<.103; three coordinates in old tail.
assert 96*103**2>1000**2 and 24*457**2>5*1000**2
assert 103000+3*263*457 < 500000
# Exact old direction cap |dot| <=1417/2500. The signed old pure poles
# with positive dot >1/2 are removed; all signed lift direction classes stay.
ca=mat(cube,p['A'].T);cb=mat(cube,p['B'].T)
dead=set()
for uu in range(len(cube)):
 for jj,(a,b) in enumerate(zip(ca[uu].tolist(),cb[uu].tolist(),strict=True)):
  U=a*a+2*b*b;V0=2*a*b
  assert nonneg(48*N*1417**2-2500**2*U,-2500**2*V0)
  if nonneg(a,b) and not nonneg(12*N-U,-V0):dead.add(jj)
# Old body head contribution sqrt(3)/8 < .216507 and old tail term
# <= one-half times the exact direction cap.
assert 64*216507**2>3*1000000**2
assert Fr(216507,1000000)+Fr(1417,5000)<Fr(1,2)
# Cube's same-tail opposite-sign pair dot is (3/4)-(1/4)=1/2.
# Cube's head amplitude 1/2 is safe at its own retained equatorial ray.

V=d['base'];E=d['extra_points'];T=d['transforms'];R=d['extra_offsets']
M=np.unpackbits(d['retained_bits'],axis=1)[:,:len(V)].astype(bool)
H=mat(G,X.T)
visited=0
for j,U in enumerate(T):
 c=mat(U,H)
 stock=(V[M[j]],E[R[j]:R[j+1]],np.array(D0.get(j,[]),dtype='i8').reshape(-1,48),Y.get(j,np.empty((0,48),dtype='i8')))
 for k,z in enumerate(stock):
  assert num(fast(z,c))<=3
  visited+=len(z)*(2 if k<2 else 1)
# New-vs-old lifts: head term <= sqrt((1/3)*(3/4))/2=1/4,
# since abs raw head dot ≤3; old tail amplitude 1/2 gives ≤1/4.
# Equator: all distinct norm-six lines have |raw dot|<=3, so safe.
# Remove both old equator rays for precisely those axes with a head squared amplitude above 1/4.
deleted=2*len(set(I[4*norm<3*square]))
n=2*(len(Z)+len(cube))
print(json.dumps(dict(status='PASS',dimension=59,prior=base,old_signed_superset=visited,axes=len(X),weighted_tails=len(Z),cube_words=len(cube),doubled_words=len(cube),new_vectors=n,removed_equator=deleted,removed_poles=len(dead),improvement=n-deleted-len(dead),lower_bound=base+n-deleted-len(dead))))

assert json.loads((P / 'claims.json').read_text())['lower_bounds']['59'] == base+n-deleted-len(dead)
