from pathlib import Path
from itertools import combinations
import json,re
import numpy as np
from sympy import Matrix
HERE=Path(__file__).resolve().parent
c=json.loads((HERE/'certificate.json').read_text())
import hashlib
from itertools import product
from math import isqrt
import theta
S=HERE/'catalog_gram.json'
assert hashlib.sha256(S.read_bytes()).hexdigest()=='1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab'
J=np.array(c['leech_frame'],dtype='i8')
n=72
vals=json.loads(S.read_text())
assert len(vals)==n*n
GG=Matrix(n,n,vals)
assert GG.det(method='domain-ge')==1 and GG.T==GG and all(GG[i,i]%2==0 for i in range(72))
G=np.array(GG.tolist(),dtype='i8')
def exact(x):
    assert all(z.denominator==1 for row in x.tolist() for z in row)
    return Matrix([[int(z) for z in row] for row in x.tolist()])
U=Matrix(c['leech_basis_change'])
Ui=exact(U.inv(method="DM")); assert abs(U.det())==1
I=Matrix([[int(i==j) for j in range(24)] for i in range(24)])
F0=exact(Matrix(G[:24,:24].tolist())/2)
A0=exact(Matrix(G[:24,48:].tolist())*F0.inv(method="DM"))
B0=I-A0
F=np.array((U*F0*U.transpose()).tolist(),dtype='i8')
A=np.array((U*A0*Ui).tolist(),dtype='i8')
B=np.eye(24,dtype='i8')-A
O=np.zeros((24,24),'i8')
ii=np.array(Ui.tolist(),dtype='i8')
t=np.array((A0*Ui).tolist(),dtype='i8')
b=np.array((B0*Ui).tolist(),dtype='i8')
T=np.block([[ii,ii,t],[O,b,b],[O,O,2*ii]])
Ti=Matrix(T.tolist()).inv(method="DM")
assert all(a.denominator<=2 for row in Ti.tolist() for a in row)
# The Leech block F is as in Nebe's exact construction (first G block divided by 2).
# Index of T is 2^36 because T diag(F,F,F) T^T = 2 G and both G and F are unimodular.
blockF=np.kron(np.eye(3,dtype='i8'),F)
assert np.array_equal(T@blockF@T.T,2*G)
assert Matrix(F.tolist()).det()==1
assert np.array_equal(J@F@J.T,8*np.eye(24,dtype='i8'))
assert all(np.max(abs(x))<10**4 for x in [J,G,F,T,A,B])
A%=2; B%=2
I=np.eye(24,dtype='i8')
assert np.array_equal((A+B)%2,I) and np.array_equal(A@A%2,A)
assert (A@F@A.T%2==0).all() and (B@F@B.T%2==0).all()
assert ((A@F@A.T).diagonal()%4==0).all() and ((B@F@B.T).diagonal()%4==0).all()
p24=1<<np.arange(24)
def wpack(v):return int(v%2@p24)
def rank(words):
 local={}
 for x in words:
  while x:
   i=x.bit_length()-1
   if i not in local:local[i]=x;break
   x^=local[i]
 return len(local)
assert rank(map(wpack,A))==rank(map(wpack,B))==12
for x,y,z in T.reshape(72,3,24):
 assert np.array_equal(x@A%2,y@A%2) and np.array_equal(x@A%2,z@A%2)
 assert ((x+y+z)@B%2==0).all()
# These three independent rank-12 constraints characterize the entire T image.
# Enumerate the complete norm-four Leech shell (integer vectors in a frame, sum of squares 32).
fwd=F@J.T
p2=np.unique(fwd%2,axis=0)
assert [r.tolist() for r in p2]==[[0]*24,[1]*24]
def pack4(r): return (wpack(r),wpack(r>>1))
def add4(a,b): return (a[0]^b[0],a[1]^b[1]^(a[0]&b[0]))
R={(0,0)}
for r in fwd:
 a=pack4(r%4); a2=add4(a,a)
 R={add4(t,x) for t in R for x in [(0,0),a,a2,add4(a2,a)]}
assert len(R)==8192 and all(lo or hi.bit_count()!=4 for lo,hi in R)
rows=[]
for i,j in combinations(range(24),2):
 for x,y in [(4,4),(4,-4),(-4,4),(-4,-4)]:
  q=np.zeros(24,dtype='i8');q[i]=x;q[j]=y;rows.append(q)
for lo,hi in R:
 if lo==0 and hi.bit_count()==8:
  idx=[i for i in range(24) if hi>>i&1]
  q=np.zeros(24,dtype='i8')
  for z in range(256):
   q[idx]=[-2 if z>>k&1 else 2 for k in range(8)]
   rows.append(q.copy())
 elif lo==(1<<24)-1:
  q=np.array([-1 if hi>>i&1 else 1 for i in range(24)],dtype='i8')
  for j in range(24):
   v=q.copy();v[j]*=-3;rows.append(v)
v=np.array(rows)@J
W=v[(v%8==0).all(axis=1)]//8
assert len(W)==len(np.unique(W,axis=0))==196560 and ((((W@F)*W).sum(1))==4).all()
M4={}
for a,b in zip((W@A%2)@p24,(W@B%2)@p24):
 M4.setdefault(int(a),set()).add(int(b))
assert len(M4)==4095 and 0 not in M4 and set(map(len,M4.values()))=={24}
def bits(i):return (i>>np.arange(24))&1
def orth(a,b):return int(bits(a)@F@bits(b))%2==0
def syndrome(r):
 x,y,z=np.array(r,dtype='i8').reshape(3,24)
 return wpack((x+z)@A),wpack((y+z)@A),wpack((x+y+z)@B)
def prove(r,cut):
 a,b,t=syndrome(r)
 if cut==4:
  assert a==b==0 and t!=0
  # A=0 has no norm-four label, and B is singular, so no single norm-six either.
  return
 assert cut==6 and a and b and a!=b
 U0=(a,b,a^b)
 assert all(orth(s,t) for s in U0)
 assert all(orth(s,x) for j in U0 for x in M4[j] for s in U0)
 # For two occupied blocks their two A labels are any of (a,b),(a,a+b),(b,a+b).
 assert all(t!=x^y for a0,b0 in combinations(U0,2) for x in M4[a0] for y in M4[b0])
 # This excludes 4+4. For 4+6, the norm-four B label and the total B label
 # annihilate U0, so the norm-six B label annihilates its own A label, impossible.
# Upper bounds from first few exact even unimodular Leech theta coefficients.
# Weight 12 form E4^3 - 720*Delta.
E4=theta.E4; De=theta.Del; Li=theta.power(E4,3)
S=[a-720*b for a,b in zip(Li,De)]
assert S[:5]==[1,0,196560,16773120,398034000]
# Every norm4 class has exactly two reps: two distinct norm4
# points x,y in the same parity imply x+y,x-y both norm>=16 unless x=±y.
assert sum(map(len,M4.values()))==98280
# Every norm6 class likewise has two reps and is anisotropic.
assert S[3]//2==2**23-2**11
# Norm8 reps can inhabit neither norm4 nor norm6 classes by the same
# parallelogram argument; in a fixed norm8 class different rays
# are orthogonal, hence at most 48 points. The exact shell count fills them all.
assert S[4]==48*(2**24-1-S[2]//2-S[3]//2)

def norm4(a,b):return b in M4.get(a,())
def norm8(a,b):return bool(a or b) and orth(a,b) and not norm4(a,b)
def contact(raw,w):
 a,b,t=syndrome(raw)
 one=[]
 # When the other two blocks are zero, the three possible placements are:
 if b==0:one.append((a,t))
 if a==0:one.append((b,t))
 if a==b:one.append((a,t))
 if w==2:return 2*sum(norm4(i,j) for i,j in one)
 if w==4:
  # For each placement with two occupied blocks distinct, pick first B and use total t for the other.
  pairs=((a,b),(a^b,b),(a^b,a))
  return 48*sum(norm8(i,j) for i,j in one)+4*sum(
    sum((s^t) in M4.get(j,()) for s in M4.get(i,())) for i,j in pairs)
 if w==6:
  prove(raw,6)
  return 136704
 if w==0:return int((a,b,t)==(0,0,0))
 raise AssertionError(w)

# The rank-four Cartan block in simple-root coordinates.
d4=np.array([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]],dtype='i8')
roots={tuple(x%2 for x in z) for z in product(range(-3,4),repeat=4) if int(z@d4@z)==2}
others=sorted(set(product(range(2),repeat=4))-roots-{(0,0,0,0)})
assert len(roots)==12 and len(others)==3
Bseed=next([a,b,*others[:2]]for a,b in combinations(sorted(roots),2)if Matrix([a,b,*others[:2]]).det()%2)
Q=np.array(Matrix(Bseed).inv_mod(2).tolist(),dtype='i8')
old=np.array(Q@c['raw_D4_seed']%2,dtype='i8')
new=np.array(Q@c['raw_D4_disjoint_plane_seed']%2,dtype='i8')
# Reverse the two occupied Leech block roles to obtain the opposite MRD orbit, keeping the same B kernel.
opp=old.reshape(4,3,24)[:,[1,0,2],:].reshape(4,72)
u1=np.array(c['raw73'],dtype='i8').reshape(3,24)
u1=u1[[1,0,2]]
u2=np.vstack([u1.ravel(),np.roll(u1,1,axis=0).ravel()])
big=np.vstack([old,opp,new])
def choices(v,cap=8):
  V=Matrix(v.tolist())
  assert V.is_positive_definite
  inv=V.inv()
  B=[isqrt(int(cap*inv[i,i]))for i in range(len(v))]
  U=[]
  for z in product(*(range(-b,b+1)for b in B)):
   z=np.array(z,dtype='i8')
   ss=int(z@v@z)
   if ss<=cap:U.append((ss,z))
  return U
def candidate(r):
 if r==1:
  return [np.array([[2]],dtype='i8')],np.array([c['raw73']],'i8')
 if r<=4:
  return [d4[:r,:r]],old[:r]
 if r<=6:
  return [d4,d4[:r-4,:r-4]],np.vstack([old,u2])[:r]
 return [d4[:min(4,r-i),:min(4,r-i)]for i in range(0,r,4)],big[:r]
def verify_short(raw,v):
 pa=syndrome(raw)
 if v==2:
  prove(raw,6)
 elif v==4:
  a,b,t=pa
  if a==b==0:prove(raw,4)
  elif a==b or not(a and b):assert t not in M4[a or b] and orth(a or b,t)
 elif v==6:assert any(pa)
 else:assert v==8
def build(r):
 V,raw=candidate(r)
 W=2*Matrix(raw.tolist())*Ti
 assert all(z.is_Integer for z in W)
 # All block-coordinate shifts lie in the integral Leech^3 ambient; thus u^2 is divisible by 4.
 assert (np.diag(2*raw@blockF@raw.T)%4==0).all()
 from collections import Counter
 D=Counter()
 for z in product(*(choices(v)for v in V)):
  s=sum(a[0]for a in z)
  if 0<s<=8:
   vz=np.concatenate([a[1]for a in z])
   D[s,tuple(vz@raw%2)]+=1
 by_v={2:0,4:0,6:0,8:0}
 verify_short_set=set()
 for (v,q),n in D.items():
  verify_short(q,v)
  by_v[v]+=n*contact(q,8-v)
 # Γ72's original norm-eight shell and all nonzero-tail layers are disjoint.
 return by_v,6218175600+sum(by_v.values())

report=[]
for r in range(1,11):
 by,total=build(r)
 row={'dimension':72+r,'count':total,'contacts_by_vertical_norm':by}
 report.append(row)
 print('PASS',row,flush=True)
import sys
if '--write' in sys.argv:
 (HERE/'results.json').write_text(json.dumps(report,indent=2)+'\n')
else:
 assert json.loads(json.dumps(report))==json.loads((HERE/'results.json').read_text())
prev=json.loads((HERE/'comparisons.json').read_text())
for previous,certified in zip(prev,report,strict=True):
 assert previous['dimension']==certified['dimension']
 assert previous['lower_bound']==certified['count']
 assert previous['gain_over_previous']==certified['count']-previous['previous_repository_bound']>0
print('PASS exact comparison arithmetic')
