"""Exact replay of Γ72 depth-six fibres using seven automorphic integral-block charts."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations, product
import hashlib,sys,runpy,json,numpy as np
from sympy import Matrix,eye,zeros,diag
assert __debug__
here=Path(__file__).resolve().parent
base=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else here.parent/'dim73-82-gamma72-depth6-lattice'
pins={
'verify.py':'35b943ae846ba62a0270f1b3d5ff4f5280c208de4f4c2c5dd5df24ce3bc787df',
'theta.py':'168d5d052ce0a6e62085800a1fa956e88791519349bab7309409dec85420f073',
'certificate.json':'e5e656d728ef6ca1a7242e2863c8c41447247a9165a7b5aa862ea15d11a1584e',
'results.json':'f4042a060aad10f6a3766aff8347e05fd5a9541a517516d3af3f00e49d1b5ccf',
'comparisons.json':'d04816fd2d1cd4fa1efcc880245a3fcd104d0254c74dfaebe2b1dcf7cdeb7e4d',
'catalog_gram.json':'1791fab8c152302db2652ec66f227eb473ce4007f934431ef6bbb2ab2f2d87ab'}
assert all(hashlib.sha256((base/f).read_bytes()).hexdigest()==h for f,h in pins.items())
sys.path.insert(0,str(base)); P=runpy.run_path(str(base/'verify.py'))
G,T,F,A,B=(Matrix(P[k].tolist()) for k in ('G','T','F','A','B'))
U=Matrix(P['c']['leech_basis_change'])
Ti=T.inv(method='DM')
I=eye(24); Z=zeros(24)
FA=U*G[:24,48:]*(G[:24,:24]/2).inv()*U.inv()
FB=I-FA
assert FA*F==F*FB.T and FA*FB==2*I and FA*F*FA.T==2*F
R=eye(72)-Matrix.vstack(I,I,FB)*Matrix.hstack(I,I,FA)/2
C=zeros(72)
for i in range(72):C[i,(i+24)%72]=1
S=diag(-I,I,I)
X=Matrix.vstack(Matrix.hstack(Z,I,Z),Matrix.hstack(I,Z,Z),Matrix.hstack(Z,Z,I))
def pack(v):return sum((int(x)&1)<<i for i,x in enumerate(v))
def rows(v):return tuple(pack(r)for r in v.tolist())
def apply(m,v):
 w=0
 while v:
  t=v&-v;v-=t;w^=m[t.bit_length()-1]
 return w
def echelon(V):
 R={}
 for u in V:
  while u:
   i=u.bit_length()-1
   if i in R:u^=R[i]
   else:R[i]=u;break
 for i in sorted(R):
  for j in R:
   if j>i and R[j]>>i&1:R[j]^=R[i]
 return tuple(R[i]for i in sorted(R))
generators=[]
for m in (R,C,S,X):
 L=T*m*Ti
 assert all(z.q==1 for z in L) and L*G*L.T==G
 generators.append(rows(L))
I0=echelon(rows(2*Ti));assert len(I0)==36
frames=[I0];words=[()]
E={I0}
for i,v in enumerate(frames):
 for j,m in enumerate(generators):
  f=echelon(apply(m,x)for x in v)
  if f not in E:
   E.add(f);frames.append(f);words.append(words[i]+(j,))
 assert len(E)<=7
assert len(frames)==7
D=json.loads((here/'witness.json').read_text())
a,b=D['triangle']
B0={0,a,b,a^b};assert len(B0)==4
raw_of=lambda a,b:[(w>>i)&1 for w in (a,b,0)for i in range(24)]
legal={}
for u,v in product(B0,repeat=2):
 raw=raw_of(u,v)
 assert P['syndrome'](raw)==(u,v,0)
 if u and v and u!=v:
  P['prove'](raw,6)
  lower=6
 elif u or v:
  # There is a one-block representation with these singular labels.
  # To beat squared depth four, it must have Leech norm four or six.
  x=u or v
  assert 0 not in P['M4'][x]
  assert P['orth'](x,0)
  lower=4
 else:lower=0
 v0=2*Matrix([raw])*Ti
 assert all(x.q==1 for x in v0)
 phase=pack(list(v0))
 for w in words:
  z=phase
  for j in w:z=apply(generators[j],z)
  assert z not in legal or legal[z]==lower
  legal[z]=lower
assert len(legal)==64 and sum(x==4 for x in legal.values())==21 and sum(x==6 for x in legal.values())==42
EE=[]
for i,j in combinations(range(8),2):
 for p,q in product((-2,2),repeat=2):
  z=[0]*8;z[i]=p;z[j]=q;EE.append(z)
for r in product((-1,1),repeat=8):
 if sum(x<0 for x in r)%2==0:EE.append(list(r))
assert len(EE)==240
def scope(z,typ):
 if typ=='E7':return z[0]==z[1]
 if typ=='E6':return z[0]==z[1]==z[2]
 if typ=='D6':return z[6]==z[7]==0
 if typ=='D5':return z[5]==z[6]==z[7]==0
 assert typ=='E8';return True
dim={'E8':8,'E7':7,'E6':6,'D6':6,'D5':5}
prev=json.loads((here/'claims.json').read_text())
pre=here/'predecessor.json'
assert hashlib.sha256(pre.read_bytes()).hexdigest()=='70492ff918b4e38416be11a53a471ec05cba1b5c3f0003dfd0bb77538da7c124'
pre=json.loads(pre.read_text())
assert set(prev)=={str(72+dim[row['vertical']])for row in D['rows']}
for row in D['rows']:
 typ=row['vertical'];d=72+dim[typ]
 E=[z for z in EE if scope(z,typ)]
 assert all(sum(a*a for a in z)==8 for z in E)
 assert Matrix(E).rank()==dim[typ]
 assert all(sum(a*b for a,b in zip(x,y))<=4 for x,y in combinations(E,2))
 W=row['rays']; seen=set()
 for phase,z in W:
  assert z in E and tuple(z) not in seen and tuple(-a for a in z)not in seen
  seen.add(tuple(z));assert legal.get(phase,-1)==6
 for (p,x),(q,y) in combinations(W,2):
  t=abs(sum(a*b for a,b in zip(x,y)))
  assert t in (0,4)
  assert legal.get(p^q,-1)>=4+t//2
 # Both signs of every ray use the same half phase. Their tail separation
 # is sqrt(8); other rays have squared tail separation >=4-2|x.y|/4.
 # All pure vertical points coexist with every tilted fibre: 2*sqrt(2)=sqrt(8).
 old=pre[str(d)]['lower_bound']
 count=6218175600+2*len(W)*136704+len(E)
 assert prev[str(d)]=={'prior':old,'gain':count-old,'lower_bound':count}
 assert count>old
 print('PASS',d,count,len(W),'fibre pairs',len(E),'pure roots')
