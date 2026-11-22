"""Independent direct exact enumeration of every changed tail-module pair."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from functools import reduce
import numpy as np,json,argparse,hashlib
P=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--siblings',type=Path,default=P.parent)
args=p.parse_args()
Gfile=args.siblings/'dim49-60-p48-pointwise-refill/source/p48p_gram.npy'
assert hashlib.sha256(Gfile.read_bytes()).hexdigest()=='01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'
G=np.load(Gfile).astype('i8')
W=json.loads((P/'witness.json').read_text())
H=np.asarray(W['frame'],dtype='i8')
S=np.asarray(W['sparse'],dtype='i8')
U={x['color']:x['tetrahedra']for x in W['upgrade']}
C=set(W['subcode'])
zero=np.zeros(48,dtype='i8')
def bits(w):return np.array([1-2*((w>>j)&1) for j in range(12)],dtype='i8')
def syn(w):return reduce(int.__xor__,[i for i in range(12)if w>>i&1],0)
A={}
def add(alpha,k,h,z):
 a=A.setdefault((alpha,k),[])
 a.append((h,z))
B={}
for w in range(4096):
 if not w.bit_count()%2:B.setdefault(syn(w),set()).add(w)
for s,V in B.items():
 if s in U:
  colors=[]
  while V:
   w=min(V);block={w^v for v in C};assert block<=V
   colors.append(block);V-=block
  assert len(colors)==8
  for pal,ts in zip(U[s],colors):
   if pal is None:
    for w in ts:
     for h in (H[s],-H[s]):add(F(1,4),12,h,bits(w))
   else:
    for w in ts:
     for h in pal:add(F(2,5),12,h,bits(w))
 else:
  for w in V:
   for h in (H[s],-H[s]):add(F(1,4),12,h,bits(w))
for off in (0,6):
 for mask in range(8):
  choices=[off+2*j+((mask>>j)&1)for j in range(3)]
  for signs in product((-1,1),repeat=3):
   v=np.zeros(12,dtype='i8');v[choices]=signs
   color=min(mask,mask^7)
   pol=int(np.prod(signs))
   for h in pol*S[color]:add(F(3,8),3,h,v)
# Reconstruct the literal old 840 tail directions independently by edge colors.
palette=[]
for j in range(12):
 for e in (-2,2):
  z=np.zeros(12,dtype='i8');z[j]=e;palette.append(z)
supports=[]
for off in (0,6):
 for i in (0,2,4):
  supports.append([off+j for j in range(6) if j not in (i,i+1)])
def tag(i,j):
 return i if j==5 else j if i==5 else (3*(i+j))%5
for i,j in product(range(6),repeat=2):
 if i>=j:continue
 for a,b in product(range(6),repeat=2):
  if a<b and tag(i,j)==tag(a,b):supports.append([i,j,a+6,b+6])
assert len(supports)==51
for s in supports:
 for signs in product((-1,1),repeat=4):
  z=np.zeros(12,dtype='i8');z[s]=signs;palette.append(z)
assert len(palette)==840
# Compare two signed radicals using two squarings in rational arithmetic.
def le(h,a,z,b):
 x=h*h*a;y=z*z*b; c=F(1,4)
 if h>=0 and z>=0:return x+y<=c and 4*x*y<=(c-x-y)**2
 if h<=0 and z<=0:return True
 if h<0:x,y=y,x
 return x<=c+y or (x-c-y)**2<=y
T=[(k,np.array([x[0]for x in v],dtype='i8'),np.array([x[1]for x in v],dtype='i8'))for k,v in A.items()]
count=0
for i,((a,k),X,V) in enumerate(T):
 assert np.all(((X@G)*X).sum(1)==6) and np.all((V*V).sum(1)==k)
 for j,((b,l),Y,Z) in enumerate(T[:i+1]):
  s=a*b/36;t=(1-a)*(1-b)/(k*l)
  table=np.asarray([[le(h,s,w,t)for w in range(-12,13)]for h in range(-6,7)])
  for d in range(0,len(V),128):
   e=min(d+128,len(V));h=X[d:e]@G@Y.T;z=V[d:e]@Z.T
   assert np.all(abs(h)<=6) and np.all(abs(z)<=12)
   ok=table[h+6,z+12]
   if i==j:
    keep=np.arange(d,e)[:,None]>np.arange(len(Y))
    assert np.all(ok[keep]);count+=sum(map(sum,keep))
   else:
    assert np.all(ok);count+=ok.size
 z=V@np.asarray(palette,dtype='i8').T
 assert np.all(z*z*(1-a).numerator<=k*(1-a).denominator)
 count+=z.size
print('Exact direct new/new and new/pole pairs:',count,'new points:',sum(len(X)for _,X,_ in T))
