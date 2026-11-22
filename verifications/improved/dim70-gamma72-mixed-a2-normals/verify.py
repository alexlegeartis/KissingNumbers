"""Exact Gamma72 mixed-normal link replay; GCC C++17, NumPy, SymPy, Node."""
import hashlib,json,subprocess,sys,tempfile,os
from collections import defaultdict
from fractions import Fraction as Fq
from itertools import combinations_with_replacement
from math import lcm
from pathlib import Path
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent
PARENT=Path(os.getenv('GAMMA_PARENT',str(HERE.parent/'dim68-91-lattice-code-continuation/mixed_shell71_20261117/leech_structure.npz'))).resolve()
PIN='3253520e074ea6866bb39c22bf05b595945b6183b48d9ae145758c65db79752e'
assert hashlib.sha256(PARENT.read_bytes()).hexdigest()==PIN
DATA=np.load(PARENT)
F,A,B,T=[DATA[k].astype('i8')for k in ('gram','a','b','t')]
S=1<<np.arange(24,dtype='i8')
def key(X):return (X%2)@S
def line(z):return ' '.join(map(str,z))+'\n'

def source():
 G=np.loadtxt(HERE/'gram72.txt',dtype='i8');FF=np.kron(np.eye(3,dtype='i8'),F)
 assert np.array_equal(T@FF@T.T,2*G)
 U=s.Matrix(T);assert abs(int(U.det()))==2**36
 assert np.array_equal(A+B,np.eye(24,dtype='i8'))
 assert np.array_equal(A@B,2*np.eye(24,dtype='i8'))
 assert np.array_equal(A@F@A.T,2*F) and np.array_equal(B@F@B.T,2*F)
 for t in T.reshape(72,3,24):
  a,b,c=t
  assert key(a@B)==key(b@B)==key(c@B) and key((a+b+c)@A)==0
 # A,B give complementary binary maps of ranks twelve. The label equations
 # have 36 independent constraints, hence by the determinant check describe
 # exactly this Gamma72 lattice rather than merely a congruence over-lattice.
 au=s.Matrix(A);bu=s.Matrix(B)
 assert abs(int(au.det()))==abs(int(bu.det()))==2**12
 def rank(m):
  q=set()
  for w in m:
   z=int(key(w))
   while z:
    bit=z.bit_length()-1
    hit=next((v for v in q if v.bit_length()-1==bit),None)
    if hit is None:q.add(z);break
    z^=hit
  return len(q)
 assert rank(A)==rank(B)==12
 return FF

def ldl():
 G=s.Matrix(F);L,D=G.LDLdecomposition(hermitian=False)
 assert L*D*L.T==G and all(D[i,i]>0 for i in range(24))
 weights=[];cols=[]
 for i in range(24):
  r=[L[j,i]for j in range(i,24)];x=lcm(*(int(v.q)for v in r));cols.append([int(x*v)for v in r]);weights.append(D[i,i]/x**2)
 m=lcm(*(int(r.q)for r in weights)); assert m<2**96
 return line([m])+''.join(line([int(m*c),*row])for c,row in zip(weights,cols))

def minlink(work,r):
 inp=ldl();exe=work/'enum'
 subprocess.run(['g++','-O3','-std=c++17',str(HERE/'enumerate.cpp'),'-o',str(exe)],check=True)
 def make(g,radius):
  f=work/'leech.bin'
  subprocess.run([str(exe),str(g),str(f),str(radius)],check=True,text=True,input=inp)
  return np.fromfile(f,dtype='<i4').reshape(-1,25).astype('i8')
 z=make(-1,4);assert len(z)==196560 and (z[:,0]==4).all()
 W=z[:,1:];V=W@F
 assert np.all((W*V).sum(1)==4) and len(np.unique(W,axis=0))==len(W)
 g=int(key(r@B)); Z=make(g,16);A0=Z[:,1:]@A//2;K0=key(A0)
 assert np.array_equal(A0*2,Z[:,1:]@A)
 assert np.array_equal((A0@F*A0).sum(1)*2,Z[:,0]) and np.all(key(A0@B)==g)
 blocks=defaultdict(list);R=[]
 for i,k in enumerate(K0):blocks[Z[i,0]//2,k].append(i)
 for i,a in enumerate(A0):
  for j in blocks.get((12-Z[i,0]//2,K0[i]^key(r)^g),[]):
   R.append(np.r_[a,A0[j],r])
 assert len(R)==12800
 p=V@r;zero=np.zeros(24,dtype='i8')
 for i in np.flatnonzero(p==2):R.append(np.r_[zero,zero,2*W[i]])
 assert len(R)==17400
 O=W[p==0]+r; G0=key(O@B); K1=key(O)
 W8=W@A
 U8=defaultdict(list)
 for i,k in enumerate(key(W8)):U8[k].append(i)
 for o in np.flatnonzero(G0==0):
  for i in U8.get(K1[o],[]):R.extend((np.r_[zero,W8[i],O[o]],np.r_[W8[i],zero,O[o]]))
 group=key(W@B);K=key(W)
 indices=np.argsort(group);U,st,ns=np.unique(group[indices],return_index=True,return_counts=True)
 ids={k:indices[p:p+n]for k,p,n in zip(U,st,ns)}
 oi=np.argsort(G0);OU,Os,On=np.unique(G0[oi],return_index=True,return_counts=True)
 for k,e,f in zip(OU,Os,On):
  if not k:continue
  p=ids[k];cand=(K[p,None]^K[p][None,:])^k
  mapping=defaultdict(list)
  for ii,jj in np.ndindex(cand.shape):mapping[cand[ii,jj]].append((p[ii],p[jj]))
  for o in oi[e:e+f]:
   for i,j in mapping.get(K1[o],[]):R.append(np.r_[W[i],W[j],O[o]])
 R=np.asarray(R,dtype='i8').reshape(-1,3,24)
 return R

def link_population():
 # Norm-eight products with a norm-eight axis are 0, ±1, ±2, ±3, ±4, ±8.
 # The spherical-design moments through degree 8 suffice.
 from math import prod
 M=6218175600
 rhs=[s.Rational(M*prod(range(1,2*k,2))*64**k,prod(72+2*i for i in range(k)))-2*8**(2*k)for k in range(5)]
 V=s.Matrix([[j**(2*k)for j in range(5)]for k in range(5)])
 y=list(V.inv()*s.Matrix(rhs))
 assert all(v.is_Integer and v>=0 for v in y)
 return int(y[-1])//2

def replay():
 FF=source()
 N=json.loads((HERE/'normals.json').read_text());r=np.array(N['pure_root'],dtype='i8')
 m=np.array(N['mixed_normal'],dtype='i8').reshape(3,24)
 z=np.zeros(24,dtype='i8');a=np.r_[z,z,2*r];b=m.ravel()
 assert int(r@F@r)==4
 n=np.stack([a,b]);assert (n@FF@n.T).tolist()==[[16,8],[8,16]]
 assert np.count_nonzero(np.any(m!=0,axis=1))>=2
 with tempfile.TemporaryDirectory()as d:R=minlink(Path(d),r)
 # Validate every actual leaf with integral arithmetic. A complete list follows
 # from the exact shell population, regardless of completeness of enumeration.
 assert len(R)==link_population()==127800
 assert len(np.unique(R.reshape(-1,72),axis=0))==len(R)
 bs=[key(R[:,i]@B)for i in range(3)]
 assert (bs[0]==bs[1]).all() and (bs[1]==bs[2]).all()
 assert (key(R.sum(1)@A)==0).all()
 assert ((R@F*R).sum(axis=(1,2))==16).all() and (R[:,2]@F@r==4).all()
 assert ((R==m[None,:,:]).all(axis=(1,2))).sum()==1
 x=int(np.count_nonzero(R.reshape(-1,72)@(FF@b)==8))
 # Same three projection fibres; certify all six unordered pair types.
 def g(p,q):return Fq(2*p[0]*q[0]+2*p[1]*q[1]-p[0]*q[1]-p[1]*q[0],12)
 U=((-1,-1),(-1,0),(-2,-1))
 for p,q in combinations_with_replacement(U,2):
  a,b=8-g(p,p),8-g(q,q);c=4-g(p,q)
  assert a>0 and b>0 and (c<=0 or 4*c*c<=a*b)
 c=json.loads(subprocess.check_output(['node',str(HERE/'affine.cjs'),str(HERE/'affine.json'),str(x)],text=True))
 vals=list(map(int,c['counts']));L=sum(vals);old=1361142400
 assert vals==[610980300-36*x,610980300-36*x,139219696+48*x] and L==int(c['lower'])
 out={'dimension':70,'previous':old,'lower_bound':L,'gain':L-old,'corner':x,'link_count':len(R),'Fibre_populations':vals}
 assert out==json.loads((HERE/'results.json').read_text()) and L>old
 print(json.dumps(out))
if __name__=='__main__':
 assert __debug__
 replay()
