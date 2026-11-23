"""Exact replay: partial minimum-depth-six Γ72 fibres in 80D."""
from pathlib import Path
from itertools import combinations,product
from functools import lru_cache
from hashlib import sha256
import sys,runpy,json
import numpy as np
from sympy import Matrix
assert __debug__
here=Path(__file__).resolve().parent
prior=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here.parent/'dim77-79-gamma72-seven-charts'
base=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else prior.parent/'dim73-82-gamma72-depth6-lattice'
pins = {
  "verify.py": "9743a247fbd10c9ec82ea1892a8ea6968fe92fa33ec2db2ad0188231ab5d7563",
  "witness.json": "0aad7a48df37bf825dcf734fa5b574602bbbb189169fb80e09b05659d0132a21",
  "predecessor.json": "70492ff918b4e38416be11a53a471ec05cba1b5c3f0003dfd0bb77538da7c124",
  "claims.json": "ec7d20f3211210f4fc6280ce18dd8f9e60d8c2c26fa91ec91885bd1a6b6d51ac"
}
assert all(sha256((prior/k).read_bytes()).hexdigest()==v for k,v in pins.items())
sys.argv=[str(prior/'verify.py'),str(base)]
R=runpy.run_path(sys.argv[0])
S=R['legal']
D0=R['D']
a,b=D0['triangle']
raw={j:{}for j in range(len(R['words']))}
for u,v in product((0,a,b,a^b),repeat=2):
 k=R['pack'](list(2*Matrix([R['raw_of'](u,v)])*R['Ti']))
 for j,w in enumerate(R['words']):
  h=k
  for i in w:
   h=R['apply'](R['generators'][i],h)
  assert h not in raw[j]
  raw[j][h]=(u,v,0)
assert len(raw)==7 and all(len(d)==16 for d in raw.values())

M=R['P']['M4']
P0=R['P']; G0=P0['F']; A0=P0['A']; B0=P0['B']; W0=P0['W']; pk=P0['p24']
R4={}
for a,b,row in zip((W0@A0%2)@pk,(W0@B0%2)@pk,W0):R4[int(a),int(b)]=row
assert len(R4)==98280
@lru_cache(None)
def norm8(z):
  # Class (z,0) is nonzero singular and not represented at norm4.
  assert 0 not in M[z] and z
  u,x=next((u,x)for u,S in M.items() if u!=z for x in S if (u^z,x)in R4)
  v=R4[u,x]+R4[u^z,x]
  assert int(v@G0@v)==8
  Z=v-2*W0[W0@G0@v==4]; Z=np.vstack([Z,v,-v])
  assert len(Z)==len(np.unique(Z,axis=0))==48
  assert np.all(((Z@G0)*Z).sum(axis=1)==8) and np.all((Z@A0%2)@pk==z) and np.all((Z@B0%2)@pk==0)
  return Z
def cutgood(V,delta,mask):
 a,b,t=delta;assert t==0 and (a or b) and (a==0 or b==0 or a==b)
 z=a or b
 for h in [0,z]:
  u=(a^h,b^h,h)
  ix=[i for i,r in enumerate(u)if r]
  if len(ix)==1:
   gg=norm8(z)@G0
   mask&=~(np.abs(V[:,ix[0],:]@gg.T)==4).any(axis=1)
  else:
   assert len(ix)==2
   j,k=ix
   gg=np.array([R4[z,x]for x in sorted(M[z])])@G0
   l=np.abs(V[:,j,:]@gg.T);r=np.abs(V[:,k,:]@gg.T)
   mask&=~((l+r==4)|(np.abs(l-r)==4)).any(axis=1)
 return mask
def make_triple(q):
 a,b,t=q
 return np.array([(R4[h^a,x],R4[h^b,y],R4[h,x^y^t])for h,H0 in M.items() if h not in (a,b)
  for x in M[h^a] for y in M[h^b] if x^y^t in H0],dtype='i8')


def bits3(P):
 Bc=[]
 k=R['P']['B']
 nb={}
 for rr in k:
  r=sum(int(x)<<i for i,x in enumerate(rr))
  while r:
   j=r.bit_length()-1
   if j in nb:r^=nb[j]
   else:nb[j]=r;break
 S=[0]
 for y in nb.values():S += [x^y for x in S]
 assert len(S)==4096
 return S
Bl=bits3(R['P'])
Bv=np.array(Bl,dtype='i8')
def orth(a,b):
 return R['P']['orth'](a,b)
@lru_cache(None)
def get_pair(a,b):
 for u,S in M.items():
  v=u^a
  if not v:continue
  for q in S:
   if (v,q^b) in R4:return R4[u,q],R4[v,q^b]
 raise AssertionError((a,b))
@lru_cache(None)
def norm_6(a,b):
 assert not orth(a,b)
 U,V=get_pair(a,b)
 d=int(U@G0@V)
 assert abs(d)==1
 w=U-V*d
 assert int(w@G0@w)==6 and int((w@A0%2)@pk)==a and int((w@B0%2)@pk)==b
 return w
@lru_cache(None)
def norm_8(a,b):
 assert (a or b) and orth(a,b) and b not in M.get(a,())
 U,V=get_pair(a,b)
 assert int(U@G0@V)==0
 n=U+V
 r=n-2*W0[W0@G0@n==4]; r=np.vstack([r,n,-n])
 assert len(r)==len(np.unique(r,axis=0))==48
 assert np.all(((r@G0)*r).sum(axis=1)==8)
 assert np.all((r@A0%2)@pk==a) and np.all((r@B0%2)@pk==b)
 cent=np.argmax(r!=0,axis=1)
 return r[r[np.arange(48),cent]>0]
def two_occupied(q):
 a,b,t=q
 ans=[]; zz=np.zeros(24,dtype='i8')
 for h in [0,a,b]:
  u=(h^a,h^b,h)
  I=[j for j,c in enumerate(u)if c]
  assert len(I)==2
  for j,k in (I,I[::-1]):
   for q0 in M[u[j]]:
    val=[zz,zz,zz];val[j]=R4[u[j],q0]
    for y in norm_8(u[k],q0^t):
     val[k]=y;ans.append(tuple(val))
  j,k=I
  for q0 in Bl:
   if not orth(u[j],q0) and not orth(u[k],q0^t):
    val=[zz,zz,zz];val[j]=norm_6(u[j],q0);val[k]=norm_6(u[k],q0^t)
    ans.append(tuple(val))
 assert len(ans)*4==26112
 return np.array(ans,dtype='i8')

W=json.loads((here/'witness.json').read_text())
stock=W['fibres']
E={tuple(v)for v in R['EE']}
assert len(E)==240
def dot(x,y):
 return sum(a*b for a,b in zip(x,y))
assert all(dot(x,x)==8 for x in E)
assert all(dot(x,y)<=4 for x,y in combinations(E,2))
p=[int(r['phase'])for r in stock];T=[tuple(r['tail2'])for r in stock]
blk=[r.get('clip',[])for r in stock]
assert all(S.get(u)==6 for u in p) and len(set(p))==len(p)
assert all(v in E and next(x for x in v if x)>0 for v in T) and len(set(T))==len(T)
assert all(len(L)==len(set(L)) and all(i!=j and 0<=i<len(T) for j in L) for i,L in enumerate(blk))
for i,j in combinations(range(len(stock)),2):
 d=abs(dot(T[i],T[j]))
 q=S[p[i]^p[j]]
 assert d in (0,4) and q in (4,6)
 if d==4 and q==4:
  assert i in blk[j] or j in blk[i]
 else:
  assert i not in blk[j] and j not in blk[i]
n=0
for i,r in enumerate(stock):
 if not blk[i]:
  count=136704
 else:
  j=r['chart']
  q=raw[j][p[i]]
  assert q[0] and q[1] and q[0]!=q[1] and q[2]==0
  V=make_triple(q)
  assert len(V)*8==110592
  VV=two_occupied(q)
  mask=np.ones(len(V),dtype='bool')
  mask2=np.ones(len(VV),dtype='bool')
  for k in blk[i]:
   d=raw[j][p[i]^p[k]]
   # Raw norm8 is the only possible short distance in this B=0 dependent coset.
   assert d[2]==0 and (d[0]or d[1]) and (not d[0]or not d[1]or d[0]==d[1])
   mask=cutgood(V,d,mask)
   mask2=cutgood(VV,d,mask2)
  count=int(mask.sum())*8+int(mask2.sum())*4
  assert count==r['count_per_sign']
 n+=count
bound=6218175600+240+2*n
assert bound==W['lower_bound'] and bound-6228106200==W['gain']
print('PASS 80D',bound,'gain',W['gain'],'full',sum(not x for x in blk),'partial',sum(bool(x)for x in blk))

