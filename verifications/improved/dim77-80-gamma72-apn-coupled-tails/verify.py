"""Exact all-root APN replay of 77D–80D."""
from pathlib import Path
from itertools import combinations
from hashlib import sha256
import sys, runpy,json
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
prev=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here.parent/'dim81-85-gamma72-apn-multipartite'
pins={'verify.py':'725066df79cf01b22204118096994714711d82869b83094e7a7e843543eebbbd','witness.json':'e2978ae01bc358549003d331edc65a0be1235fa01d49ac43500e09c56bd3a4a7'}
assert all(sha256((prev/f).read_bytes()).hexdigest()==h for f,h in pins.items())
base=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else prev.parent/'dim73-82-gamma72-depth6-lattice'
sys.argv=[str(prev/'verify.py'),str(base)]
R=runpy.run_path(sys.argv[0])
globals().update({k:R[k]for k in ('P','H','K','M','Or','Blist','RA','RB','T','MD0')})
a,b=K[:2]
HH={t:i for i,t in enumerate(H)}
W=json.loads((here/'witness.json').read_text())


from functools import lru_cache
F=P['F']; A0=P['A'];Bmat=P['B'];WW=P['W'];pk=P['p24']
R4={}
for aa,bb,v in zip((WW@A0%2)@pk,(WW@Bmat%2)@pk,WW): R4[int(aa),int(bb)]=v
assert len(R4)==98280
@lru_cache(None)
def getpair(a,b):
 for u,S0 in M.items():
  v=u^a
  if not v:continue
  for q in S0:
   if (v,q^b) in R4:return R4[u,q],R4[v,q^b]
 raise AssertionError((a,b))
@lru_cache(None)
def norm6(a,b):
 x,y=getpair(a,b);d=int(x@F@y);assert abs(d)==1
 z=x-d*y
 assert int(z@F@z)==6 and int((z@A0%2)@pk)==a and int((z@Bmat%2)@pk)==b
 return z
QS8=set();QS6=set(); layouts={}
for t in H:
 U48=[];U66=[]
 for h in (0,a,b):
  x,y=[z for z in (h^a,h^b,h) if z]
  for u,v in ((x,y),(y,x)):
   for j in sorted(M[u]):
    l=j^t
    assert l not in M[v] and Or(v,l)
    U48.append((u,j,v,l));QS8.add((v,l))
  for j in Blist:
   l=j^t
   if not Or(x,j) and not Or(y,l):
    U66.append((x,j,y,l));QS6.update(((x,j),(y,l)))
 assert 24*len(U48)+len(U66)==6528,(len(U48),len(U66))
 layouts[t]=U48,U66
A8=sorted(QS8|{(0,d)for d in H[1:]});I8={p:i for i,p in enumerate(A8)}
A6=sorted(QS6);I6={p:i for i,p in enumerate(A6)}
assert 24*int(np.max(np.abs(WW)))*int(np.max(np.abs(WW@F)))*2<2**53
V8=np.empty((len(A8),24,24),dtype='i8')
for k in range(0,len(A8),64):
 C=A8[k:k+64]; pairs=[getpair(*s)for s in C]
 xx,yy=np.array(pairs,dtype='i8').transpose(1,0,2)
 assert np.all(((xx@F)*yy).sum(axis=1)==0)
 nn=xx+yy
 D=(nn@F).astype('f8')@WW.T.astype('f8')
 for j,(aa,bb) in enumerate(C):
  z=nn[j]-2*WW[D[j]==4];z=np.vstack([z,nn[j],-nn[j]])
  assert len(z)==len(np.unique(z,axis=0))==48
  assert np.all(((z@F)*z).sum(axis=1)==8)
  assert np.all((z@A0%2)@pk==aa) and np.all((z@Bmat%2)@pk==bb)
  i=np.argmax(z!=0,axis=1)
  good=z[z[np.arange(48),i]>0];assert len(good)==24
  V8[k+j]=good
V6=np.array([norm6(*u) for u in A6],dtype='i8')
V8F=V8.reshape(-1,24).astype('f8');V6F=V6.astype('f8')
assert 24*max(int(np.max(np.abs(V8))),int(np.max(np.abs(V6))))*int(np.max(np.abs(V8.reshape(-1,24)@F)))<2**53
BAD8=np.zeros((128,len(A8)),dtype='u4');BAD6=np.zeros((128,len(A6)),dtype='bool')
def pack24(e):
 q=np.packbits(e,axis=1,bitorder='little').astype('u4')
 return q[:,0]|(q[:,1]<<8)|(q[:,2]<<16)
for i,d in enumerate(H[1:],1):
 g=(V8[I8[(0,d)]]@F).T.astype('f8')
 X=np.abs(V8F@g)==4
 BAD8[i]=pack24(X.any(axis=-1).reshape(-1,24))
 BAD6[i]=(np.abs(V6F@g)==4).any(axis=-1)
for t,(U,V) in layouts.items():
 layouts[t]=(np.array([(RA[x],RB[j],I8[(y,l)])for x,j,y,l in U]),np.array([(I6[(x,j)],I6[(y,l)])for x,j,y,l in V]))
print('READY full shell classes',len(A8),len(A6),flush=True)
def shell(t,bs):
 U,V=layouts[t]; r1,r2,r3=U.T
 bad4=np.zeros(len(U),dtype='bool');bad8=np.zeros(len(U),dtype='u4');bad6=np.zeros(len(V),dtype='bool')
 for delta in bs:
  ii=HH[delta]
  bad4|=MD0[r1,r2^RB[delta]]
  bad8|=BAD8[ii,r3]
  bad6|=BAD6[ii,V].any(axis=1)
 good48=sum(24-int(v).bit_count()for v in bad8[~bad4])
 good66=len(V)-int(bad6.sum())
 return 4*good48,4*good66

from itertools import product
def tail(m):
 D=[]
 if m==5:
  for i,j in combinations(range(5),2):
   for c in (-2,2):
    x=[0]*5;x[i]=2;x[j]=c;D.append(tuple(x))
  return sorted(D)
 for i,j in combinations(range(8),2):
  for c in (-2,2):
   x=[0]*8;x[i]=2;x[j]=c;D.append(tuple(x))
 for s in product((-1,1),repeat=7):
  z=(1,*s)
  if s.count(-1)%2==0:D.append(z)
 if m==6:D=[x for x in D if x[5]==x[6]==x[7]]
 if m==7:D=[x for x in D if x[6]==x[7]]
 return sorted(D)
cmp={77:6223198744,78:6225814200,79:6226634478,80:6228431712}
def tailmat(m):
 s=tail(m);M=np.array(s,dtype='i8');G=M@M.T
 assert len(s)=={5:20,6:36,7:63,8:120}[m]
 assert len(set(s))==len(s) and all(next(v for v in x if v)>0 for x in s)
 assert np.all(np.diag(G)==8)
 G=G-np.diag(np.diag(G));assert np.isin(G,[-4,0,4]).all()
 if m==6:assert all(x[5]==x[6]==x[7] for x in s)
 if m==7:assert all(x[6]==x[7] for x in s)
 return s,G!=0

def corekeys(A,B):
 return (A[:,2].astype('i8')<<36)|(B[:,0].astype('i8')<<24)|(B[:,1].astype('i8')<<12)|B[:,2].astype('i8')
KEY={}
for t in H:
 A,B=T[RB[t]]
 K0=corekeys(A,B); idx=np.argsort(K0);KEY[t]=(K0,K0[idx],idx)
@lru_cache(None)
def tripleedges(t,u):
 assert t!=u
 A,B=T[RB[t]];K0,_,_=KEY[t];_,target,idx=KEY[u]
 d=RB[t^u];ans=[]
 for b in range(3):
  test=K0^(d<<(24-12*b))
  p=np.searchsorted(target,test)
  p=np.minimum(p,len(target)-1)
  k=np.flatnonzero(target[p]==test)
  ans.append((k,idx[p[k]]))
 x,y=np.concatenate([p for p,q in ans]),np.concatenate([q for p,q in ans])
 return x,y
def reclaim(Ds,AD,expected):
 n=len(Ds); NN=13824;keep=[np.ones(NN,dtype='bool')for _ in range(n)]
 links=[]
 for j,(u,cu,tu)in enumerate(Ds):
  for i,(v,cv,tv)in enumerate(Ds[:j]):
   if cu==cv and AD[u,v]:
    x,y=tripleedges(H[tv],H[tu]);links.append((i,j,x,y));keep[j][y]=False
 orig=sum(int(x.sum())for x in keep)
 assert orig==expected
 for j in reversed(range(n)):
  forbid=np.zeros(NN,dtype='bool')
  for ii,jj,x,y in links:
   if j==jj:forbid[y[keep[ii][x]]]=True
   if j==ii:forbid[x[keep[jj][y]]]=True
  keep[j]|=~forbid
 for i,j,x,y in links:assert not np.any(keep[i][x]&keep[j][y])
 return (sum(int(x.sum())for x in keep)-orig)*16

out={}
for dk in sorted(W):
 d=int(dk); m=d-72; cur=W[dk]
 assert m in (5,6,7,8)
 rays,AD=tailmat(m);n=len(rays);Ds=cur['rays']
 assert len(Ds)==n and {v[0]for v in Ds}==set(range(n))
 assert len({(c,t)for j,c,t in Ds})==n and all(c in (0,1,2) and 0<=t<128 for j,c,t in Ds)
 s=0;partial=0;triple=0
 for i,(j,c,t) in enumerate(Ds):
  bl=[H[t^v]for k,col,v in Ds[:i] if col==c and AD[j,k]]
  if not bl:s+=136704;triple+=13824;continue
  A,B=T[RB[H[t]]];bad=np.zeros(len(A),dtype='bool')
  for dd in bl:
   assert dd
   bad|=MD0[A,B^RB[dd]].any(axis=1)
  tr=len(A)-int(bad.sum());triple+=tr
  s+=8*tr+sum(shell(H[t],bl));partial+=1
 bound=6218175600+2*s+2*n
 assert d in cmp and bound==cur['lower_bound'] and bound-cmp[d]==cur['gain']>0
 extra=reclaim(Ds,AD,triple)
 assert extra==cur['reclaimed']
 assert bound+extra==cur['final'] and cur['gain']+extra==cur['final_gain']
 out[dk]={'lower_bound':bound+extra,'prior':cmp[d],'gain':cur['final_gain']}
 print('PASS coupled tail roots',d,bound+extra,'gain',cur['final_gain'],'rays',n,'reclaimed',extra,flush=True)
assert set(map(str,range(77,81)))==set(out)

# Exact norm-four additional tails in distinct B-only shifts. The
# old roots of D5 are padded with three zeros for this ambient check.
groups=runpy.run_path(str(here/'q4tail.py'))['groups']
assert H[0]==0 and len(H)==len(set(H))==128
res={}
for m,W4 in groups():
 d=str(m+72);assert len(W4)=={5:26,6:30,7:84,8:240}[m]
 root=np.asarray(tail(m),dtype='i8')
 if m==5:root=np.pad(root,((0,0),(0,3)))
 V=np.asarray([z for z,h in W4],dtype='i8')
 assert np.max(np.abs(V@root.T))<=8
 assert np.all(np.sum(V*V,axis=1)==16)
 for a in {h for z,h in W4}:
  assert 0<a<len(H)
  raw=[(z>>i)&1 for z in (H[a],0,0)for i in range(24)]
  assert P['syndrome'](raw)==(0,0,H[a])
  P['prove'](raw,4)
  assert P['contact'](raw,4)==144
 res[d]=out[d]['lower_bound']+144*len(W4)
 print('PASS with norm4 tails',d,res[d],len(W4),flush=True)
assert {'status':'PASS','lower_bounds':res}==json.loads((here/'claims.json').read_text())
