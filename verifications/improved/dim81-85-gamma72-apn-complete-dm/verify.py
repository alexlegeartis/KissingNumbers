"""Exact all-root APN replay of 81D–85D."""
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

def conflict(x,y):
 return x[3]==y[3] and len(set(x[:2])&set(y[:2]))==1
out={}
for rec in R['W']['rows']:
 d=rec['dimension']; m=d-72; cur=W[str(d)]
 old=rec['rays']; Ds=[[i,j,e,c,t,n,z[:]]for i,j,e,c,t,n,z in old]
 for v in cur['new']:
  i,j,e,c,t=v
  bl=[k for k,q in enumerate(Ds)if conflict(v,q)]
  Ds.append([*v,None,bl])
 assert {(i,j,e)for i,j,e,*_ in Ds}=={(i,j,e)for i,j in combinations(range(m),2)for e in (0,1)}
 assert len(Ds)==m*(m-1)
 assert len({(x[3],x[4])for x in Ds})==len(Ds)
 assert all(x[3]in (1,2,3) and x[4]in H for x in Ds)
 assert all(conflict(x,y)==(j in x[6] or i in y[6]) for i,j in combinations(range(len(Ds)),2)for x,y in [(Ds[i],Ds[j])])
 full=0;extra=0
 for i,v in enumerate(Ds):
  _,_,_,_,t,oldn,bl=v
  if not bl:n=136704
  else:
   A,B=T[RB[t]];bad=np.zeros(len(A),dtype='bool')
   delta=[t^Ds[k][4]for k in bl]
   for z in delta:bad|=MD0[A,B^RB[z]].any(axis=1)
   n=8*(len(A)-int(bad.sum()))
   extra+=sum(shell(t,delta))
  assert i>=len(old) or oldn==n
  full+=n
 bound=6218175600+2*(full+extra)+2*m*(m-1)
 assert 2*extra==cur['two_block_gain']
 assert bound==cur['lower_bound'] and bound-rec['lower_bound']==cur['gain']>0
 out[str(d)]={'lower_bound':bound,'prior':rec['lower_bound'],'gain':cur['gain']}
 print('PASS all roots',d,bound,'gain',cur['gain'],'two-block gain',2*extra,flush=True)
assert set(out)==set(W)==set(map(str,range(81,86)))
