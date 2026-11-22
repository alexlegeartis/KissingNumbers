"""Exact replay against the pinned harmonic latitude package."""
import hashlib,json,sys,os
from fractions import Fraction as F
from pathlib import Path
from math import ceil
from collections import defaultdict
import numpy as np
HERE=Path(__file__).resolve().parent
PARENT=Path(os.getenv('HARMONIC_PARENT',str(HERE.parent/'dim74-85-gamma72-harmonic-latitudes'))).resolve()
PINS={'harmonic.py':'ce8c00a60ba9d0e9aef651aa8ecfec3769120f4f27fef735c27ffdf21d45b432','model.py':'81cdf9df5b5b121e9437bc3df6ee1d021ba1682d7fea817b9190f4c553e83a61','results.json':'e06e365655d8f888dda38892d158231bc7acab6e65b61888356cafc0159faddd'}
PINS.update({'ico.py':'d2152e7f8fa77a54215f56eeb25398810b64d2ff9dd2a90a75930913b144a693','rays.txt':'4dc4147c9696dbf284017fcac18ef9e5201824fbeed9f391a1c89128f38fd3a0','polygon13.txt':'832ad51fcd367e160611d01352448735ed7534915ccd0501a5f3acfb6d0c10d4'})
for f,s in PINS.items():
 assert hashlib.sha256((PARENT/f).read_bytes()).hexdigest()==s
sys.path.insert(0,str(PARENT))
import model
from harmonic import system,pol

J=(4,5,6);N={j:model.P[j]//2 for j in J}
def std(a,b):return (a,b) if (a,b)>=(0,0) else (-a,-b)
def key(j,s,a,b):return (j,s,*std(a,b))

def rooted_system(witness):
 W=witness['affine'];K=[tuple(map(int,k.split(',')))for k in W['K']];pos={k:i for i,k in enumerate(K)};assert len(K)==len(pos)==328
 specs={};all_rows=[]
 for j in J:
  P,E,k,Z=system(j); assert set(K0 for K0 in K if K0[0]==j)=={key(j,*x)for x in k}
  locations=[pos[key(j,*x)]for x in k];specs[j]=(P,E)
  for z in Z:
   v={i:c for i,c in zip(locations,z[:-1])if c}
   all_rows.append((v,z[-1]))
  # M_2=0: norm-eight harmonic degree fourteen vanishes too.
  assert not [(a,b)for a in range(2)for b in range(2)if 4*a+6*b==2]
  for a in range(15):
   H=pol(a,14-a,j,8)
   def h(x,y):return sum(c*x**r*y**s for (r,s),c in H.items())
   v=defaultdict(F)
   for x,y in P[8]:v[pos[key(j,8,x,y)]]+=h(x,y)
   c=-sum(h(x,y)for x,y in E[8])
   from math import lcm
   den=lcm(F(c).denominator,*(z.denominator for z in v.values()))
   v={i:int(x*den)for i,x in v.items()if x};c=int(c*den)
   if v or c:all_rows.append((v,c))
 for j in J:
  for a in J:
   if j>=a:continue
   for b in range(-10,11):
    v={};c=0
    for p,q,z in [(j,a,N[j]),(a,j,-N[a])]:
     k=key(p,10,q,b)
     if k in pos:v[pos[k]]=z
     elif (q,b) in specs[p][1][10]:c-=z
    if v or c:all_rows.append((v,c))
 den=int(W['den']);verts=[tuple(map(int,v))for v in W['v']]
 assert den>0 and len(verts)==len(K) and all(len(v)==4 for v in verts)
 free=[pos[tuple(map(int,k.split(',')))]for k in W['free']]
 assert len(free)==len(set(free))==3 and all(verts[k]==tuple(den*int(i==j)for i in range(4))for j,k in enumerate(free))
 mod=49999
 mat=np.zeros((len(all_rows),len(K)),dtype='i8')
 for ri,(v,c) in enumerate(all_rows):
  assert all(sum(w*verts[i][j]for i,w in v.items())==(c*den if j==3 else 0)for j in range(4))
  for i,w in v.items():mat[ri,i]=w%mod
 assert rank_mod(mat,mod)==325
 return K,pos,den,verts,specs

def rank_mod(A,p):
 r=0
 for i in range(A.shape[1]):
  nz=np.flatnonzero(A[r:,i])
  if not len(nz):continue
  t=r+int(nz[0]);A[[t,r]]=A[[r,t]]
  A[r,i:]=A[r,i:]*pow(int(A[r,i]),-1,p)%p
  A[r+1:,i:]=(A[r+1:,i:]-A[r+1:,i:i+1]*A[r:r+1,i:])%p
  r+=1
 return r

def objective(d,K,pos,specs):
 k=d-72;shape,X=model.PALETTES[d];ts=model.tails_of(k,shape)
 assert len(ts)==len(set(map(tuple,ts)))
 assert set(map(tuple,ts))==set(map(tuple,-ts))
 root=np.array((*shape,*([0]*(k-len(shape)))),dtype='i8')
 m=int(root@root); A=X.numerator*m; B=8*X.denominator-10*X.numerator; C=4*X.denominator*m
 Q=8-10*X
 assert Q>0 and 25*X<=16 and 3*X+Q<=4
 assert (np.sum(ts*ts,axis=1)==m).all()
 assert 4*(shape[0]+shape[1])**2*Q/m<=16
 def gate(j,x):return abs(A*j+B*x)>C
 a=ts@root
 same=int(gate(10,a).sum())-1
 D=same+sum(int(gate(j,a).sum())*N[j]for j in J)
 obj=[0]*len(K); c=2*same*D
 for j in J:
  P,E=specs[j];bins=P[10]+E[10]
  groups={}
  start={x:gate(x,a)for x,y in bins}
  for t in ts[gate(j,a)]:
   g=tuple(tuple(sorted(t[root==n]))if n else tuple(sorted(abs(t[root==0])))for n in sorted(set(root)))
   if g not in groups:groups[g]=[t,0]
   groups[g][1]+=1
  h=dict.fromkeys(bins,0)
  for t,w in groups.values():
   b=ts@t;end={y:gate(y,b)for x,y in bins}
   for x,y in bins:h[x,y]+=w*int((start[x]&end[y]).sum())
  assert all(h[x,y]==h[-x,-y]for x,y in bins)
  for x,y in P[10]:obj[pos[key(j,10,x,y)]]+=N[j]*h[x,y]
  c+=N[j]*sum(h[x,y]for x,y in E[10])
 return len(ts),D,c,obj

def replay():
 model.independently_replay_parent_counts()
 W=json.loads((HERE/'witness.json').read_text())
 state=rooted_system(W)
 K,pos,den,verts,specs=state
 previous={r['dimension']:r['lower_bound']for r in json.loads((PARENT/'results.json').read_text())}
 import oriented
 import gzip
 with gzip.open(HERE/'oriented.json.gz','rt')as f:rows=json.load(f)
 assert [r['d']for r in rows]==[74,75,76]
 assert all(r['free']==W['affine']['free'] for r in rows)
 out=[]
 for row in rows:
  d=row['d'];T,S,L=oriented.replay(PARENT,state,row,key,model)
  out.append({'dimension':d,'previous':previous[d],'lower_bound':L,'gain':L-previous[d],'triangles':T,'selected':S})
  print(d,L,L-previous[d],flush=True)
 for row in W['results']:
  d=row['n'];k=d-72
  M,D,c,O=objective(d,K,pos,specs)
  V=[sum(w*z[i]for w,z in zip(O,verts)if w)+(c*den if i==3 else 0)for i in range(4)]
  nu=list(map(int,row['w']['num']));de=int(row['w']['den'])
  ids=[pos[tuple(map(int,z.split(',')))]for z in row['w']['ids']]
  assert len(ids)==len(nu)==3 and len(set(ids))==3 and min(nu)>=0 and de>0
  target=[de*x+sum(w*verts[t][i]for w,t in zip(nu,ids))for i,x in enumerate(V)]
  assert target[:3]==[0,0,0] and target[3]>0
  T=ceil(F(target[3],4*de*den))
  S=2*ceil(F(model.N//2*M,D+1)*model.low_factor(D,T))
  L=model.N8+S+2*k*(k-1)
  assert (M,D,T,S,L)==tuple(int(row[x])for x in ('M','D','T','S','lower'))
  out.append({'dimension':d,'previous':previous[d],'lower_bound':L,'gain':L-previous[d],'triangles':T,'selected':S})
  print(d,L,L-previous[d],flush=True)
 return out
if __name__=='__main__':
 assert __debug__
 out=replay()
 if '--record' in sys.argv:(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 else:assert json.loads((HERE/'results.json').read_text())==out
