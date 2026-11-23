#!/usr/bin/env python3
"""Exact Barnes–Wall-type equality sections as Γ72 tail stocks."""
import json,sys,runpy,subprocess,struct,tempfile
from itertools import combinations,product
from pathlib import Path
from hashlib import sha256
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
helper=here.parent/'dim81-85-gamma72-half-sign'
helper_pins={'two_profiles.py': '0159dc533e090920cc3371b4744437e9371d5868ecbe3dcab103e29a88ca70d4', 'component.cpp': '3559c988348f3064a8910709563835b02bf985115de1099470110ff2eed68e15'}
assert all(sha256((helper/n).read_bytes()).hexdigest()==s for n,s in helper_pins.items())
sys.path.insert(0,str(helper))
from two_profiles import improve as recover_profiles
parent=here.parent/'dim81-85-gamma72-apn-complete-dm'
pins={'verify.py':'f9888cb194f2b16624bbabbccdedb802800b13c266a06a12710aef66fc190e33','witness.json':'5f89656cd4916c2534f9ed4c372ae625cf95d6fe82bfca50dd91af1b911ee334'}
assert all(sha256((parent/n).read_bytes()).hexdigest()==s for n,s in pins.items())
sys.argv=[str(parent/'verify.py')]
R=runpy.run_path(sys.argv[0])
H=R['H'];rb=R['RB'];ts=R['T'];md0=R['MD0']
out={}
hh=np.array([rb[t]for t in H],dtype='u2')
quo=np.min(np.arange(4096,dtype='u2')[:,None]^hh,axis=1)
vc,qids=np.unique(quo,return_inverse=True);assert len(vc)==32
dtype=np.dtype([('g','<u4'),('r','<u2'),('b0','<u2'),('b1','<u2'),('b2','<u2'),('i','<u4')])
assert dtype.itemsize==16
temp=tempfile.TemporaryDirectory()
exe=Path(temp.name)/'component'
subprocess.run(['g++','-O3','-std=c++17',str(helper/'component.cpp'),'-o',str(exe)],check=True)
def optimize(rows,active):
 n=len(rows);k=13824
 V=np.array([u for u,_,_ in rows]);col=np.array([a for _,a,_ in rows])
 edges=(col[:,None]==col[None,:])&(V@V.T!=0)
 np.fill_diagonal(edges,False)
 W=np.empty((n,k),dtype)
 for i,(_,c,t) in enumerate(rows):
  A,B=ts[rb[t]]
  W[i]['g']=((c*4096+A[:,2].astype('u4'))<<10)|(qids[B[:,0]]<<5)|qids[B[:,1]]
  W[i]['r']=i
  for j in range(3):W[i][f'b{j}']=B[:,j]
  W[i]['i']=i*k+np.arange(k)
 seed=np.asarray(active,dtype='u1')
 inp=struct.pack('<2I',n,n*k)+edges.astype('u1').tobytes()+W.tobytes()+seed.tobytes()
 q=subprocess.run([str(exe)],input=inp,check=True,capture_output=True)
 print(q.stderr.decode().strip(),flush=True)
 optimum=np.frombuffer(q.stdout,dtype='u1').reshape((n,k))
 assert np.max(optimum)<=1
 more=int(optimum.sum())-int(seed.sum())
 assert more>=0
 if __import__('os').environ.get('EXPORT_TRIPLE'):
  np.save(Path(__import__('os').environ['EXPORT_TRIPLE'])/f"triple{n}.npy",optimum)
 return more
# A triple row is determined by A-label 2 and B-labels 0 and 1.
def keys(t):
 A,B=ts[rb[t]]
 x=(A[:,2].astype('u8')<<24)|(B[:,0].astype('u8')<<12)|B[:,1].astype('u8')
 assert len(np.unique(x))==len(x)
 ix=np.argsort(x);return x,ix,x[ix]
KK={t:keys(t)for t in H}
def rawmask(t,bl):
 A,B=ts[rb[t]];bad=np.zeros(len(A),bool)
 for z in set(bl):bad|=md0[A,B^rb[z]].any(axis=1)
 return ~bad
def recover(rows,present,order):
 Vs=np.array([u for u,_,_ in rows]); labels=np.array([a for _,a,_ in rows])
 dots=Vs@Vs.T
 edges=(labels[:,None]==labels[None,:])&(dots!=0)
 np.fill_diagonal(edges,False)
 ns=[np.flatnonzero(row) for row in edges]
 add=0
 for i in order:
  _,_,t=rows[i];X=KK[t][0]
  eligible=~present[i]
  for j in ns[i]:
   _,_,u=rows[j];k=rb[t^u];assert k
   I=np.flatnonzero(eligible)
   if not len(I):break
   _,perm,S0=KK[u];Q=X[I][None,:]^np.array([0,k,k<<12],dtype='u8')[:,None]
   J=np.minimum(np.searchsorted(S0,Q),len(S0)-1)
   eligible[I]&=~((S0[J]==Q)&present[j][perm[J]]).any(axis=0)
  present[i]|=eligible
  add+=int(eligible.sum())
 return add
PART={10:[5250,256,8,16384,10304,512,4,1,16,32768,32],
      11:[10240,32768,128,512,20,40,5120,2,64,1,16384,256],
      13:[0,256,16,130,1024,4096,2048,65,8,33280,16384,8192,32,4]}
def tails(m):
 p=PART[m];z,g=p[0],p[1:]
 assert len(g)==m and sum(g)+z==65535 and (z|int(np.bitwise_or.reduce(g)))==65535
 S=set()
 for i,j in combinations(range(16),2):
  for e in [-2,2]:
   x=[0]*16;x[i]=2;x[j]=e;S.add(tuple(x))
 for a in range(1,16):
  for b in (0,1):
   p=[j for j in range(16)if (a&j).bit_count()%2==b]
   for s in product([-1,1],repeat=7):
    if s.count(-1)%2==0:
     x=[0]*16
     for i,v in zip(p,(1,*s)):x[i]=v
     S.add(tuple(x))
 def ok(v):
  if any(v[i]for i in range(16)if (z>>i)&1):return False
  return all(len({v[i]for i in range(16)if (a>>i)&1})==1 for a in g)
 V=np.array(sorted(v for v in S if ok(v)),dtype='i8')
 assert np.all(np.diag(V@V.T)==8) and np.max(np.abs(V@V.T-8*np.eye(len(V),dtype='i8')))<=4
 return V

phases=json.loads((here/'phases.json').read_text())
for dm,raymap in sorted(phases.items()):
 m=int(dm)-72;V=tails(m)
 assert len(V)=={10:168,11:219,13:453}[m]
 assert len({j for j,c,t in raymap})==len(raymap)
 assert len({(c,t) for j,c,t in raymap})==len(raymap)
 assert all(0<=j<len(V) and c in range(3) and t in range(128) for j,c,t in raymap)
 rows=[];clips=[];present=[]
 for j,c,t in raymap:
  v=V[j]; a=H[t]
  bl=[a^u for x,b,u in rows if b==1+c and np.dot(v,x)!=0]
  clips.append(bl);present.append(rawmask(a,bl));rows.append((v,1+c,a))
 raw=sum(int(x.sum())for x in present)
 order=sorted(range(len(rows)),key=lambda i:(int(present[i].sum()),i))
 more=recover(rows,present,order)
 k=optimize(rows,present)
 a,b,mat,stats=recover_profiles(R,rows,clips,dm)
 out[dm]=6218175600+16*(raw+more+k)+8*b+2*len(rows)
 print('PASS BW tail section',dm,out[dm],'raw triple',raw,'recoupled',more,'component gain',k,'two-profile counts',a,b,'rays',len(rows),flush=True)
claim=here/'claims.json'
if claim.exists():assert json.loads(claim.read_text())['lower_bounds']==out
else:print(json.dumps({'lower_bounds':out},indent=2))
