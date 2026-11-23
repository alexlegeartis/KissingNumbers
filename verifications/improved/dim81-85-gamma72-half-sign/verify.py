#!/usr/bin/env python3
"""Exact signed half-support lift with all compatible raw norm-twelve profiles."""
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import json,sys,runpy, subprocess,struct, tempfile
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here))
from two_profiles import improve as recover_profiles
parent=here.parent/'dim81-85-gamma72-norm4-layer'
pins={'verify.py':'ee22054c256686598c8b3429d92740478d3cc82a4abc99944ec08a6b2fb58c27',
      'claims.json':'597b65f9a695807b17fc719c3258f0b41a687e548f119081ced535606c141ee6'}
assert all(sha256((parent/n).read_bytes()).hexdigest()==s for n,s in pins.items())
sys.argv=[str(parent/'verify.py')]
old=runpy.run_path(sys.argv[0])
R=old['R']; RR=R['R']
H=R['H']; rb=R['RB']; ts=R['T']; md0=R['MD0']; S=R['shell']
assert len(H)==len(set(H))==128 and H[0]==0 and all(a^b in H for a in H for b in H)
phase=json.loads((here/'phases.json').read_text())
alphabet='0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-_'
out={}
hh=np.array([rb[t]for t in H],dtype='u2')
quo=np.min(np.arange(4096,dtype='u2')[:,None]^hh,axis=1)
vc,qids=np.unique(quo,return_inverse=True);assert len(vc)==32
dtype=np.dtype([('g','<u4'),('r','<u2'),('b0','<u2'),('b1','<u2'),('b2','<u2'),('i','<u4')])
assert dtype.itemsize==16
temp=tempfile.TemporaryDirectory()
exe=Path(temp.name)/'component'
subprocess.run(['g++','-O3','-std=c++17',str(here/'component.cpp'),'-o',str(exe)],check=True)
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
old_counts={'81':688235,'82':643182,'83':599585,'84':1533245,'85':1498461}
for rec in RR['W']['rows']:
 d=str(rec['dimension']); m=int(d)-72
 a=[tuple(range(8))]
 if m>=12: a.extend([(0,1,2,3,8,9,10,11),(4,5,6,7,8,9,10,11)])
 def spin(s,x):
  v=[0]*m;v[s[0]]=1
  for i,j in enumerate(s[1:]):v[j]=1-2*((x>>i)&1)
  return v
 U=np.array([spin(s,x) for s in a for x in range(128) if x.bit_count()%2==0],dtype='i8')
 assert np.all(np.sum(U*U,axis=1)==8)
 assert len(U)==len(set(map(tuple,U))) and np.max(np.abs(U@U.T-8*np.eye(len(U),dtype='i8')))<=4
 Ds=[z[:5] for z in rec['rays']]+[z[:5] for z in R['W'][d]['new']]
 rows=[]; present=[]; clips=[]
 for r,(i,j,e,c,t) in enumerate(Ds):
  v=[0]*m;v[i]=2;v[j]=2 if e else -2
  if r<len(rec['rays']):
   bl=[t^Ds[k][4]for k in rec['rays'][r][6]]
  else:
   bl=[t^s for u,b,s in rows if c==b and np.dot(v,u)!=0]
  present.append(rawmask(t,bl)); clips.append(bl)
  rows.append((v,c,t))
 assert np.all(np.abs(np.array([q[0]for q in rows])@U.T)<=4)
 want=phase[d];assert len(want)==3*len(U)
 ID=set(); tri=0; two=0; total_profiles=[0,0]
 assert sum(int(v.sum()) for v in present) <= len(present)*13824
 shifted={(c,t) for _,c,t in rows}
 for i in range(0,len(want),3):
  x=sum(alphabet.index(ch)<<(6*j)for j,ch in enumerate(want[i:i+3]))
  k,x=divmod(x,384);c,t=divmod(x,128);c+=1;t=H[t]
  assert 0<=k<len(U) and k not in ID and (c,t)not in shifted
  ID.add(k);shifted.add((c,t))
  v=U[k]; bl=[t^s for u,b,s in rows if c==b and np.dot(v,u)]
  assert 0 not in bl
  keep=rawmask(t,bl);present.append(keep);tri+=int(keep.sum());clips.append(bl)
  vals=S(t,bl);two+=sum(vals)
  for j,n in enumerate(vals):total_profiles[j]+=n
  rows.append((v,c,t))
 assert len(ID)==len(U)
 # Old pure tails (norm 8) have dot ≤4 by Cauchy–Schwarz.
 # q4 old tails: ||w||_1≤4 and |v_i/2|≤1/2, so dot≤2.
 # Their B-only phase differs from each nonzero colour by depth ≥6.
 # Same-ray antipodes and old Γ72 are checked by norm/depth (8 and 6+2).
 assert tri==old_counts[d]
 order=list(range(len(rows)))
 # Most heavily clipped rays get first choice; tie-break without any search.
 order.sort(key=lambda i:(int(present[i].sum()),i))
 more=recover(rows,present,order)
 optimized=optimize(rows,present)
 profile0,profile1,nmask,stats=recover_profiles(R,rows,clips,d)
 gain=16*(tri+more+optimized)+2*two+8*(profile1-profile0)
 print("two-profile replay",d,profile0,profile1,nmask,stats,flush=True)
 out[d]=old['res'][d]+gain
 print('PASS half-sign',d,out[d],'gain',gain,'triples',tri,'recoupled',more,'component gain',optimized,'two-block',2*two,
       'profile counts',total_profiles,flush=True)
cf=here/'claims.json'
if cf.exists():
 claim=json.loads(cf.read_text())["lower_bounds"];assert {d:claim[d] for d in out}==out
else:print(json.dumps({'lower_bounds':out},indent=2))
