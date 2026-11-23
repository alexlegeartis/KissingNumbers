#!/usr/bin/env python3
"""Exact coupled horizontal profile replay on previously published 77–80D tails."""
import json,sys,runpy, subprocess,struct, tempfile
from pathlib import Path
from hashlib import sha256
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here))
from two_profiles import improve as recover_profiles
parent=here.parent/'dim77-80-gamma72-apn-coupled-tails'
pins={'claims.json': '9fa9b91c6017cc725f8b3feb5e4c02ebd3b89f29683e9926a656ddab36ae02cc', 'q4tail.py': '7604b007b6a8a660afff10034292515091a1d61f77ef201e0a4d1ec1da04c325', 'witness.json': '33ee860295e2deb8ecb45f162da5280c438393ea6f051b7c9ed24c07e6c12884', 'verify.py': '72126ebd1834a88bd65119d5ed74ca159fcf3962c1a262068522c6f3280825dc'}
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

for d, w in sorted(R['W'].items()):
 m=int(d)-72;V=R['tail'](m)
 rows=[];clips=[];present=[]
 for j,c,t in w['rays']:
  v=V[j]; a=H[t]
  bl=[a^u for x,b,u in rows if b==1+c and np.dot(v,x)!=0]
  clips.append(bl);present.append(rawmask(a,bl));rows.append((v,1+c,a))
 reclaimed=recover(rows,present,list(reversed(range(len(rows)))))
 assert 16*reclaimed==w['reclaimed']
 k=optimize(rows,present)
 a,b,mat,stats=recover_profiles(R,rows,clips,d)
 out[d]=R['res'][d]+16*k+8*(b-a)
 print('PASS coupled profiles',d,out[d],'new triples',k,'profile counts',a,b,'new profile contacts',8*(b-a),flush=True)
cf=here/'claims.json'
if cf.exists():
 claim=json.loads(cf.read_text())["lower_bounds"];assert {d:claim[d] for d in out}==out
else:print(json.dumps({'lower_bounds':out},indent=2))
