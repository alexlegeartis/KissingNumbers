import os
import numpy as np,time
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
p='research/mixed_tier17_eval_subgroup_20261105/cache/'
r='research/mixed_tier17_eval_subgroup_20261105/seed/'
N=3968
V=np.fromfile(p+'mixed17_eval4_orbits.bin',np.int8).reshape(-1,16).astype(np.int16)
H=np.fromfile(p+'mixed17_eval4_group.bin',np.uint16)
F=np.fromfile(p+'mixed17_eval4_flat_reps.bin',np.uint16)
H0=np.fromfile(r+'mixed17_rank4_group.bin',np.uint16)
S=lambda H:np.array([[1-2*(int(h)>>i&1) for i in range(16)] for h in H],np.int16)
S0,S1=S(H0),S(H)
T=np.fromfile(r+'mixed17_seed_tier.bin',np.int8).reshape(-1,16).astype(np.int16)
neg_t={tuple(y) for t in T for y in t*S0}
flat0=np.fromfile(r+'mixed17_seed_flatminus.bin',np.uint16)
neg_f={int(w)^int(h) for w in flat0 for h in H0}
neg=[]
for i,w in enumerate(F):
 orb={int(w)^int(h) for h in H}
 if orb<=neg_f:neg.append(i);neg_f-=orb
for i,v in enumerate(V):
 orb={tuple(y) for y in v*S1}
 if orb<=neg_t:neg.append(128+i);neg_t-=orb
assert not neg_t and not neg_f and len(neg)==44
pos=np.fromfile(p+'baseline_seed.bin',np.int32)
base=np.array(list(pos)+[N+x for x in neg],np.int32)
E=np.fromfile(p+'graph_edges.bin',np.int32).reshape(-1,2)
A=np.zeros((2*N,2*N),bool)
A[E[:,0],E[:,1]]=1;A[E[:,0]+N,E[:,1]+N]=1
opp=np.fromfile(p+'mixed17_eval4_opposite.bin',np.uint8).reshape(3840,128)
ti,fl=np.where(opp)
A[128+ti,N+fl]=1; A[fl,N+128+ti]=1
A[np.arange(128),N+np.arange(128)]=1
A|=A.T
assert not A[np.ix_(base,base)].any(), 'baseline caps conflict'
np.savetxt(p+'negative_seed.txt',neg,fmt='%d')
# Pair each subgroup orbit with its other half under full sign negation.
fm={int(w)^int(h):i for i,w in enumerate(F) for h in H}
tm={tuple(v*s):128+i for i,v in enumerate(V) for s in S1}
partner=np.array([fm[int(w)^65535] for w in F]+[tm[tuple(-v)] for v in V],np.int32)
assert np.all(partner[partner]==np.arange(N))
def mate(v):return int(v//N*N+partner[v%N])
# Search by neutral and improving simultaneous replacements on both caps.
rng=np.random.default_rng(int(os.environ.get('EVAL17_SEED','171448'))); start=time.time(); states=set(); moves=0; best=len(base)
for epoch in range(int(os.environ.get('EVAL17_EPOCHS','1000'))):
 m=[]
 for row in A[:,base]:
  x=0
  for j in np.flatnonzero(row):x|=1<<int(j)
  m.append(x)
 isbase=np.zeros(2*N,bool);isbase[base]=1
 broken=[int(v) for v in base if not isbase[mate(int(v))]]
 cand=np.flatnonzero(~isbase)
 cost=np.array([v.bit_count() for v in m]);candidates=cand[cost[cand]<=14]
 for trial in range(int(os.environ.get('EVAL17_TRIALS','6000'))):
  seed=int(rng.choice(candidates));vict=m[seed]
  mateid=mate(seed)
  for j,b in enumerate(base):
   if b==mateid:vict|=1<<j
  for _ in range(int(rng.integers(1,8))):
   sample=rng.choice(candidates,size=min(200,len(candidates)),replace=False)
   choices=[v for v in sample if 0<(m[int(v)]|vict).bit_count()-vict.bit_count()<=3 and (m[int(v)]|vict).bit_count()<=28]
   if not choices:break
   vict|=m[int(rng.choice(choices))]
  # Preserve existing broken halves and forbid reuniting them; this
  # prevents the walk from merely moving a single missing half around.
  if len(base)<88 and broken:
   locked=broken[:min(2,len(broken))]
   if any(vict>>int(np.flatnonzero(base==v)[0])&1 for v in locked):continue
  else:locked=[]
  r=vict.bit_count();ids=np.array([int(v) for v in cand if m[int(v)] & ~vict==0])
  if len(ids)<r-1 or len(ids)>200:continue
  D=A[np.ix_(ids,ids)];i,j=np.where(np.triu(D,1))
  M=coo_matrix((np.ones(2*len(i)),(np.repeat(np.arange(len(i)),2),np.c_[i,j].ravel())),shape=(len(i),len(ids))).tocsr()
  lo=np.zeros(len(ids));hi=np.ones(len(ids));lo[np.flatnonzero(ids==seed)]=1
  hi[np.flatnonzero(ids==mateid)]=0
  for v in locked:hi[np.flatnonzero(ids==mate(v))]=0
  if np.any(lo>hi):continue
  res=milp(-np.ones(len(ids)),integrality=np.ones(len(ids)),bounds=Bounds(lo,hi),constraints=LinearConstraint(M,-np.inf,1),options={'time_limit':.5})
  if res.fun is None:continue
  chosen=ids[res.x>.5];gain=len(chosen)-r
  if gain>=0 or (gain == -1 and len(base) == 88 and rng.random()<0.10):
   kept=np.array([int(b) for j,b in enumerate(base) if not(vict>>j&1)])
   new=np.sort(np.r_[kept,chosen]).astype(np.int32); key=tuple(new)
   if key in states:continue
   assert not A[np.ix_(new,new)].any()
   states.add(key);base=new;moves+=1
   print('move',moves,'gain',gain,'cap orbit count',len(base),'removed',r,'trial',trial,'time',round(time.time()-start,1),flush=True)
   np.save(p+'two_cap_walk.npy',base)
   if len(base)>best:
    np.save(p+'two_cap_POTENTIAL.npy',base);print('POTENTIAL, requires full exact verification',flush=True);raise SystemExit
   break
 else:
  print('no neutral move epoch',epoch,'time',round(time.time()-start,1),flush=True);break
