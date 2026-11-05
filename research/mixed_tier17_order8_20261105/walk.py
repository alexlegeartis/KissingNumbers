"""Joint 8-point orbit exchanges in both caps; candidates require exact replay."""
import numpy as np,time,os
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
N=7936;P=2*N
oldpos=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/baseline_seed.bin',np.int32);oldneg=np.loadtxt('research/mixed_tier17_order8_20261105/cache/negative_seed.txt',dtype=np.int32);assert len(oldpos)==len(oldneg)==44
pos=np.array([2*int(v)+j for v in oldpos for j in (0,1)],np.int32);neg=np.array([N+2*int(v)+j for v in oldneg for j in (0,1)],np.int32);base=np.r_[pos,neg];assert len(base)==176
E=np.fromfile('research/mixed_tier17_order8_20261105/cache/graph_edges.bin',np.int32).reshape(-1,2);X=np.fromfile('research/mixed_tier17_order8_20261105/cache/cross_edges.bin',np.int32).reshape(-1,2);A=np.zeros((P,P),bool)
for offset in (0,N):
 A[E[:,0]+offset,E[:,1]+offset]=1;A[E[:,1]+offset,E[:,0]+offset]=1
A[X[:,0],X[:,1]+N]=1;A[X[:,1]+N,X[:,0]]=1
assert not A[np.ix_(base,base)].any(),'baseline conflict in subgroup graph'
# Traverse neutral or one-orbit-loss states instead of restarting from the seed.
rng=np.random.default_rng(170820);seen=set();start=time.time();original=len(base)
for epoch in range(int(os.getenv('EVAL3_EPOCHS','500'))):
 masks=[]
 for row in A[:,base]:
  q=0
  for j in np.flatnonzero(row):q|=1<<int(j)
  masks.append(q)
 isbase=np.zeros(P,bool);isbase[base]=1;cand=np.flatnonzero(~isbase);cost=np.array([q.bit_count() for q in masks]);candidates=cand[cost[cand]<=22];state=tuple(sorted(map(int,base)));seen.add(state);found=False
 print('epoch',epoch,'current',len(base),'elapsed',round(time.time()-start,1),flush=True)
 for trial in range(int(os.getenv('EVAL3_WALK_TRIALS','8000'))):
  seed=int(rng.choice(candidates));vict=masks[seed]
  for _ in range(int(rng.integers(2,14))):
   choices=[int(v) for v in rng.choice(candidates,size=min(400,len(candidates)),replace=False) if 0<(masks[int(v)]|vict).bit_count()-vict.bit_count()<=5 and (masks[int(v)]|vict).bit_count()<=48]
   if not choices:break
   vict|=masks[int(rng.choice(choices))]
  rem=vict.bit_count();ids=np.array([int(v) for v in cand if masks[int(v)] & ~vict==0]);
  if len(ids)<rem-1 or len(ids)>350:continue
  D=A[np.ix_(ids,ids)];i,j=np.where(np.triu(D,1));M=coo_matrix((np.ones(2*len(i)),(np.repeat(np.arange(len(i)),2),np.c_[i,j].ravel())),shape=(len(i),len(ids))).tocsr();res=milp(-np.ones(len(ids)),integrality=np.ones(len(ids)),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,1),options={'time_limit':.5})
  if res.fun is None:continue
  chosen=ids[res.x>.5];kept=np.array([int(b) for j,b in enumerate(base) if not (vict>>j&1)]);new=np.r_[kept,chosen];state=tuple(sorted(map(int,new)))
  if len(new)<max(original-1,len(base)-1) or state in seen:continue
  assert not A[np.ix_(new,new)].any();seen.add(state);delta=len(new)-len(base);base=new;found=True;print('MOVE',epoch,'gain',delta,'total',len(base),'removed',rem,'trial',trial,flush=True)
  if len(base)>original:
   np.savez('research/mixed_tier17_order8_20261105/cache/walk_POTENTIAL.npz',chosen=base,kept=np.array([],np.int32));print('POTENTIAL',4322+8*len(base),'requires full exact replay',flush=True)
  break
 if not found:print('no admissible move in sampled neighborhoods',epoch,flush=True);break
