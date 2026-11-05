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
# Neighbor masks use arbitrary-precision integers: 176 incumbent orbit bits.
m=[]
for row in A[:,base]:
 q=0
 for j in np.flatnonzero(row):q|=1<<int(j)
 m.append(q)
isbase=np.zeros(P,bool);isbase[base]=1;cand=np.flatnonzero(~isbase);cost=np.array([q.bit_count() for q in m]);candidates=cand[cost[cand]<=20];print('seed',len(base),'candidate',len(candidates),'min blockers',cost[candidates].min(),flush=True)
rng=np.random.default_rng(170819);best=-999;seen=set();solved=0;t=time.time()
for trial in range(int(os.getenv('EVAL3_TRIALS','30000'))):
 seed=int(rng.choice(candidates));vict=m[seed]
 for _ in range(int(rng.integers(1,10))):
  options=[int(v) for v in rng.choice(candidates,size=min(350,len(candidates)),replace=False) if 0<(m[int(v)]|vict).bit_count()-vict.bit_count()<=4 and (m[int(v)]|vict).bit_count()<=36]
  if not options:break
  vict|=m[int(rng.choice(options))]
 if vict in seen:continue
 seen.add(vict);removed=vict.bit_count();ids=np.array([int(v) for v in cand if m[int(v)] & ~vict==0]);
 if len(ids)<=removed or len(ids)>300:continue
 D=A[np.ix_(ids,ids)];i,j=np.where(np.triu(D,1));M=coo_matrix((np.ones(2*len(i)),(np.repeat(np.arange(len(i)),2),np.c_[i,j].ravel())),shape=(len(i),len(ids))).tocsr();res=milp(-np.ones(len(ids)),integrality=np.ones(len(ids)),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,1),options={'time_limit':.5});solved+=1
 if res.fun is not None:
  val=int(round(-res.fun));gain=val-removed
  if gain>best:best=gain;print('best gain',gain,'orbits','victims',removed,'added',val,'eligible',len(ids),'trial',trial,'secs',round(time.time()-t,1),flush=True)
  if gain>0:
   chosen=ids[res.x>.5];assert len(chosen)-removed==gain;kept=np.array([int(b) for j,b in enumerate(base) if not (vict>>j&1)]);assert not A[np.ix_(np.r_[chosen,kept],np.r_[chosen,kept])].any();np.savez('research/mixed_tier17_order8_20261105/cache/POTENTIAL.npz',chosen=chosen,kept=kept,gain=gain);print('POTENTIAL',5730+8*gain,'requires independent full-coordinate verification',flush=True);break
 if trial%1000==0:print('trial',trial,'solved',solved,'best',best,'secs',round(time.time()-t,1),flush=True)
print('DONE best',best,'solved',solved,flush=True)
