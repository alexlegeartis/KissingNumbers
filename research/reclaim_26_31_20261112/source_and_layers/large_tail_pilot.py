"""Joint coherent-fiber pilot for previously absent tail radii in dimension31.

Uses exact bounded integer dot products. The default finite grammar uses core
supports in the normal octad. All incumbent points are removable. Negative
outcomes apply only to the saved finite pool and stated blocker cap.
"""
import argparse,itertools,json,time,subprocess,sys
from pathlib import Path
import numpy as np
from source_paths import PACKINGSTAR_DEFAULT, packingstar_file
from scipy.sparse import coo_matrix, save_npz
from scipy.optimize import milp,Bounds,LinearConstraint

ap=argparse.ArgumentParser();ap.add_argument('--core',choices=['octad','all'],default='octad');ap.add_argument('--cap',type=int,default=96);ap.add_argument('--solver-seconds',type=float,default=90);ap.add_argument('--max-groups',type=int,default=4000);ap.add_argument('--packingstar',type=Path,default=PACKINGSTAR_DEFAULT);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--wall-seconds',type=float,default=120);args=ap.parse_args()
out=args.output;out.mkdir(parents=True,exist_ok=True)
f=np.load(packingstar_file(args.packingstar,31),allow_pickle=False);den=int(f['coordinate_denominator']);X=(f['points']*(2//den)).astype(np.int16);n,d=X.shape
assert np.all((X.astype(np.int32)**2).sum(1)==128)
core=list(range(24)) if args.core=='all' else [0,2,4,5,6,10,11,23]
t0=time.time();blocks=[];meta=[];sizes=[]
for r in [3,2,1]:
 for ts in itertools.combinations(range(24,31),8-r):
  for st in itertools.product([-4,4],repeat=8-r):
   for hs in itertools.combinations(core,r):
    if r==3:
     hsets=[[(a,b,c) for a,b,c in itertools.product([-4,4],repeat=3) if (a*b*c>0)==parity] for parity in [False,True]]
    elif r==2:
     hsets=[[(-4,-4),(4,4)],[(-4,4),(4,-4)]]
    else:hsets=[[(-4,)],[(4,)]]
    for hsigns in hsets:
     z=np.zeros((len(hsigns),31),np.int16);z[:,ts]=st;z[:,hs]=hsigns
     assert np.all((z.astype(np.int32)**2).sum(1)==128)
     blocks.append(z);sizes.append(len(z));meta.append((r,ts,st,hs,hsigns))
Y=np.concatenate(blocks);sizes=np.array(sizes,np.int32);starts=np.r_[0,np.cumsum(sizes)];owner=np.repeat(np.arange(len(blocks)),sizes)
del blocks
np.savez_compressed(out/'grammar.npz',points=Y,group_starts=starts,group_sizes=sizes)
print(json.dumps({'stage':'generated','groups':len(sizes),'points':len(Y),'seconds':time.time()-t0}),flush=True)
# Every candidate clears all Leech-core points, pure D7 tails, and flat
# norm96 core midpoint bodies. Check sufficient l1 gates against actual data.
iscore=np.all(X[:,24:]==0,axis=1)
ispure=np.all(X[:,:24]==0,axis=1)
isflat=np.all(np.abs(X[:,:24])==2,axis=1)&((X[:,24:].astype(np.int32)**2).sum(1)==32)
assert np.max(np.sort(np.abs(X[iscore,:24]),axis=1)[:,-3:].sum(1))*4<=64
assert np.max(np.abs(X[ispure,24:]).sum(1))*4<=64
assert np.all((X[isflat,24:]!=0).sum(1)==2)
oldids=np.flatnonzero(~(iscore|ispure|isflat))
# Candidates vanish outside these core coordinates, so identical projections
# have identical blockers. Retain multiplicities as exact deletion weights.
cols=core+list(range(24,31));proj,pinverse,pweights=np.unique(X[oldids][:,cols],axis=0,return_inverse=True,return_counts=True)
P=np.zeros((len(proj),31),np.int16);P[:,cols]=proj
rng=np.random.default_rng(20261112);order=rng.permutation(len(P))
np.savez_compressed(out/'profiles.npz',points=P,weights=pweights,old_indices=oldids,profile_indices=pinverse)
active=np.arange(len(Y));cost=np.zeros(len(Y),np.int32);live=np.ones(len(sizes),bool)
for at in range(0,len(order),256):
 if len(active)==0:break
 ids=order[at:at+256]
 for off in range(0,len(active),8192):
  ii=active[off:off+8192]
  bad=Y[ii].astype(np.float32)@P[ids].astype(np.float32).T>64
  cost[ii]+=(bad*pweights[ids]).sum(1,dtype=np.int32)
 dead=owner[active[cost[active]>args.cap]];live[dead]=False;active=active[live[owner[active]]]
 if at%1024==0:print(json.dumps({'stage':'gate','profiles_checked':at+len(ids),'profiles_total':len(P),'old_total':len(oldids),'active_points':len(active),'live_groups':int(live.sum()),'seconds':time.time()-t0}),flush=True)
gids=np.flatnonzero(live);admitted=[];blockers=[];pf=P.astype(np.float32)
for at in range(0,len(gids),256):
 gg=gids[at:at+256];z=np.concatenate([Y[starts[g]:starts[g+1]] for g in gg]);ss=np.r_[0,np.cumsum(sizes[gg])]
 bad=z.astype(np.float32)@pf.T>64
 for j,gid in enumerate(gg):
  blocked=np.flatnonzero(bad[ss[j]:ss[j+1]].any(axis=0))
  if pweights[blocked].sum()<=args.cap:admitted.append(gid);blockers.append(blocked)
report={'source_points':n,'core_scope':args.core,'total_groups':len(sizes),'total_points':len(Y),'blocker_cap':args.cap,'admitted_groups':len(admitted),'relevant_old_points':len(oldids),'weighted_profiles':len(P),'seconds':time.time()-t0,'proof_status':'search pool only'}
print(json.dumps({'stage':'repriced',**report}),flush=True)
gids=np.array(admitted,np.int32);indptr=np.r_[0,np.cumsum(list(map(len,blockers)))];indices=np.concatenate(blockers) if blockers else np.empty(0,np.int32)
np.savez_compressed(out/'pool.npz',points=Y,group_starts=starts,group_sizes=sizes,admitted_group_ids=gids,blocker_indptr=indptr,blocker_indices=indices)
if len(gids)==0:
 (out/'report.json').write_text(json.dumps(report,indent=2));raise SystemExit
# Keep a few alternative core fibers for every tail direction before adding
# more alternatives for one tail, so the finite model retains all radii.
buckets={}
for i,gid in enumerate(gids):
 r,ts,st,hs,hsigns=meta[int(gid)];buckets.setdefault((r,ts,st),[]).append(i)
for bucket in buckets.values():
 bucket.sort(key=lambda i:(pweights[blockers[i]].sum()/sizes[gids[i]],pweights[blockers[i]].sum(),int(gids[i])))
order=[];rank=0
while len(order)<min(args.max_groups,len(gids)):
 available=[bucket[rank] for bucket in buckets.values() if len(bucket)>rank]
 if not available:break
 available.sort(key=lambda i:(pweights[blockers[i]].sum()/sizes[gids[i]],pweights[blockers[i]].sum(),int(gids[i])))
 order.extend(available[:args.max_groups-len(order)]);rank+=1
gids=gids[order];blockers=[blockers[i] for i in order];ng=len(gids)
Z=np.concatenate([Y[starts[g]:starts[g+1]] for g in gids]);zowner=np.repeat(np.arange(ng),sizes[gids])
edges=set()
for at in range(0,len(Z),256):
 bad=Z[at:at+256].astype(np.float32)@Z.astype(np.float32).T>64
 ii,jj=np.nonzero(bad)
 for i,j in zip(zowner[at+ii],zowner[jj]):
  if i<j:edges.add((int(i),int(j)))
victims=np.unique(np.concatenate(blockers));lookup={int(v):i for i,v in enumerate(victims)};rr=[];cc=[];vv=[];ub=[]
for i,blocked in enumerate(blockers):
 for v in blocked:
  row=len(ub);rr.extend((row,row));cc.extend((i,ng+lookup[int(v)]));vv.extend((1,-1));ub.append(0)
for i,j in sorted(edges):
 row=len(ub);rr.extend((row,row));cc.extend((i,j));vv.extend((1,1));ub.append(1)
A=coo_matrix((vv,(rr,cc)),shape=(len(ub),ng+len(victims))).tocsr();c=np.r_[-sizes[gids].astype(float),pweights[victims].astype(float)]
print(json.dumps({'stage':'solve','groups':ng,'victims':len(victims),'edges':len(edges),'seconds':time.time()-t0}),flush=True)
save_npz(out/'model_constraints.npz',A)
np.savez_compressed(out/'model_data.npz',c=c,ub=np.array(ub),group_ids=gids,victims=victims)
report.update(model_groups=ng,model_victims=len(victims),model_cross_edges=len(edges))
try:
 worker=Path(__file__).with_name('solve_model.py')
 completed=subprocess.run([sys.executable,str(worker),'--directory',str(out),'--seconds',str(args.solver_seconds)],timeout=args.wall_seconds,capture_output=True,text=True)
 (out/'solver.log').write_text(completed.stdout+completed.stderr)
 if completed.returncode!=0:raise RuntimeError('Solver subprocess failed; see solver.log')
 solver=json.loads((out/'solver_result.json').read_text());report.update(solver)
 solution=np.load(out/'solver_solution.npy',allow_pickle=False) if solver['has_incumbent'] else None
except subprocess.TimeoutExpired:
 report.update(solver_status='external wall timeout',solver_message='Unresolved: solver returned no result before the external wall limit',has_incumbent=False)
 solution=None
report['seconds']=time.time()-t0
if solution is not None:
 chosen=np.flatnonzero(solution[:ng]>.5);removed_profiles=np.unique(np.concatenate([blockers[i] for i in chosen])) if len(chosen) else np.empty(0,int)
 remove=oldids[np.isin(pinverse,removed_profiles)]
 added=np.concatenate([Y[starts[gids[i]]:starts[gids[i]+1]] for i in chosen]) if len(chosen) else np.empty((0,31),np.int16)
 gain=len(added)-len(remove);report.update(gain=int(gain),added=len(added),removed=len(remove))
 if gain>0:
  keep=np.ones(n,bool);keep[remove]=False;points=np.concatenate([X[keep],added]);assert len(np.unique(points,axis=0))==len(points)
  # Direct int64 changed-pair verification, independent of the support gates.
  for i,p in enumerate(added):
   dots=np.sum(points.astype(np.int64)*p.astype(np.int64),axis=1);dots[n-len(remove)+i]=-128
   assert dots.max()<=64
  np.savez_compressed(out/f'candidate31_{len(points)}.npz',points=points,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
  report['proof_status']='exact changed-pair proof passes; full final geometry still required'
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
