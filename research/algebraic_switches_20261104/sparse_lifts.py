"""Break signed support orbits into individual six-core/two-tail choices.

The core oracle is exact incidence/sign arithmetic; all remaining blockers
are explicit coordinate products. Every incumbent point remains removable.
"""
from pathlib import Path
from itertools import combinations,product
import argparse,hashlib,importlib.util,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('dimension',type=int);ap.add_argument('--blockers',type=int,default=2)
 ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=120)
 ap.add_argument('--source',type=Path);ap.add_argument('--core-weight',type=int,choices=(5,6),default=6)
 ap.add_argument('--body-limit',type=int);args=ap.parse_args();start=time.monotonic();weight=args.core_weight
 if args.body_limit is not None and args.body_limit<=0:raise ValueError('Positive body limit')
 args.out.mkdir(parents=True,exist_ok=False)
 source=args.source or next((ROOT/'verifications/improved/dim27-31-joint-refills/data').glob(f'coordinates{args.dimension}_*.npz'))
 X=np.load(source,allow_pickle=False)['points'].astype(np.int16);n,d=X.shape
 if d!=args.dimension or not np.all(np.sum(X.astype(np.int64)**2,axis=1)==128):raise ValueError('Source norms/dimension')
 shape=np.all((X==0)|(abs(X)==4),axis=1)&((X[:,:24]!=0).sum(axis=1)==6)&((X[:,24:]!=0).sum(axis=1)==2)
 six_supports={tuple(map(int,np.flatnonzero(row[:24]))) for row in X[shape]}
 supports=sorted({part for support in six_supports for part in combinations(support,weight)})
 core_ids=np.flatnonzero(~np.any(X[:,24:],axis=1));aux_ids=np.flatnonzero(np.any(X[:,24:],axis=1))
 core=X[core_ids,:24];octads=np.all((core==0)|(abs(core)==4),axis=1)&((core!=0).sum(axis=1)==8)
 other=core[~octads]
 if np.any(np.sort(abs(other),axis=1)[:,-weight:].sum(axis=1)>16):raise ValueError('Uncovered core shape')
 core_octads=core[octads];octad_ids=core_ids[octads]
 masks=np.sum((core_octads!=0).astype(np.uint32)*(np.uint32(1)<<np.arange(24,dtype=np.uint32)),axis=1,dtype=np.uint32)
 signs=np.sum((core_octads<0).astype(np.uint32)*(np.uint32(1)<<np.arange(24,dtype=np.uint32)),axis=1,dtype=np.uint32)
 bodies=[];body_blockers=[]
 for support in supports:
  mask=np.uint32(sum(1<<j for j in support));intersection=masks&mask
  ids=np.flatnonzero(np.bitwise_count(intersection)>=5)
  for word in range(1<<weight):
   sign=np.uint32(sum(1<<j for t,j in enumerate(support) if word>>t&1))
   blocked=octad_ids[ids[((signs[ids]^sign)&intersection[ids])==0]]
   if len(blocked)>args.blockers:continue
   row=np.zeros(d,np.int16);row[list(support)]=[-4 if word>>t&1 else 4 for t in range(weight)]
   bodies.append(row);body_blockers.append(blocked.tolist())
 admitted_bodies=len(bodies)
 if args.body_limit is not None and len(bodies)>args.body_limit:
  rng=np.random.default_rng(20261104);order=rng.permutation(len(bodies)).tolist()
  order.sort(key=lambda i:len(body_blockers[i]));order=order[:args.body_limit]
  bodies=[bodies[i] for i in order];body_blockers=[body_blockers[i] for i in order]
 print(json.dumps({'supports':len(supports),'core_weight':weight,'signed_body_shapes':len(supports)*(1<<weight),'core_admitted_bodies':admitted_bodies,'tested_bodies':len(bodies),'core_points':len(core),'auxiliary_points':len(aux_ids)}),flush=True)
 known={row.tobytes() for row in X};candidates=[];blockers=[]
 for body,blocked in zip(bodies,body_blockers):
  for tail in combinations(range(24,d),8-weight):
   for tail_signs in product((-4,4),repeat=8-weight):
    row=body.copy();row[list(tail)]=tail_signs
    if row.tobytes() in known:continue
    candidates.append(row);blockers.append(blocked[:])
 Y=np.array(candidates,np.int16).reshape((-1,d));active=np.arange(len(Y));checks=0
 if d*int(abs(X).max(initial=0))*4>=2**24:raise ValueError('Exact accumulation bound')
 rng=np.random.default_rng(20261104);order=rng.permutation(aux_ids)
 for at in range(0,len(order),128):
  if not len(active):break
  ids=order[at:at+128];bad=Y[active].astype(np.float32)@X[ids].astype(np.float32).T>64
  checks+=len(active)*len(ids);counts=np.array([len(blockers[int(i)]) for i in active])+bad.sum(axis=1)
  good=counts<=args.blockers;retained=active[good]
  rows,cols=np.nonzero(bad[good])
  for i,j in zip(rows,cols):blockers[int(retained[i])].append(int(ids[j]))
  for i in active[~good]:blockers[int(i)]=[]
  active=retained
 print(json.dumps({'new_candidates':len(Y),'admitted_candidates':len(active),'dot_products':checks}),flush=True)
 Z=Y[active].astype(np.int64);selected_blockers=[blockers[int(i)] for i in active]
 np.savez_compressed(args.out/'pool.npz',candidates=Z,blocker_indptr=np.r_[0,np.cumsum(list(map(len,selected_blockers)))],
   blocker_indices=np.array([j for row in selected_blockers for j in row],np.int64))
 spec=importlib.util.spec_from_file_location('_sparse_lift_assembly',ROOT/'research/sol_discovery_20261104/contact_campaign.py')
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 result=module.optimize(X.astype(np.int64),np.zeros((n,d),np.int64),Z,np.zeros_like(Z),selected_blockers,128,args.seconds) if len(Z) else {'gain':0,'added':[],'removed':[]}
 if result['gain']>0:
  keep=np.ones(n,bool);keep[result['removed']]=False;final=np.r_[X[keep],Z[result['added']]].astype(np.int16)
  # Exact direct integer replay of every changed inequality, independently
  # of the incidence oracle and the numerical optimization model.
  lookup={row.tobytes():i for i,row in enumerate(final)};maximum=-128
  for row in Z[result['added']]:
   values=np.zeros(len(final),np.int64)
   for j in np.flatnonzero(row):values+=final[:,j].astype(np.int64)*int(row[j])
   values[lookup[row.astype(np.int16).tobytes()]]=-128
   if values.max()>64:raise ValueError('Direct inserted-point angle')
   maximum=max(maximum,int(values.max()))
  result['direct_integer_maximum']=maximum
  np.savez_compressed(args.out/'candidate.npz',points=final,squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
 report={'dimension':d,'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
  'baseline':n,'supports':len(supports),'core_weight':weight,'signed_body_shapes':len(supports)*(1<<weight),'core_admitted_bodies':admitted_bodies,
  'tested_bodies':len(bodies),'body_limit':args.body_limit,'body_generation_complete':len(bodies)==admitted_bodies,
  'candidate_pool':len(Y),'admitted_candidates':len(Z),'blocker_cutoff':args.blockers,'dot_products':checks,
  'optimization':result,'seconds':time.monotonic()-start,'new_sota':False,
  'scope':'Seeded core/tail signed-support grammar with stated weight and body budget; individually selectable points, all incumbent layers mutable. Full final geometry replay is required for a lower bound.'}
 (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
