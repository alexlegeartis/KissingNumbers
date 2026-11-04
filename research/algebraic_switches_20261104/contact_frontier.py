"""Generate root sums from the new frontier and unpaired incumbent points."""
from pathlib import Path
import argparse,hashlib,importlib.util,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True)
 ap.add_argument('--previous',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
 ap.add_argument('--blockers',type=int,default=2);ap.add_argument('--seconds',type=float,default=120)
 ap.add_argument('--anchor-limit',type=int,default=64);args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False)
 start=time.monotonic();X=np.load(args.source)['points'].astype(np.int64);previous=np.load(args.previous)['points'].astype(np.int64)
 known={tuple(row) for row in X};old={tuple(row) for row in previous}
 frontier=[i for i,row in enumerate(X) if tuple(row) not in old]
 unpaired=[i for i,row in enumerate(X) if tuple(-row) not in known]
 anchors=list(dict.fromkeys(frontier+unpaired))[:args.anchor_limit]
 candidates=[];origins=[];seen=set(known)
 for i in anchors:
  dots=X@X[i]
  for j in np.flatnonzero(dots==-64):
   row=X[i]+X[j];key=tuple(row)
   if key in seen:continue
   if row@row!=128:raise ValueError('Contact sum norm')
   seen.add(key);candidates.append(row);origins.append([int(i),int(j)])
 Y=np.array(candidates,np.int64).reshape((-1,X.shape[1]));print(json.dumps({'anchors':len(anchors),'new_frontier':len(frontier),'unpaired':len(unpaired),'generated':len(Y)}),flush=True)
 spec=importlib.util.spec_from_file_location('_contact_frontier',ROOT/'research/sol_discovery_20261104/contact_campaign.py')
 m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 zero=np.zeros_like(X);Z,_,bs,ids,hist=m.screen(X,zero,Y,np.zeros_like(Y),128,args.blockers)
 np.savez_compressed(args.out/'pool.npz',candidates=Z,source_candidate_indices=np.array(ids,np.int64),
  blocker_indptr=np.r_[0,np.cumsum(list(map(len,bs)))],blocker_indices=np.array([j for row in bs for j in row],np.int64))
 print(json.dumps({'retained':len(Z),'blocker_histogram':hist}),flush=True)
 result=m.optimize(X,zero,Z,np.zeros_like(Z),bs,128,args.seconds) if len(Z) else {'gain':0,'added':[],'removed':[]}
 if result['gain']>0:
  keep=np.ones(len(X),bool);keep[result['removed']]=False;final=np.r_[X[keep],Z[result['added']]]
  lookup={tuple(row):i for i,row in enumerate(final)}
  for row in Z[result['added']]:
   dots=np.sum(final*row,axis=1);dots[lookup[tuple(row)]]=-128
   if dots.max()>64:raise ValueError('Exact changed angle')
  np.savez_compressed(args.out/'candidate.npz',points=final.astype(np.int16),squared_norm=128,inner_product_limit=64,coordinate_denominator=2)
 np.savez_compressed(args.out/'origins.npz',pairs=np.array(origins,np.int32).reshape((-1,2)))
 report={'source':str(args.source),'source_sha256':hashlib.sha256(args.source.read_bytes()).hexdigest(),
  'baseline':len(X),'anchors':anchors,'generated':len(Y),'candidate_origins':{'file':'origins.npz',
    'sha256':hashlib.sha256((args.out/'origins.npz').read_bytes()).hexdigest(),'count':len(origins)},
  'retained':len(Z),'blocker_cutoff':args.blockers,'histogram':hist,'optimization':result,
  'seconds':time.monotonic()-start,'new_sota':False,'scope':'Finite contact sums from new frontier and unpaired points; all source points mutable. Fresh full output replay required.'}
 (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ('candidate_origins','anchors','histogram')},indent=2),flush=True)
if __name__=='__main__':main()
