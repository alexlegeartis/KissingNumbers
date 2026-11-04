"""Change all112 published cap orbits jointly while proposing extra head orbits."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from scipy.optimize import minimize
from build19 import build,codes,H
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.numerical_shapes import orbit_shape_objective
OUT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=20261104);ap.add_argument('--extra',type=int,default=1);ap.add_argument('--restarts',type=int,default=2);ap.add_argument('--iterations',type=int,default=700);args=ap.parse_args()
 body,caps,words,g=build(mark=746);E,octads,U0=codes();rng=np.random.default_rng(args.seed)
 signs=np.array([[1-2*(sum((m>>i&1)*row[j] for i,row in enumerate(g))%2) for j in range(19)] for m in range(16)])
 def canonical(row):return min(tuple(map(int,r)) for r in signs*np.array(row))
 def headkey(row):return min(tuple(map(int,r[:16])) for r in signs*np.array(row))
 capreps=sorted({canonical(row) for row in caps});oldheads={headkey(row) for row in caps}
 allheads={headkey([1-2*(c>>j&1) for j in H]+[0,0,0]) for c in E}
 missing=sorted(allheads-{tuple(1 if x>0 else -1 for x in h) for h in oldheads})
 if len(capreps)!=112 or len(missing)!=16:raise ValueError(('orbit count',len(capreps),len(missing)))
 rng.shuffle(missing);initial=[np.array(row,float)/np.linalg.norm(row) for row in capreps]
 for head in missing[:args.extra]:
  tail=rng.normal(size=3);tail*=np.sqrt(2.)/np.linalg.norm(tail)
  row=np.r_[head,tail];initial.append(row/np.linalg.norm(row))
 shapes=np.array(initial);context=RunContext(ROOT,{'max_gram_entries':10000000,'max_rank_work':10000000,'max_materialized_points':20000})
 objective,audit=orbit_shape_objective(g,shapes.tolist(),body,context)
 history=[];best=None;started=time.monotonic()
 for run in range(args.restarts):
  z=shapes.copy()
  if run:z+=rng.normal(0,.02+.02*(run%4),z.shape)
  bounds=[]
  for row in shapes:
   for j,x in enumerate(row):
    bounds.append((.025,.6) if j<16 and x>0 else ((-.6,-.025) if j<16 else (-.8,.8)))
  result=minimize(objective,z.ravel(),jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':args.iterations,'ftol':1e-15,'gtol':1e-9,'maxls':30})
  report=audit(result.x);score=max(report['max_pair_cosine'],report['max_base_cosine'])
  row={k:v for k,v in report.items() if k!='shapes'}
  row.update({'restart':run,'loss':float(result.fun),'iterations':int(result.nit),'success':bool(result.success),'message':str(result.message),'seconds':time.monotonic()-started})
  print(json.dumps(row),flush=True);history.append(row)
  if best is None or score<best[0]:best=(score,report['shapes'])
 result={'seed':args.seed,'extra_head_orbits':args.extra,'added_head_keys':missing[:args.extra],'initial_cap_orbits':len(capreps),'orbit_size':16,'base_points':len(body),'proposed_total':len(body)+len(shapes)*16,'history':history,'best_shapes':best[1],'new_sota':False,'verification_status':'numerical','scope':'All112 old cap shapes can change continuously, plus extra missing-head orbits. Body fixed. No exact feasible improvement or exclusion is inferred.'}
 (OUT/f'continuous19-extra{args.extra}-s{args.seed}.json').write_text(json.dumps(result)+'\n')
if __name__=='__main__':main()
