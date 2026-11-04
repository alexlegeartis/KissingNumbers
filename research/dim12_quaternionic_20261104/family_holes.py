"""Look for an insertable antipodal pair in the movable exact 840 family."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[key]='1'
import argparse,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from balanced_family import Family
from search import ROOT,normalize,metrics

_family=None


def initialize():
    global _family
    threadpool_limits(1);_family=Family()


def objective(flat):
    raw=flat[:12];norm=np.linalg.norm(raw);z=raw/norm;a=flat[12:]
    core=_family.core(a);x=np.vstack([core,_family.heads]);v=np.maximum(x@z-(.5-1e-8),0)
    value=v@v;gz=2*v@x;ga=np.zeros(3)
    for r in range(3):
        sl=slice(408+128*r,408+128*(r+1));gc=2*v[sl,None]*z
        ga[r]=(np.sum(gc[:,4:8]*core[sl,4:8])/(2*a[r])
               -np.sum(gc[:,8:]*core[sl,8:])/(2*(.75-a[r])))
    gz-=(gz@z)*z
    return float(value),np.r_[gz/norm,ga]


def trial(spec):
    index,out=spec;out=Path(out);rng=np.random.default_rng(42042120261104+index)
    z=rng.normal(size=12);z/=np.linalg.norm(z)
    a=np.ones(3)*.375 if index%3==0 else rng.uniform(.25,.5,size=3)
    opt=minimize(objective,np.r_[z,a],jac=True,method='L-BFGS-B',
                 bounds=[(None,None)]*12+[(.25,.5)]*3,
                 options=dict(maxiter=1200,maxcor=10,ftol=1e-16,gtol=1e-10,maxls=30))
    z=opt.x[:12]/np.linalg.norm(opt.x[:12]);a=opt.x[12:]
    x=np.vstack([_family.core(a),_family.heads]);peak=float(np.max(x@z))
    result=dict(index=index,parameters=a.tolist(),direction=z.tolist(),max_cos_to_family=peak,
                energy=float(opt.fun),iterations=opt.nit)
    if peak<.5:
        np.save(out/f'pair_candidate_{index:05d}.npy',np.vstack([x,z,-z]))
        (out/f'pair_candidate_{index:05d}.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=8192)
    p.add_argument('--workers',type=int,default=16);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True);start=time.monotonic();results=[]
    with ProcessPoolExecutor(max_workers=a.workers,initializer=initialize) as pool:
        for f in as_completed([pool.submit(trial,(i,str(a.out.resolve()))) for i in range(a.trials)]):
            r=f.result();results.append(r)
            if len(results)%256==0:
                best=min(results,key=lambda x:x['max_cos_to_family'])
                print(json.dumps(dict(completed=len(results),best=best)),flush=True)
                (a.out/'progress.json').write_text(json.dumps(dict(completed=len(results),best=best),indent=2)+'\n')
    report=dict(trials=a.trials,workers=a.workers,elapsed_seconds=time.monotonic()-start,
                insertable_pairs_found=sum(r['max_cos_to_family']<.5 for r in results),
                best=min(results,key=lambda x:x['max_cos_to_family']))
    (a.out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
