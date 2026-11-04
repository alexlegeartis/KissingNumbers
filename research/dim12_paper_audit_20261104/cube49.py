"""Search 49 cube-constrained R6 points: with 12 axes this would give 61."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import argparse,itertools,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.optimize import minimize
from numba import njit
from threadpoolctl import threadpool_limits
from audit import clebsch_coefficients,FACTORS
from archive_exact import check
ROOT=Path(__file__).resolve().parent


@njit(cache=True,fastmath=True)
def value_gradient(raw,target):
    n,d=raw.shape;x=np.empty_like(raw);norm=np.empty(n);grad=np.zeros_like(raw);value=0.
    for i in range(n):
        norm[i]=np.sqrt(np.sum(raw[i]*raw[i]));x[i]=raw[i]/norm[i]
        for k in range(d):
            v=abs(x[i,k])-target
            if v>0:value+=v*v;grad[i,k]+=2*v*(1 if x[i,k]>0 else -1)
    for i in range(n):
        for j in range(i):
            dot=0.
            for k in range(d):dot+=x[i,k]*x[j,k]
            v=dot-target
            if v>0:
                value+=v*v
                for k in range(d):grad[i,k]+=2*v*x[j,k];grad[j,k]+=2*v*x[i,k]
    for i in range(n):
        radial=np.sum(grad[i]*x[i]);grad[i]=(grad[i]-radial*x[i])/norm[i]
    return value,grad


def trial(spec):
    i,steps,out=spec;out=Path(out);threadpool_limits(1);rng=np.random.default_rng(496120261104+i);start=time.monotonic()
    valid=sorted((ROOT/'runs'/'archive_valid').glob('*.npy'))
    if i%8<7:x=np.load(valid[i%7])
    else:
        a,b=clebsch_coefficients();x=(a[:48,:6]+np.sqrt(2)*b[:48,:6])/12
    x=x/np.linalg.norm(x,axis=1,keepdims=True)
    extra=np.array(list(itertools.product([-1,1],repeat=6)))[(i//8)%64]/np.sqrt(6)
    x=np.vstack([x,extra]);x+=rng.normal(scale=[0.,.003,.015,.05][(i//64)%4],size=x.shape)
    best=None;energy=float('inf');trace=[]
    def fun(flat,target=.5):
        v,g=value_gradient(flat.reshape(49,6),target);return float(v),g.ravel()
    for r in range(3):
        opt=minimize(fun,x.ravel(),jac=True,method='L-BFGS-B',options=dict(maxiter=steps,maxcor=12,ftol=1e-16,gtol=1e-11,maxls=30))
        x=opt.x.reshape(49,6);x/=np.linalg.norm(x,axis=1,keepdims=True)
        trace.append(float(opt.fun))
        if opt.fun<energy:energy=float(opt.fun);best=x.copy()
        if energy<1e-20:break
        x=best.copy();chosen=rng.choice(49,size=5,replace=False);x[chosen]+=rng.normal(scale=.04,size=(5,6))
    g=best@best.T;np.fill_diagonal(g,-np.inf)
    result=dict(index=i,energy=energy,max_cos=float(g.max()),max_coordinate=float(abs(best).max()),
                trace=trace,elapsed_seconds=time.monotonic()-start,certified61=False)
    if energy<1e-15:
        # Aim for positive margin before attempting the exact predicate.
        opt=minimize(lambda z:fun(z,.5-1e-9),best.ravel(),jac=True,method='L-BFGS-B',
                     options=dict(maxiter=10000,ftol=1e-20,gtol=1e-12,maxls=40))
        best=opt.x.reshape(49,6)
        certificate=check(best);result['head_certificate']=certificate
        if certificate['passed']:
            head=np.vstack([best,np.eye(6),-np.eye(6)]);bridges=[]
            for matching in FACTORS:
                for (a,b),(c,d) in itertools.product(matching,repeat=2):
                    for s in itertools.product([-1,1],repeat=4):
                        z=np.zeros(12);z[[a,b,c+6,d+6]]=s;bridges.append(z)
            whole=np.vstack([np.pad(head,((0,0),(0,6))),np.pad(head,((0,0),(6,0))),bridges])
            path=out/f'842_candidate_{i:04d}.npy';np.save(path,whole)
            sys.path.insert(0,str(ROOT.parent/'dim12_quaternionic_20261104'))
            from verify import verify
            result['whole_certificate']=verify(path);assert result['whole_certificate']['passed'];result['certified61']=True
    np.save(out/f'head49_{i:04d}.npy',best)
    (out/f'head49_{i:04d}.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=1024)
    p.add_argument('--workers',type=int,default=24);p.add_argument('--steps',type=int,default=3500)
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'cube49');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    value_gradient(np.ones((4,6)),.5);results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(trial,(i,a.steps,str(a.out.resolve()))) for i in range(a.trials)]):
            r=f.result();results.append(r)
            if len(results)%32==0 or r['certified61']:
                print(json.dumps(dict(completed=len(results),best=min(results,key=lambda x:x['energy']))),flush=True)
                (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
    (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
