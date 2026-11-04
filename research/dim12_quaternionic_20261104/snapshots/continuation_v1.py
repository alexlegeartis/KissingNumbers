"""Escape contact basins by changing ambient rank or temporarily decompressing."""
import os
for _key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[_key]='1'

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from search import DATA, ROOT, normalize, metrics, seed_configuration, relax


class PowerHinge:
    def __init__(self, shape, power=2, target=.5, extra=0, penalty=0.):
        self.shape=shape
        self.power=power
        self.target=target
        self.extra=extra
        self.penalty=penalty
        self.calls=0

    def __call__(self, flat):
        self.calls+=1
        raw=flat.reshape(self.shape)
        norms=np.maximum(np.linalg.norm(raw,axis=1,keepdims=True),1e-100)
        x=raw/norms
        violation=np.maximum(x@x.T-self.target,0)
        np.fill_diagonal(violation,0)
        peak=violation.max()
        if peak>0:
            scaled=violation/peak
            s=.5*np.sum(scaled**self.power)
            value=peak*peak*s**(2/self.power)
            weights=2*peak*s**(2/self.power-1)*scaled**(self.power-1)
            grad=weights@x
        else:
            value=0.
            grad=np.zeros_like(x)
        if self.extra:
            value+=self.penalty*np.sum(x[:,-self.extra:]**2)
            grad[:,-self.extra:]+=2*self.penalty*x[:,-self.extra:]
        grad-=np.sum(grad*x,axis=1,keepdims=True)*x
        return float(value),(grad/norms).ravel()


def optimize(x, steps, power=2, target=.5, extra=0, penalty=0.):
    fun=PowerHinge(x.shape,power,target,extra,penalty)
    opt=minimize(fun,x.ravel(),jac=True,method='L-BFGS-B',options=dict(
        maxiter=steps,maxfun=3*steps,maxcor=12,ftol=1e-16,gtol=1e-10,maxls=30))
    return normalize(opt.x.reshape(x.shape)),dict(nit=opt.nit,calls=fun.calls,
                                                  objective=float(opt.fun))


def rotate_pca(x):
    _,v=np.linalg.eigh(x.T@x)
    return x@v[:,::-1]


def trial(spec):
    threadpool_limits(1)
    method,source,index,seed,steps,cycles,out=spec
    out=Path(out)
    rng=np.random.default_rng(seed)
    start=time.monotonic()
    if source=='published842':
        x=np.load(DATA/'published842.npy').copy()
    else:
        x,_=seed_configuration(source,rng,index)
        x,_=relax(x,steps)
    initial=metrics(x)
    best=x.copy()
    bm=initial
    trace=[]
    prefix=f'{method}_{source}_{index:04d}'
    for cycle in range(cycles):
        x=best.copy()
        if method=='lift':
            extra=[1,2,4][index%3]
            z=rng.normal(size=(len(x),extra))
            if index%4==0 and source in ['quaternion26','record26']:
                z[:816]*=.1
            scale=[.05,.10,.20,.35][(index//3)%4]
            x=normalize(np.column_stack([x,scale*z]))
            for penalty in [0.,1e-6,1e-5,1e-4,1e-3,.01,.1,1.]:
                x=rotate_pca(x)
                x,opt=optimize(x,steps,extra=extra,penalty=penalty)
                trace.append(dict(cycle=cycle,stage='rank',extra=extra,penalty=penalty,
                                  extra_mass=float(np.sum(x[:,-extra:]**2)),**metrics(x),**opt))
            x=normalize(rotate_pca(x)[:,:12])
        elif method=='decompress':
            maximum=[.505,.51,.52,.54][index%4]
            scale=[.005,.02,.05,.10][(index//4)%4]
            x=normalize(x+rng.normal(scale=scale,size=x.shape))
            for target in np.linspace(maximum,.5,7):
                x,opt=optimize(x,steps,target=float(target))
                trace.append(dict(cycle=cycle,stage='decompress',target=float(target),**metrics(x),**opt))
        else:
            raise ValueError(method)
        x,opt=optimize(x,steps*2)
        m=metrics(x)
        trace.append(dict(cycle=cycle,stage='12d',**m,**opt))
        # Energy controls feasibility. Preserve a separate minimax checkpoint below.
        if m['energy']<bm['energy']:
            best=x.copy();bm=m
            np.save(out/f'{prefix}.npy',best)
        if m['max_cos']<.5-1e-10:
            break
    x=best.copy()
    mm=metrics(x)
    minimax=x.copy()
    for power in [4,8,16,32]:
        x,opt=optimize(x,steps,power=power)
        m=metrics(x)
        trace.append(dict(stage='minimax',power=power,**m,**opt))
        if m['max_cos']<mm['max_cos']:
            minimax=x.copy();mm=m
    np.save(out/f'{prefix}.npy',best)
    np.save(out/f'{prefix}_minimax.npy',minimax)
    result=dict(method=method,source=source,index=index,seed=seed,steps=steps,cycles=cycles,
                initial=initial,best=bm,minimax=mm,trace=trace,
                elapsed_seconds=time.monotonic()-start)
    (out/f'{prefix}.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--methods',nargs='+',default=['lift','decompress'])
    p.add_argument('--sources',nargs='+',default=['quaternion26','record26','published842','prune844'])
    p.add_argument('--trials',type=int,default=24)
    p.add_argument('--workers',type=int,default=64)
    p.add_argument('--steps',type=int,default=500)
    p.add_argument('--cycles',type=int,default=2)
    p.add_argument('--seed',type=int,default=42120261104)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    config=vars(a).copy();config['out']=str(a.out.resolve())
    config['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    config['search_sha256']=hashlib.sha256((ROOT/'search.py').read_bytes()).hexdigest()
    (a.out/'config.json').write_text(json.dumps(config,indent=2)+'\n')
    specs=[(method,source,i,a.seed+mi*1000000+si*10000+i,a.steps,a.cycles,str(a.out.resolve()))
           for i in range(a.trials) for mi,method in enumerate(a.methods) for si,source in enumerate(a.sources)]
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(trial,s) for s in specs]):
            r=future.result();results.append(r)
            print(json.dumps(dict(method=r['method'],source=r['source'],index=r['index'],
                                  best=r['best'],minimax=r['minimax'],
                                  seconds=round(r['elapsed_seconds'],2))),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':
    main()
