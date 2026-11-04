"""Separate discovery control: direct feasibility optimization from exact 840s.

This is not the authors' Riesz schedule. A strict target and scaled residual
avoid treating a floating-point near-zero objective as a valid 841 code.
"""
from common import *
import argparse
import time
from concurrent.futures import ProcessPoolExecutor,as_completed
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from search import Hinge
sys.path.insert(0,str(AUDIT))
from audit import clebsch_coefficients
from alternative_bridges import coefficients as alternative_coefficients

def scaled_relax(x,steps=8000,target=.5-2e-8):
    fun=Hinge(x.shape,target=target)
    def objective(raw):
        value,grad=fun(raw);return value*1e6,grad*1e6
    out=minimize(objective,normalize(x).ravel(),jac=True,method='L-BFGS-B',
                 options=dict(maxiter=steps,maxfun=steps*3,maxcor=24,ftol=1e-15,gtol=1e-9,maxls=40))
    return normalize(out.x.reshape(x.shape)),dict(nit=int(out.nit),calls=fun.calls,scaled_energy=float(out.fun),message=str(out.message))

def trial(spec):
    index,count,steps,out=spec;threadpool_limits(1);out=Path(out);start=time.monotonic();rng=np.random.default_rng(8412026110500+index)
    alternative=bool(index%2);a,b=alternative_coefficients() if alternative else clebsch_coefficients()
    base=(a+np.sqrt(2)*b)/12;extra=rng.choice([-1.,1.],size=(count-840,12))/np.sqrt(12)
    x=np.vstack([base,extra]);noise=[0.,.0001,.001,.01][(index//2)%4]
    x=normalize(x+rng.normal(scale=noise,size=x.shape));initial=metrics(x)
    x,first=relax(x,steps//2);middle=metrics(x);x,last=scaled_relax(x,steps)
    m=metrics(x);path=out/f'candidate_{index:03d}.npy';np.save(path,x)
    r=dict(index=index,alternative=alternative,noise=noise,initial=initial,middle=middle,final=m,
           first_optimizer=first,last_optimizer=last,elapsed_seconds=time.monotonic()-start)
    if m['max_cos']<=.5+1e-12:r['exact_certificate']=verify(path)
    save_json(out/f'candidate_{index:03d}.json',r);return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=128);p.add_argument('--workers',type=int,default=24)
    p.add_argument('--count',type=int,default=841);p.add_argument('--steps',type=int,default=8000)
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'direct841');a=p.parse_args();reserve();a.out.mkdir(parents=True,exist_ok=True)
    save_json(a.out/'config.json',dict(trials=a.trials,workers=a.workers,count=a.count,steps=a.steps,source_sha256=sha(__file__)))
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(trial,(i,a.count,a.steps,str(a.out.resolve()))) for i in range(a.trials)]):
            r=future.result();results.append(r);save_json(a.out/'summary.json',results)
            print(json.dumps({k:r[k] for k in ['index','alternative','noise','final','elapsed_seconds']}),flush=True)
