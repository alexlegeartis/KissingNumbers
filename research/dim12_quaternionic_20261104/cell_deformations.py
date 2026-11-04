"""Try two new points from each alternative quaternionic 840-point seed."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,json,time
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from search import ROOT,normalize,metrics,relax
from continuation import optimize


def trial(spec):
    index,seed,steps,rounds,out=spec
    threadpool_limits(1)
    rng=np.random.default_rng(seed);start=time.monotonic();out=Path(out)
    source=ROOT/'runs'/'quaternion_cells'/f'code_{index%256:04d}.npy'
    old=np.load(source)
    mass=np.array([np.sum(old[:,4*j:4*j+4]**2) for j in range(3)])
    factor=int(np.argmin(mass))
    exceptional=np.arange(24*factor,24*(factor+1))
    core=np.delete(old,exceptional,axis=0)
    free=np.zeros((26,12));free[:24]=old[exceptional]
    if index%4==2:
        q=normalize(rng.normal(size=(26,4)))
        q,_=relax(q,1000)
        free[:,4*factor:4*factor+4]=q
    else:
        free[24:,4*factor:4*factor+4]=normalize(rng.normal(size=(2,4)))
    free+=rng.normal(scale=[.01,.03,.06,.12][index%4],size=free.shape)
    free,_=relax(normalize(free),steps,fixed=core)
    x=np.vstack([core,free]);best=x.copy();bm=metrics(x);trace=[]
    for r in range(rounds):
        x,opt=relax(x,steps)
        m=metrics(x);trace.append(dict(round=r,**m,**opt))
        if m['energy']<bm['energy']:
            best=x.copy();bm=m;np.save(out/f'cells_{index:04d}.npy',best)
        # Reconfigure the 26-point region with nearby core points released.
        x=best.copy()
        chosen=rng.choice(np.arange(816,842),size=int(rng.choice([3,6,12,26])),replace=False)
        scale=float(rng.choice([.03,.08,.15]))
        x[chosen]+=rng.normal(scale=scale,size=(len(chosen),12))
        x=normalize(x)
    x=best.copy();mm=metrics(x);minimax=x.copy()
    for p in [4,8,16,32]:
        x,opt=optimize(x,700,power=p);m=metrics(x)
        trace.append(dict(stage='minimax',power=p,**m,**opt))
        if m['max_cos']<mm['max_cos']:minimax=x.copy();mm=m
    np.save(out/f'cells_{index:04d}.npy',best)
    np.save(out/f'cells_{index:04d}_minimax.npy',minimax)
    result=dict(index=index,seed=seed,steps=steps,rounds=rounds,source=str(source.relative_to(ROOT)),
                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                block_mass=mass.tolist(),free_factor=factor,best=bm,minimax=mm,trace=trace,
                elapsed_seconds=time.monotonic()-start)
    (out/f'cells_{index:04d}.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=128)
    p.add_argument('--workers',type=int,default=16);p.add_argument('--steps',type=int,default=3000)
    p.add_argument('--rounds',type=int,default=3);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    cfg=vars(a).copy();cfg['out']=str(a.out.resolve())
    cfg['sources_sha256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['cell_deformations.py','search.py','continuation.py','kernels.py']}
    (a.out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    specs=[(i,84220261104+i,a.steps,a.rounds,str(a.out.resolve())) for i in range(a.trials)]
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(trial,s) for s in specs]):
            r=f.result();results.append(r)
            print(json.dumps({k:r[k] for k in ['index','block_mass','best','minimax','elapsed_seconds']}),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
