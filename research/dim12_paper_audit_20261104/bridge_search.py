"""Exploit the paper audit: alternative bridges and independently movable 48 blocks."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import argparse,hashlib,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from threadpoolctl import threadpool_limits
from audit import clebsch_coefficients
from alternative_bridges import coefficients as other_coefficients
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'dim12_quaternionic_20261104'))
from search import normalize,metrics,relax
from continuation import optimize
from verify import verify


def initialize(index,rng):
    alt=bool(index%2);a,b=other_coefficients() if alt else clebsch_coefficients()
    bridges=(a[120:]+np.sqrt(2)*b[120:])/12
    fixed=np.vstack([np.eye(12),-np.eye(12),bridges])
    files=sorted((ROOT/'runs'/'archive_valid').glob('*.npy'))
    left=normalize(np.load(files[(index//2)%7]));right=normalize(np.load(files[(index//14)%7]))
    left=left[:,rng.permutation(6)];right=right[:,rng.permutation(6)]
    free=np.vstack([np.pad(left,((0,0),(0,6))),np.pad(right,((0,0),(6,0)))])
    mode=(index//98)%3
    if mode==0:
        extra=np.zeros((2,12));extra[0,:6]=rng.choice([-1,1],size=6)/np.sqrt(6)
        extra[1,6:]=rng.choice([-1,1],size=6)/np.sqrt(6);free=np.vstack([free,extra])
    elif mode==1:
        # Split one point in each R6 factor; all 98 points may leave the factors.
        extra=[]
        for j in [int(rng.integers(48)),int(rng.integers(48,96))]:
            z=free[j].copy();v=rng.normal(size=12);v-=(v@z)*z;v/=np.linalg.norm(v)
            theta=float(rng.choice([.25,.4,.55]));free[j]=np.cos(theta)*z+np.sin(theta)*v
            extra.append(np.cos(theta)*z-np.sin(theta)*v)
        free=np.vstack([free,extra])
    else:
        z=rng.normal(size=12);z/=np.linalg.norm(z);free=np.vstack([free,z,-z])
    return fixed,normalize(free),alt,mode


def trial(spec):
    index,steps,out=spec;out=Path(out);threadpool_limits(1);start=time.monotonic();seed=720842260618984+index
    rng=np.random.default_rng(seed);fixed,free,alternative,mode=initialize(index,rng)
    free,info=relax(free,steps,fixed=fixed)
    x=np.vstack([fixed,free]);frozen=metrics(x);trace=[];best=x.copy();bm=frozen
    for r in range(2):
        x,opt=relax(x,steps*2);m=metrics(x);trace.append(dict(stage='released',round=r,**m,**opt))
        if m['energy']<bm['energy']:best=x.copy();bm=m
        if bm['max_cos']<.5-1e-10:break
        x=best.copy();chosen=rng.choice(np.arange(744,842),size=12,replace=False)
        x[chosen]+=rng.normal(scale=.04,size=(12,12));x=normalize(x)
    x=best.copy();minimax=x.copy();mm=bm
    for power in [4,8,16,32,64]:
        x,opt=optimize(x,600,power=power);m=metrics(x)
        if m['max_cos']<mm['max_cos']:minimax=x.copy();mm=m
    path=out/f'bridge_{index:04d}.npy';np.save(path,minimax)
    np.save(out/f'bridge_{index:04d}_energy.npy',best)
    result=dict(index=index,seed=seed,alternative_bridges=alternative,seed_mode=mode,
                frozen=frozen,best=bm,minimax=mm,trace=trace,elapsed_seconds=time.monotonic()-start)
    if mm['max_cos']<=.5+1e-12:result['exact_certificate']=verify(path)
    (out/f'bridge_{index:04d}.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=294)
    p.add_argument('--workers',type=int,default=32);p.add_argument('--steps',type=int,default=2000)
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'bridge_search');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    cfg=vars(a).copy();cfg['out']=str(a.out.resolve());cfg['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (a.out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n');results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(trial,(i,a.steps,str(a.out.resolve()))) for i in range(a.trials)]):
            r=f.result();results.append(r)
            print(json.dumps({k:r[k] for k in ['index','alternative_bridges','frozen','best','minimax']}),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
