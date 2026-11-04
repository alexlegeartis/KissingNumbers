"""Reflected Gram feasibility search with a periodically corrected rank projection.

The angle/diagonal projection is exact. Rank-12 projection uses warm subspace
iteration between full symmetric eigensolves. This is a heuristic optimizer;
only verify.py may certify a resulting point set.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[key]='1'
import argparse,hashlib,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits
from search import ROOT,DATA,normalize,metrics,relax
from continuation import optimize


def rank_projection(matrix,q,exact):
    n=len(matrix)
    if exact or q is None:
        w,q=eigh(matrix,subset_by_index=[n-12,n-1],check_finite=False)
    else:
        # Shift preserves eigenvectors and avoids selecting negative eigenvalues
        # by magnitude. A full solve every 40 iterations corrects the subspace.
        shifted=matrix+np.eye(n)*10
        for _ in range(4):q,_=np.linalg.qr(shifted@q,mode='reduced')
        w,v=np.linalg.eigh(q.T@matrix@q);q=q@v
    factor=q*np.sqrt(np.maximum(w,0))
    return factor@factor.T,q,factor


def trial(spec):
    i,steps,out=spec;out=Path(out);threadpool_limits(1)
    seed=12084220261104+i;rng=np.random.default_rng(seed);start=time.monotonic()
    if i%4==0:x=np.load(DATA/'published842.npy')
    else:
        route=['quaternion26','record26','record_global'][i%4-1]
        records=json.loads((ROOT/'runs'/'ensemble_fast'/'summary.json').read_text())
        r=min((r for r in records if r['route']==route),key=lambda r:r['best']['energy'])
        x=np.load(ROOT/'runs'/'ensemble_fast'/f"{route}_{r['trial']:04d}.npy")
    original=metrics(x);best=x.copy();bm=original;z=x@x.T
    noise=[0.,.001,.01,.04][(i//4)%4]
    a=rng.normal(scale=noise,size=z.shape);z+=(a+a.T)/2
    beta=[.3,.6,.9,1.2][(i//16)%4];q=None;trace=[]
    for step in range(steps):
        a=np.minimum(z,.5);np.fill_diagonal(a,1.)
        b,q,factor=rank_projection(2*a-z,q,step%40==0)
        z+=beta*(b-a)
        if step%20==0:
            candidate=normalize(factor);m=metrics(candidate)
            trace.append(dict(step=step,projection_gap=float(np.linalg.norm(b-a)),**m))
            if m['max_cos']<bm['max_cos']:
                best=candidate.copy();bm=m;np.save(out/f'gram_{i:04d}.npy',best)
            if m['max_cos']<.5-1e-10:break
    # Release all Gram auxiliary variables for a final direct polishing pass.
    x,_=relax(best,1800);energy=metrics(x)
    np.save(out/f'gram_{i:04d}_energy.npy',x)
    for power in [4,8,16,32,64]:
        x,_=optimize(x,600,power=power);m=metrics(x)
        if m['max_cos']<bm['max_cos']:best=x.copy();bm=m
    np.save(out/f'gram_{i:04d}.npy',best)
    r=dict(index=i,seed=seed,steps=steps,beta=beta,gram_noise=noise,
           original=original,best=energy,minimax=bm,trace=trace,
           elapsed_seconds=time.monotonic()-start)
    (out/f'gram_{i:04d}.json').write_text(json.dumps(r,indent=2)+'\n')
    return r


def main():
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=32)
    p.add_argument('--workers',type=int,default=32);p.add_argument('--steps',type=int,default=1200)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    cfg=vars(a).copy();cfg['out']=str(a.out.resolve())
    cfg['sources_sha256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['gram_reflection.py','search.py','continuation.py','kernels.py']}
    (a.out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(trial,(i,a.steps,str(a.out.resolve()))) for i in range(a.trials)]):
            r=f.result();results.append(r)
            print(json.dumps({k:r[k] for k in ['index','best','minimax','elapsed_seconds']}),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
