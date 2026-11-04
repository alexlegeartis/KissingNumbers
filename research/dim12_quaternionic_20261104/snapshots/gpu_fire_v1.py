"""Bounded RTX A4000 diameter annealing, followed by double-precision polish.

FIRE trajectories use float32 only for exploration. Point diameters return to
one during the final 40 percent of GPU time. All acceptance uses verify.py.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[key]='1'
import argparse,hashlib,json,multiprocessing,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from threadpoolctl import threadpool_limits
from search import ROOT,DATA,normalize,metrics,relax
from continuation import optimize


def initial_batch(count):
    rng=np.random.default_rng(12684220261104);sources=[]
    for batch in ['ensemble_fast','continuation_fast','balanced_family','orbit_refine']:
        report=json.loads((ROOT/'runs'/batch/'summary.json').read_text())
        report=sorted(report,key=lambda r:r['best']['energy'])[:24]
        for r in report:
            if batch=='ensemble_fast':stem=f"{r['route']}_{r['trial']:04d}"
            elif batch=='continuation_fast':stem=f"{r['method']}_{r['source']}_{r['index']:04d}"
            elif batch=='balanced_family':stem=f"balanced_{r['index']:04d}"
            else:stem=f"pool_{r['index']:04d}"
            sources.append(ROOT/'runs'/batch/(stem+'.npy'))
    result=[];provenance=[]
    for i in range(count):
        path=DATA/'published842.npy' if i%8==0 else sources[i%len(sources)]
        x=np.load(path)
        if i>=len(sources):x=normalize(x+rng.normal(scale=.004,size=x.shape))
        result.append(x);provenance.append(str(path.relative_to(ROOT)))
    return np.array(result),provenance


def polish(spec):
    i,path,out=spec;threadpool_limits(1);out=Path(out);start=time.monotonic()
    x=normalize(np.load(path)[i]);initial=metrics(x)
    x,info=relax(x,5000);best=x.copy();bm=metrics(x);trace=[dict(stage='hinge',**bm,**info)]
    mm=bm;minimax=x.copy()
    for power in [4,8,16,32,64]:
        x,info=optimize(x,600,power=power);m=metrics(x);trace.append(dict(stage='minimax',power=power,**m,**info))
        if m['max_cos']<mm['max_cos']:minimax=x.copy();mm=m
    np.save(out/f'gpu_{i:04d}.npy',best);np.save(out/f'gpu_{i:04d}_minimax.npy',minimax)
    r=dict(index=i,initial=initial,best=bm,minimax=mm,trace=trace,elapsed_seconds=time.monotonic()-start)
    (out/f'gpu_{i:04d}.json').write_text(json.dumps(r,indent=2)+'\n');return r


def main():
    p=argparse.ArgumentParser();p.add_argument('--batch',type=int,default=128)
    p.add_argument('--seconds',type=float,default=120);p.add_argument('--workers',type=int,default=24)
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'gpu_fire')
    p.add_argument('--polish-out',type=Path,default=ROOT/'runs'/'gpu_polish')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);a.polish_out.mkdir(parents=True,exist_ok=True)
    if 127 in os.sched_getaffinity(0):raise RuntimeError('Reserve CPU 127 with taskset -c 0-126.')
    import cupy as cp
    cp.cuda.Device(0).use();device=cp.cuda.runtime.getDeviceProperties(0)['name'].decode()
    if device!='NVIDIA RTX A4000':raise RuntimeError(f'Expected RTX A4000 at device 0; got {device}')
    cfg=dict(batch=a.batch,seconds=a.seconds,cpu_polish_workers=a.workers,device=device,
             cupy=cp.__version__,cpu_affinity=sorted(os.sched_getaffinity(0)),
             method='FIRE with mean-zero diameter annealing, equal diameters for final 40 percent',
             seed=12684220261104,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (a.out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (a.polish_out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    initial,provenance=initial_batch(a.batch);rng=np.random.default_rng(12684220261104)
    noise=np.clip(rng.normal(size=(a.batch,842)),-2,2);noise-=noise.mean(axis=1,keepdims=True)
    scales=np.array([.01,.02,.04,.06])[np.arange(a.batch)%4,None]
    noise=cp.asarray((noise*scales).astype(np.float32))
    x=cp.asarray(initial.astype(np.float32));v=cp.zeros_like(x)
    dt=cp.full((a.batch,1,1),.03,dtype=cp.float32);alpha=cp.full_like(dt,.1)
    positive=cp.zeros((a.batch,1,1),dtype=cp.int32)
    diagonal=cp.arange(842)
    overlap=cp.ElementwiseKernel('float32 dot, float32 ri, float32 rj','float32 out',
        'float s=ri+rj; out=fmaxf(dot-(1.0f-s*s*0.125f),0.0f);','diameter_overlap')
    best=initial.copy();best_energy=np.array([metrics(t)['energy'] for t in best]);trace=[]
    # Independent CPU comparison of the CUDA force before the long batch.
    control=cp.asarray(initial[:1].astype(np.float32));cg=control@cp.swapaxes(control,1,2)
    cg=cp.maximum(cg-.5,0);cg[:,diagonal,diagonal]=0;cf=2*cg@control
    host=cp.asnumpy(control[0]).astype(np.float64);gg=np.maximum(host@host.T-.5,0);np.fill_diagonal(gg,0)
    error=float(np.max(abs(cp.asnumpy(cf[0])-2*gg@host)))
    if error>2e-5:raise RuntimeError(f'CUDA force control failed: {error}')
    cp.cuda.Stream.null.synchronize();start=time.monotonic();step=0
    while True:
        elapsed=time.monotonic()-start
        if elapsed>=a.seconds:break
        amplitude=max(0.,1-elapsed/(a.seconds*.6))
        radii=1+amplitude*noise
        g=x@cp.swapaxes(x,1,2);overlap(g,radii[:,:,None],radii[:,None,:],g);g[:,diagonal,diagonal]=0
        force=-2*g@x;force-=cp.sum(force*x,axis=2,keepdims=True)*x
        v+=dt*force;power=cp.sum(v*force,axis=(1,2),keepdims=True);good=power>0
        positive=cp.where(good,positive+1,0)
        dt=cp.where(good,cp.where(positive>5,cp.minimum(dt*1.1,.25),dt),cp.maximum(dt*.5,1e-5))
        alpha=cp.where(good,cp.where(positive>5,alpha*.99,alpha),.1)
        v=cp.where(good,v,0)
        fn=cp.sqrt(cp.sum(force*force,axis=(1,2),keepdims=True));vn=cp.sqrt(cp.sum(v*v,axis=(1,2),keepdims=True))
        v=(1-alpha)*v+alpha*vn/cp.maximum(fn,1e-20)*force
        x+=dt*v;x/=cp.sqrt(cp.sum(x*x,axis=2,keepdims=True));step+=1
        if step%500==0:
            host=cp.asnumpy(x).astype(np.float64)
            m=[metrics(t) for t in host];energy=np.array([t['energy'] for t in m]);improved=energy<best_energy
            best[improved]=host[improved];best_energy=np.minimum(best_energy,energy)
            r=dict(step=step,seconds=time.monotonic()-start,diameter_amplitude=amplitude,
                   best_equal_diameter_energy=float(best_energy.min()),current_best_max_cos=min(t['max_cos'] for t in m))
            trace.append(r);print(json.dumps(r),flush=True)
    host=cp.asnumpy(x).astype(np.float64)
    for i,t in enumerate(host):
        if metrics(t)['energy']<best_energy[i]:best[i]=t
    np.save(a.out/'best_energy_states.npy',best)
    (a.out/'summary.json').write_text(json.dumps(dict(**cfg,steps=step,trace=trace,
        elapsed_seconds=time.monotonic()-start,cuda_force_max_error=error,source_paths=provenance),indent=2)+'\n')
    del x,v,g,force,noise;cp.get_default_memory_pool().free_all_blocks();cp.cuda.Stream.null.synchronize()
    results=[];ctx=multiprocessing.get_context('spawn')
    with ProcessPoolExecutor(max_workers=a.workers,mp_context=ctx) as pool:
        futures=[pool.submit(polish,(i,str((a.out/'best_energy_states.npy').resolve()),str(a.polish_out.resolve()))) for i in range(a.batch)]
        for future in as_completed(futures):
            r=future.result();results.append(r);print(json.dumps({k:r[k] for k in ['index','best','minimax']}),flush=True)
            (a.polish_out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
