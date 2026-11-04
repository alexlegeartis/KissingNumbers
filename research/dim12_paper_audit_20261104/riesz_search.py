"""Paper-specific search: broaden block placement, rotations, and 48-system seeds."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import argparse,hashlib,json,multiprocessing,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy.linalg import expm
from threadpoolctl import threadpool_limits
from audit import clebsch_coefficients,riesz
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent/'dim12_quaternionic_20261104'
sys.path.insert(0,str(OLD))
from search import metrics,normalize,relax
from continuation import optimize


def rotation(rng):
    skew=rng.normal(size=(4,4));skew=skew-skew.T;skew/=np.linalg.norm(skew)
    low=0.;high=1.;bound=3/(2*np.sqrt(2))
    for _ in range(40):
        mid=(low+high)/2;a=expm(skew*mid)
        if abs(a).sum(axis=1).max()<=bound:low=mid
        else:high=mid
    a=expm(skew*low*rng.uniform(.5,1.));assert abs(a).sum(axis=1).max()<=bound+1e-13
    return a


def initialize(batch):
    rng=np.random.default_rng(260618984842);aa,bb=clebsch_coefficients();base=(aa+np.sqrt(2)*bb)/12
    original48=base[:48,:6];bridges=base[120:];valid=sorted((ROOT/'runs'/'archive_valid').glob('*.npy'))
    pairs=[(i,j) for i in range(6) for j in range(6) if i!=j];xs=[];info=[]
    for i in range(batch):
        mode=['control','placement','rotation','archive'][i%4]
        left=original48.copy();right=original48.copy();first,last=0,5
        if mode=='rotation':
            left[:,1:5]=left[:,1:5]@rotation(rng).T;right[:,1:5]=right[:,1:5]@rotation(rng).T
        if mode in ['placement','rotation']:
            first,last=pairs[(i//4)%30];order=[first,*[j for j in range(6) if j not in [first,last]],last]
            placed=np.zeros_like(right);placed[:,order]=right;right=placed
        if mode=='archive':
            left=np.load(valid[(i//4)%7]);right=np.load(valid[(3*(i//4)+1)%7])
            left=normalize(left)[:,rng.permutation(6)];right=normalize(right)[:,rng.permutation(6)]
        l60=np.vstack([left,np.eye(6),-np.eye(6)]);r60=np.vstack([right,np.eye(6),-np.eye(6)])
        core=np.vstack([np.pad(l60,((0,0),(0,6))),np.pad(r60,((0,0),(6,0))),bridges])
        extras=rng.choice([-1.,1.],size=(2,12))/np.sqrt(12)
        while extras[0]@extras[1]>.5:extras[1]=rng.choice([-1.,1.],size=12)/np.sqrt(12)
        if mode in ['rotation','archive']:extras=normalize(extras+rng.normal(scale=.02,size=extras.shape))
        x=normalize(np.vstack([core,extras]));xs.append(x)
        info.append(dict(index=i,mode=mode,right_special_positions=[first,last],
                         initial_core=metrics(core),initial=metrics(x)))
    return np.array(xs),info


def polish(spec):
    i,path,out=spec;threadpool_limits(1);out=Path(out);start=time.monotonic()
    x=np.load(path)[i];x=normalize(x);initial=metrics(x);best=x.copy();bm=initial
    x,opt=relax(x,5500);m=metrics(x)
    if m['energy']<bm['energy']:best=x.copy();bm=m
    mm=m;minimax=x.copy();trace=[]
    for power in [4,8,16,32,64,128]:
        x,opt=optimize(x,800,power=power);m=metrics(x);trace.append(dict(power=power,**m,**opt))
        if m['max_cos']<mm['max_cos']:minimax=x.copy();mm=m
    np.save(out/f'candidate_{i:04d}.npy',minimax);np.save(out/f'energy_{i:04d}.npy',best)
    result=dict(index=i,initial=initial,best=bm,minimax=mm,trace=trace,elapsed_seconds=time.monotonic()-start)
    if mm['max_cos']<=.5+1e-12:
        from verify import verify
        result['exact_certificate']=verify(out/f'candidate_{i:04d}.npy')
    (out/f'candidate_{i:04d}.json').write_text(json.dumps(result,indent=2)+'\n');return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--batch',type=int,default=128)
    p.add_argument('--stage-steps',type=int,default=800);p.add_argument('--workers',type=int,default=24)
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'riesz');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    if 127 in os.sched_getaffinity(0):raise RuntimeError('Reserve CPU 127 using taskset.')
    import cupy as cp
    cp.cuda.Device(0).use();name=cp.cuda.runtime.getDeviceProperties(0)['name'].decode()
    if name!='NVIDIA RTX A4000':raise RuntimeError(name)
    initial,metadata=initialize(a.batch);np.save(a.out/'initial.npy',initial)
    schedule=[(8,.005),(16,.003),(32,.002),(64,.001),(128,.0005),(256,.0002),
              (512,.0001),(1024,.00005),(2048,.00001),(4096,.00001),(16384,.000003),(65536,.000001)]
    cfg=dict(batch=a.batch,stage_steps=a.stage_steps,schedule=schedule,workers=a.workers,device=name,
             seed=260618984842,cpu_affinity=sorted(os.sched_getaffinity(0)),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (a.out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    (a.out/'initial_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    diagonal=cp.arange(842)
    def gpu_gradient(x,s):
        dots=x@cp.swapaxes(x,1,2);dist=cp.maximum(2-2*dots,1e-12)
        logterms=-.5*s*cp.log(dist);logterms[:,diagonal[:x.shape[1]],diagonal[:x.shape[1]]]=-cp.inf
        peak=cp.max(logterms,axis=(1,2),keepdims=True);weights=cp.exp(logterms-peak)
        weights/=cp.sum(weights,axis=(1,2),keepdims=True)
        grad=2*s*(weights/dist)@x;grad-=cp.sum(grad*x,axis=2,keepdims=True)*x
        return grad
    # Double-precision independent CPU/CUDA force comparison on a small control.
    rng=np.random.default_rng(31);small=normalize(rng.normal(size=(17,12)));s=16.
    cpu=riesz(small,s)[1]*s;gpu=cp.asnumpy(gpu_gradient(cp.asarray(small[None]),s)[0])
    error=float(np.max(abs(cpu-gpu)));assert error<1e-10,error
    x=cp.asarray(initial.astype(np.float32));moment=cp.zeros_like(x);variance=cp.zeros_like(x)
    best=initial.copy();best_peak=np.array([metrics(v)['max_cos'] for v in initial]);step=0;trace=[];start=time.monotonic()
    for exponent,lr in schedule:
        for _ in range(a.stage_steps):
            grad=gpu_gradient(x,float(exponent));step+=1
            moment=.9*moment+.1*grad;variance=.999*variance+.001*grad*grad
            update=(moment/(1-.9**step))/(cp.sqrt(variance/(1-.999**step))+1e-8)
            update-=cp.sum(update*x,axis=2,keepdims=True)*x
            x-=lr*update;x/=cp.sqrt(cp.sum(x*x,axis=2,keepdims=True))
        host=cp.asnumpy(x).astype(np.float64);values=[metrics(v) for v in host]
        peaks=np.array([v['max_cos'] for v in values]);improved=peaks<best_peak
        best[improved]=host[improved];best_peak=np.minimum(peaks,best_peak)
        r=dict(exponent=exponent,step=step,elapsed_seconds=time.monotonic()-start,
               overall_best_max_cos=float(best_peak.min()),
               group_best={mode:float(peaks[[j for j,z in enumerate(metadata) if z['mode']==mode]].min())
                           for mode in ['control','placement','rotation','archive']})
        trace.append(r);print(json.dumps(r),flush=True)
        np.save(a.out/'best.npy',best)
    np.save(a.out/'final.npy',host)
    (a.out/'gpu_summary.json').write_text(json.dumps(dict(**cfg,trace=trace,cuda_gradient_error=error,
        elapsed_seconds=time.monotonic()-start),indent=2)+'\n')
    del x,moment,variance,grad,update;cp.get_default_memory_pool().free_all_blocks()
    results=[];ctx=multiprocessing.get_context('spawn')
    with ProcessPoolExecutor(max_workers=a.workers,mp_context=ctx) as pool:
        for f in as_completed([pool.submit(polish,(i,str((a.out/'best.npy').resolve()),str(a.out.resolve()))) for i in range(a.batch)]):
            r=f.result();r['mode']=metadata[r['index']]['mode'];results.append(r)
            print(json.dumps({k:r[k] for k in ['index','mode','best','minimax']}),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
