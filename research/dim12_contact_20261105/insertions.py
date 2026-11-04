"""841 -> 842 searches with all blockers included in a growing free region."""
from common import *
import argparse
import time
from concurrent.futures import ProcessPoolExecutor,as_completed
from threadpoolctl import threadpool_limits
from scipy.optimize import minimize

def hole(x,rng,sparse=False):
    z=normalize(rng.normal(size=(1,12)))[0]
    def fun(raw):
        norm=np.linalg.norm(raw);z=raw/norm;excess=np.maximum(x@z-.5,0)
        grad=2*excess@x;grad-=(grad@z)*z
        return float(excess@excess),grad/norm
    out=minimize(fun,z,jac=True,method='L-BFGS-B',options=dict(maxiter=300,ftol=1e-15,gtol=1e-10))
    z=out.x/np.linalg.norm(out.x)
    if sparse:
        # Reweighted least squares favors a small set of blockers. The cap
        # term prevents the trivial solution of duplicating an existing point.
        for eps in [.05,.02,.01,.005,.002]:
            old=np.maximum(x@z-.5,0);weight=eps*eps/(old*old+eps*eps)
            def sparse_fun(raw):
                norm=np.linalg.norm(raw);v=raw/norm;dots=x@v
                excess=np.maximum(dots-.5,0);cap=np.maximum(dots-.75,0)
                grad=2*(weight*excess+4*cap)@x;grad-=(grad@v)*v
                return float((weight*excess)@excess+4*(cap@cap)),grad/norm
            out=minimize(sparse_fun,z,jac=True,method='L-BFGS-B',options=dict(maxiter=150,ftol=1e-15,gtol=1e-10))
            z=out.x/np.linalg.norm(out.x)
    dots=x@z
    return z,dict(blockers=np.flatnonzero(dots>.5).tolist(),max_dot=float(dots.max()),energy=float(out.fun))

def trial(spec):
    index,steps,out,sparse,diverse=spec;threadpool_limits(1);out=Path(out);start=time.monotonic()
    rng=np.random.default_rng(8418422026110500+index);basepath=OLD/'data'/'record841.npy'
    if diverse:
        rows=json.loads((ROOT/'runs'/'repair841'/'summary.json').read_text())
        passed=sorted(z['index'] for z in rows if z.get('exact_certificate',{}).get('passed'))
        basepath=ROOT/'runs'/'repair841'/f'candidate_{passed[index%len(passed)]:03d}.npy'
    base=normalize(np.load(basepath));candidates=[hole(base,rng,sparse) for _ in range(8)]
    z,hm=min(candidates,key=lambda row:(len(row[1]['blockers']),row[1]['energy']))
    k=[16,32,64,128][index%4];blockers=np.array(hm['blockers'],dtype=int)
    w=np.maximum(base@base.T-.499,0);np.fill_diagonal(w,0)
    priority=np.zeros(841);priority[blockers]=1/max(len(blockers),1)
    restart=priority.copy();degree=np.maximum(w.sum(axis=1),1e-30)
    for _ in range(8):priority=.4*restart+.6*w@(priority/degree)
    priority[blockers]=1e6
    ids=np.argsort(priority)[-max(k,len(blockers)):];fixed=np.setdiff1d(np.arange(841),ids)
    free=np.vstack([base[ids],z]);mode=(index//4)%3
    if mode:
        free[:-1]=normalize(free[:-1]+rng.normal(scale=[0,.025,.075][mode],size=free[:-1].shape))
    initial=metrics(np.vstack([base,z]));free,info=relax(free,steps,fixed=base[fixed])
    local=metrics(np.vstack([base[fixed],free]));x=base.copy();x[ids]=free[:-1];x=np.vstack([x,free[-1]])
    # Release a larger boundary before the global solve, preserving which point
    # is the inserted one and the original 841-point row correspondence.
    larger=np.r_[np.argsort(priority)[-min(2*k,841):],841]
    fixed2=np.setdiff1d(np.arange(842),larger);y,info=relax(x[larger],steps,fixed=x[fixed2]);x[larger]=y
    x,info=relax(x,steps*2);energy=x.copy();em=metrics(x);best=x.copy();bm=em;trace=[]
    for p in [4,16,64]:
        x,info=optimize(x,600,power=p);m=metrics(x);trace.append(dict(power=p,**m))
        if m['max_cos']<bm['max_cos']:best=x.copy();bm=m
    path=out/f'insert_{index:04d}.npy';np.save(path,best);np.save(out/f'energy_{index:04d}.npy',energy)
    r=dict(index=index,patch_size=k,mutation=mode,hole=hm,sparse=sparse,base_path=str(basepath),selected=ids.tolist(),initial=initial,local=local,
           energy=em,best=bm,trace=trace,elapsed_seconds=time.monotonic()-start)
    if bm['max_cos']<=.5+1e-12:r['exact_certificate']=verify(path)
    save_json(out/f'insert_{index:04d}.json',r);return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=192);p.add_argument('--workers',type=int,default=32)
    p.add_argument('--steps',type=int,default=2000);p.add_argument('--out',type=Path,default=ROOT/'runs'/'insertions842')
    p.add_argument('--sparse',action='store_true');p.add_argument('--diverse',action='store_true')
    a=p.parse_args();reserve();a.out.mkdir(parents=True,exist_ok=True);results=[]
    save_json(a.out/'config.json',dict(trials=a.trials,workers=a.workers,steps=a.steps,sparse=a.sparse,diverse=a.diverse,source_sha256=sha(__file__)))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(trial,(i,a.steps,str(a.out.resolve()),a.sparse,a.diverse)) for i in range(a.trials)]):
            r=future.result();results.append(r);save_json(a.out/'summary.json',results)
            print(json.dumps({k:r[k] for k in ['index','patch_size','best','elapsed_seconds']}),flush=True)
