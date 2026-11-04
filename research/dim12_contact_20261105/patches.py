"""Rebuild neighborhoods selected from measured near-contact stresses."""
from common import *
import argparse
import time
from concurrent.futures import ProcessPoolExecutor,as_completed
from scipy.linalg import eigh,expm,null_space
from threadpoolctl import threadpool_limits

def choose_patch(x,data,k,rng,trial):
    n=len(x);score=data['score'];w=np.zeros((n,n))
    w[data['i'],data['j']]=data['weight'];w+=w.T
    seed=int(rng.choice(np.argsort(score)[-32:]))
    mode=(trial//16)%3
    if mode==0:
        walk=np.zeros(n);walk[seed]=1;restart=walk.copy()
        degree=np.maximum(w.sum(axis=1),1e-30)
        for _ in range(12):walk=.35*restart+.65*(w@(walk/degree))
        priority=walk/np.sqrt(np.maximum(degree,1e-12))
    elif mode==1:
        priority=np.abs(x@x[seed])+.03*score/max(score.max(),1e-30)
    else:
        priority=score*rng.lognormal(0,.15,size=n)
    selected=np.argsort(priority)[-k:]
    return selected,dict(seed_point=seed,selection_mode=mode,selected_stress=float(score[selected].sum()/2))

def soft_motion(x,ids,rng):
    """Small eigenmodes of the patch's tangent near-contact Jacobian."""
    k=len(ids);d=x.shape[1];q=d-1;position={int(v):i for i,v in enumerate(ids)}
    basis=np.array([null_space(x[v:v+1]).copy() for v in ids])
    h=np.zeros((k*q,k*q));g=x@x.T
    for local,global_i in enumerate(ids):
        block=slice(local*q,(local+1)*q)
        for j in np.flatnonzero(g[global_i]>.498):
            if j==global_i:continue
            a=x[j]@basis[local];h[block,block]+=np.outer(a,a)
            if int(j) in position:
                other=position[int(j)];b=x[global_i]@basis[other];op=slice(other*q,(other+1)*q)
                h[block,op]+=np.outer(a,b)
    assert np.max(abs(h-h.T))<1e-10
    ev,vectors=eigh(h,subset_by_index=(0,min(7,k*q-1)))
    v=vectors@rng.normal(size=vectors.shape[1]);v=v.reshape(k,q)
    motion=np.einsum('kij,kj->ki',basis,v);motion/=max(np.sqrt(np.mean(np.sum(motion*motion,axis=1))),1e-20)
    return motion,ev.tolist()

def trial(spec):
    index,steps,stress_path,out=spec;threadpool_limits(1);out=Path(out);start=time.monotonic()
    rng=np.random.default_rng(20261105842000+index);data=dict(np.load(stress_path));x=data['coordinates'].copy()
    k=[16,32,64,128][index%4];ids,meta=choose_patch(x,data,k,rng,index)
    mode=(index//4)%4;scale=[.04,.10,.22,.40][(index//48)%4];eigenvalues=None
    if mode==0:
        motion,eigenvalues=soft_motion(x,ids,rng);x[ids]=normalize(x[ids]+scale*motion)
    elif mode==1:
        skew=rng.normal(size=(12,12));skew=skew-skew.T;skew/=np.linalg.norm(skew)/np.sqrt(12)
        x[ids]=x[ids]@expm(scale*skew)
    elif mode==2:
        x[ids]=normalize(x[ids]+rng.normal(scale=scale/np.sqrt(12),size=(k,12)))
    else:
        # Replace a highly stressed point by splitting another selected point.
        victim=int(ids[np.argmax(data['score'][ids])]);donor=int(rng.choice(ids[ids!=victim]))
        z=x[donor].copy();v=rng.normal(size=12);v-=(v@z)*z;v/=np.linalg.norm(v)
        angle=.3+scale;x[victim]=np.cos(angle)*z+np.sin(angle)*v;x[donor]=np.cos(angle)*z-np.sin(angle)*v
    initial=metrics(x);fixed=np.setdiff1d(np.arange(len(x)),ids)
    local,local_opt=relax(x[ids],min(steps,1000),fixed=x[fixed]);x[ids]=local
    patch_metrics=metrics(x);x,info=relax(x,steps);energy=x.copy();em=metrics(x)
    best=x.copy();bm=em;trace=[]
    for p in [4,16,64]:
        x,info=optimize(x,600,power=p);m=metrics(x);trace.append(dict(power=p,**m))
        if m['max_cos']<bm['max_cos']:best=x.copy();bm=m
    path=out/f'patch_{index:04d}.npy';np.save(path,best);np.save(out/f'energy_{index:04d}.npy',energy)
    result=dict(index=index,patch_size=k,mutation_mode=mode,scale=scale,selected=ids.tolist(),**meta,
                small_eigenvalues=eigenvalues,initial=initial,patch=patch_metrics,energy=em,best=bm,
                trace=trace,elapsed_seconds=time.monotonic()-start)
    if bm['max_cos']<=.5+1e-12:result['exact_certificate']=verify(path)
    save_json(out/f'patch_{index:04d}.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=192);p.add_argument('--workers',type=int,default=32)
    p.add_argument('--steps',type=int,default=4000);p.add_argument('--stress',type=Path,default=ROOT/'runs'/'stress842'/'stress_1.npz')
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'patches842');a=p.parse_args();reserve();a.out.mkdir(parents=True,exist_ok=True)
    save_json(a.out/'config.json',dict(trials=a.trials,workers=a.workers,steps=a.steps,stress=str(a.stress),
                                      stress_sha256=sha(a.stress),source_sha256=sha(__file__)))
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(trial,(i,a.steps,str(a.stress.resolve()),str(a.out.resolve()))) for i in range(a.trials)]):
            r=future.result();results.append(r);save_json(a.out/'summary.json',results)
            print(json.dumps({k:r[k] for k in ['index','patch_size','mutation_mode','best','elapsed_seconds']}),flush=True)
