"""Force selected high-stress contacts open, then release the added constraints."""
from common import *
from search import Hinge
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor,as_completed
from threadpoolctl import threadpool_limits
import time
import argparse

def force_open(x,i,j,target,strength,steps):
    f=Hinge(x.shape)
    def objective(raw):
        val,grad=f(raw);rr=raw.reshape(x.shape);norm=np.linalg.norm(rr,axis=1,keepdims=True);z=rr/norm
        dot=np.einsum('ij,ij->i',z[i],z[j]);gap=np.maximum(dot-target,0);val+=strength*(gap@gap)
        force=np.zeros_like(z);np.add.at(force,i,2*strength*gap[:,None]*z[j]);np.add.at(force,j,2*strength*gap[:,None]*z[i])
        force-=(force*z).sum(axis=1,keepdims=True)*z
        return float(val),grad+(force/norm).ravel()
    opt=minimize(objective,x.ravel(),jac=True,method='L-BFGS-B',
                 options=dict(maxiter=steps,maxfun=steps*3,maxcor=20,ftol=1e-16,gtol=1e-10,maxls=40))
    return normalize(opt.x.reshape(x.shape)),dict(nit=int(opt.nit),message=str(opt.message))

def trial(spec):
    index,steps,out=spec;threadpool_limits(1);out=Path(out);start=time.monotonic();rng=np.random.default_rng(20261105842900+index)
    data=np.load(ROOT/'runs'/'stress842'/'stress_1.npz');x=data['coordinates'].copy();weight=data['weight'];rank=np.argsort(weight)[::-1]
    count=[1,2,4,8][index%4];gap=[.02,.05,.10,.20][(index//4)%4]
    mode=(index//16)%3
    if mode==0:selected=rank[(index//48)%16:(index//48)%16+count]
    elif mode==1:selected=rng.choice(rank[:128],size=count,replace=False,p=weight[rank[:128]]/weight[rank[:128]].sum())
    else:
        vertex=int(rng.choice(np.argsort(data['score'])[-24:]));incident=rank[(data['i'][rank]==vertex)|(data['j'][rank]==vertex)]
        selected=incident[:count]
    i=data['i'][selected];j=data['j'][selected];trace=[]
    for strength in [.5,5.]:
        x,info=force_open(x,i,j,.5-gap,strength,steps)
        trace.append(dict(stage='forced',strength=strength,**metrics(x),forced_max=float(np.max(np.einsum('ij,ij->i',x[i],x[j]))),**info))
    forced=x.copy();x,info=relax(x,steps*2);energy=x.copy();em=metrics(x);best=x.copy();bm=em
    for power in [4,16,64]:
        x,info=optimize(x,800,power=power);m=metrics(x);trace.append(dict(stage='released',power=power,**m))
        if m['max_cos']<bm['max_cos']:best=x.copy();bm=m
    path=out/f'candidate_{index:03d}.npy';np.save(path,best);np.save(out/f'energy_{index:03d}.npy',energy)
    np.save(out/f'forced_{index:03d}.npy',forced)
    r=dict(index=index,edges=list(zip(i.tolist(),j.tolist())),requested_gap=gap,selection_mode=mode,
           energy=em,best=bm,trace=trace,elapsed_seconds=time.monotonic()-start)
    if bm['max_cos']<=.5+1e-12:r['exact_certificate']=verify(path)
    save_json(out/f'candidate_{index:03d}.json',r);return r

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=144);p.add_argument('--workers',type=int,default=24)
    p.add_argument('--steps',type=int,default=2500);p.add_argument('--out',type=Path,default=ROOT/'runs'/'break_contacts')
    a=p.parse_args();reserve();a.out.mkdir(parents=True,exist_ok=True);results=[]
    save_json(a.out/'config.json',dict(trials=a.trials,workers=a.workers,steps=a.steps,source_sha256=sha(__file__)))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(trial,(i,a.steps,str(a.out.resolve()))) for i in range(a.trials)]):
            r=f.result();results.append(r);save_json(a.out/'summary.json',results)
            print(json.dumps({k:r[k] for k in ['index','requested_gap','best','elapsed_seconds']}),flush=True)
