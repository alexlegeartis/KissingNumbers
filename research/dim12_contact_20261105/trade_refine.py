"""Continuously release low-conflict individual-ray trade candidates."""
from common import *
from concurrent.futures import ProcessPoolExecutor,as_completed
from threadpoolctl import threadpool_limits
import time
from direct_control import scaled_relax

def trial(spec):
    index,pool,path,count,out=spec;threadpool_limits(1);out=Path(out);start=time.monotonic()
    data=np.load(pool/'rays.npz');tokens=path.read_text().split();ids=np.array([int(z) for z in tokens[1:]])
    x=data['coordinates'][ids];rng=np.random.default_rng(8418427000+index)
    if count==841:
        g=np.maximum(x@x.T-.5,0);np.fill_diagonal(g,0);scores=(g*g).sum(axis=1)
        x=np.delete(x,int(np.argmax(scores)),axis=0)
    noise=[0.,.002,.01,.03][index%4];x=normalize(x+rng.normal(scale=noise,size=x.shape));initial=metrics(x)
    x,info=relax(x,6000);energy=x.copy();em=metrics(x);best=x.copy();bm=em
    if count==841 and em['max_cos']<.501:
        x,info=scaled_relax(x,16000,target=.5-1e-8);m=metrics(x)
        if m['max_cos']<bm['max_cos']:best=x.copy();bm=m
    else:
        for power in [4,16,64]:
            x,info=optimize(x,1000,power=power);m=metrics(x)
            if m['max_cos']<bm['max_cos']:best=x.copy();bm=m
    dest=out/f'candidate_{index:03d}.npy';np.save(dest,best);np.save(out/f'energy_{index:03d}.npy',energy)
    r=dict(index=index,pool=str(pool),source=str(path),count=count,noise=noise,
           discrete_conflicts=int(tokens[0]),initial=initial,energy=em,best=bm,elapsed_seconds=time.monotonic()-start)
    if bm['max_cos']<=.5+1e-12:r['exact_certificate']=verify(dest)
    save_json(dest.with_suffix('.json'),r);return r

if __name__=='__main__':
    reserve();out=ROOT/'runs'/'trade_refine';out.mkdir(parents=True,exist_ok=True);specs=[]
    for poolname,runname in [('trade_pool','trades842'),('trade_pool_extended','trades_extended842')]:
        files=sorted((ROOT/'runs'/runname).glob('thread_*_cost_6.txt'))
        for path in files:
            for count in [841,842]:specs.append((len(specs),ROOT/'runs'/poolname,path,count,out))
    save_json(out/'config.json',dict(trials=len(specs),workers=16,source_sha256=sha(__file__)))
    results=[]
    with ProcessPoolExecutor(max_workers=16) as pool:
        for f in as_completed([pool.submit(trial,s) for s in specs]):
            r=f.result();results.append(r);save_json(out/'summary.json',results)
            print(json.dumps({k:r[k] for k in ['index','count','best','elapsed_seconds']}),flush=True)
