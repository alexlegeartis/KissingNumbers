"""Recovery control around a known feasible 841, distinct from discovery."""
from common import *
from direct_control import scaled_relax
from concurrent.futures import ProcessPoolExecutor,as_completed
from threadpoolctl import threadpool_limits
import time

def trial(index):
    threadpool_limits(1);rng=np.random.default_rng(841505000+index);start=time.monotonic()
    x=np.load(OLD/'data'/'record841.npy');scale=[.0001,.001,.01,.03][index%4]
    x=normalize(x+rng.normal(scale=scale,size=x.shape));initial=metrics(x)
    x,info=scaled_relax(x,16000,target=.5-1e-8);m=metrics(x)
    path=ROOT/'runs'/'repair841'/f'candidate_{index:03d}.npy';np.save(path,x)
    r=dict(index=index,scale=scale,initial=initial,final=m,optimizer=info,elapsed_seconds=time.monotonic()-start)
    if m['max_cos']<=.5+1e-12:r['exact_certificate']=verify(path)
    save_json(path.with_suffix('.json'),r);return r

if __name__=='__main__':
    reserve();out=ROOT/'runs'/'repair841';out.mkdir(parents=True,exist_ok=True);results=[]
    with ProcessPoolExecutor(max_workers=16) as pool:
        for f in as_completed([pool.submit(trial,i) for i in range(16)]):
            r=f.result();results.append(r);save_json(out/'summary.json',results)
            print(json.dumps(r),flush=True)
