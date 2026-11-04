"""Resume independent members of a batched Adam control on CPU workers.

The mean loss gradient is divided by the original batch size. Raw coordinates,
Adam moments, step counter, and the remaining source schedule are preserved.
Changing BLAS/backend changes roundoff, so this is not a bitwise GPU replay.
"""
from common import *
from control import SCHEDULE
import argparse
import importlib.util
import multiprocessing as mp
import time
from concurrent.futures import ProcessPoolExecutor,as_completed

def trial(index):
    start=time.monotonic();x=CHECKPOINT['x'][index:index+1].clone().requires_grad_(True)
    opt=torch.optim.Adam([x],lr=SCHEDULE[FIRST][2]);old=CHECKPOINT['optimizer']['state'][0]
    opt.state[x]={k:(v[index:index+1].clone() if v.ndim else v.clone()) for k,v in old.items()}
    trace=[];best=np.inf;bestx=None
    for stage in range(FIRST,len(SCHEDULE)):
        s,steps,lr=SCHEDULE[stage]
        for group in opt.param_groups:group['lr']=lr
        for k in range(steps):
            loss,_,_=SOURCE.riesz_log_loss_and_stats(x,s);loss=loss/BATCH
            opt.zero_grad(set_to_none=True);loss.backward();opt.step()
        with torch.no_grad():z=(x/x.norm(dim=2,keepdim=True))[0].numpy().copy()
        m=metrics(z)
        if m['max_cos']<best:best=m['max_cos'];bestx=z.copy()
        row=dict(stage=stage,s=s,steps=steps,**m,elapsed_seconds=time.monotonic()-start);trace.append(row)
        np.save(OUT/f'candidate_{index:03d}.npy',bestx);np.save(OUT/f'final_{index:03d}.npy',z)
        save_json(OUT/f'candidate_{index:03d}.json',dict(index=index,best_max=best,trace=trace))
    return dict(index=index,best_max=best,trace=trace,elapsed_seconds=time.monotonic()-start)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--workers',type=int,default=64)
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'control841_cpu');a=p.parse_args();reserve();OUT=a.out;OUT.mkdir(parents=True,exist_ok=True)
    import torch
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    CHECKPOINT=torch.load(a.checkpoint,map_location='cpu',weights_only=True)
    FIRST=int(CHECKPOINT['completed_stage'])+1;BATCH=len(CHECKPOINT['x'])
    spec=importlib.util.spec_from_file_location('source_search',AUDIT/'inputs'/'search_841_riesz.py')
    SOURCE=importlib.util.module_from_spec(spec);spec.loader.exec_module(SOURCE)
    save_json(OUT/'config.json',dict(checkpoint=str(a.checkpoint),checkpoint_sha256=sha(a.checkpoint),
                                    batch=BATCH,first_stage=FIRST,workers=a.workers,torch_version=torch.__version__,
                                    loss_divisor=BATCH,source_sha256=sha(__file__),schedule=SCHEDULE[FIRST:]))
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers,mp_context=mp.get_context('fork')) as pool:
        for future in as_completed([pool.submit(trial,i) for i in range(BATCH)]):
            r=future.result();results.append(r);save_json(OUT/'summary.json',results)
            print(json.dumps(dict(index=r['index'],best_max=r['best_max'],elapsed_seconds=r['elapsed_seconds'])),flush=True)
