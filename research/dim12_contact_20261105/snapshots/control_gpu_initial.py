"""Run the original PyTorch search functions with their complete default schedule.

Only checkpointing/logging and a selectable base construction are added. The
optimizer, raw-coordinate parametrization, loss, and default rates are unchanged.
The separate 750,000-step published polishing schedule is not run by this file.
"""
from common import *
import argparse
import importlib.util
import time

SCHEDULE=[(8,1000,.005),(16,1000,.003),(32,1000,.002),(64,2000,.001),
          (128,2000,.0005),(256,2000,.0002),(512,2000,.0001),
          (1024,4000,.00005),(2048,4000,.00001),(4096,4000,.00001),
          (10000,4000,.000005),(20000,4000,.000001),(40000,4000,.000001)]

def main():
    p=argparse.ArgumentParser();p.add_argument('--batch',type=int,default=64)
    p.add_argument('--seed',type=int,default=2026110501);p.add_argument('--alternative',action='store_true')
    p.add_argument('--benchmark',type=int,default=0);p.add_argument('--out',type=Path,default=ROOT/'runs'/'control841')
    a=p.parse_args();reserve();a.out.mkdir(parents=True,exist_ok=True)
    import torch
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    assert torch.cuda.is_available()
    torch.cuda.set_device(0);device=torch.cuda.get_device_name(0)
    assert device=='NVIDIA RTX A4000',device
    spec=importlib.util.spec_from_file_location('authors_search',AUDIT/'inputs'/'search_841_riesz.py')
    source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
    torch.manual_seed(a.seed);np.random.seed(a.seed%(2**32))
    base=source.build_clebsch_840('cuda',torch.float64)
    if a.alternative:
        sys.path.insert(0,str(AUDIT));from alternative_bridges import coefficients
        aa,bb=coefficients();base=torch.as_tensor((aa+np.sqrt(2)*bb)/12,device='cuda',dtype=torch.float64)
        base=base/base.norm(dim=1,keepdim=True)
    x=source.init_batch(base,a.batch,'cuda',torch.float64)
    opt=torch.optim.Adam([x],lr=SCHEDULE[0][2])
    config=dict(batch=a.batch,seed=a.seed,alternative=a.alternative,schedule=SCHEDULE,
                torch_version=torch.__version__,device=device,dtype='float64',
                original_source_sha256=sha(AUDIT/'inputs'/'search_841_riesz.py'),
                wrapper_sha256=sha(__file__),cpu_affinity=sorted(os.sched_getaffinity(0)))
    save_json(a.out/'config.json',config)
    np.save(a.out/'initial.npy',x.detach().cpu().numpy())
    total=0;trace=[];started=time.monotonic();best=np.full(a.batch,np.inf);best_x=None
    for stage,(s,steps,lr) in enumerate(SCHEDULE):
        for group in opt.param_groups:group['lr']=lr
        for k in range(steps):
            loss,maxdot,_=source.riesz_log_loss_and_stats(x,s)
            opt.zero_grad(set_to_none=True);loss.backward();opt.step();total+=1
            if total%250==0 or (a.benchmark and total==a.benchmark):
                torch.cuda.synchronize()
                row=dict(step=total,stage=stage,s=s,lr=lr,elapsed_seconds=time.monotonic()-started,
                         min_pre_step=float(maxdot.min()),mean_pre_step=float(maxdot.mean()),
                         max_pre_step=float(maxdot.max()))
                trace.append(row);print(json.dumps(row),flush=True);save_json(a.out/'trace.json',trace)
            if a.benchmark and total>=a.benchmark:return
        with torch.no_grad():
            z=x/x.norm(dim=2,keepdim=True);host=z.cpu().numpy()
        values=np.array([metrics(v)['max_cos'] for v in host])
        if best_x is None:best_x=host.copy()
        improved=values<best;best_x[improved]=host[improved];best=np.minimum(best,values)
        np.save(a.out/'best.npy',best_x);np.save(a.out/'final.npy',host)
        torch.save(dict(x=x.detach(),optimizer=opt.state_dict(),completed_stage=stage),a.out/'checkpoint.pt')
        save_json(a.out/'stage_summary.json',dict(completed_stage=stage,best=best.tolist(),final=values.tolist(),
                                                elapsed_seconds=time.monotonic()-started))
    save_json(a.out/'summary.json',dict(**config,steps=total,best=best.tolist(),final=values.tolist(),
                                     elapsed_seconds=time.monotonic()-started))

if __name__=='__main__':main()
