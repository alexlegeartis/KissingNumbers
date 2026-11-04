#!/usr/bin/env python3
"""Polish the near-feasible 31561-point candidate with quadratic and higher-power penalties."""
import argparse,json,time
import numpy as np
import torch
from baselines import DATA,baseline


def evaluate(t,B,power,scale):
    with torch.no_grad():
        norms=torch.linalg.vector_norm(t,dim=1);U=t/norms[:,None]
        grad=torch.empty_like(U);loss=torch.zeros((),dtype=t.dtype,device=t.device)
        maximum=torch.full((),-1.,dtype=t.dtype,device=t.device)
        for a in range(0,len(U),1024):
            g=torch.zeros_like(U[a:a+1024])
            for V,weight in ((U,1.),(B,2.)):
                H=U[a:a+1024]@V.T
                if weight==1.:
                    ii=torch.arange(len(H),device=t.device);H[ii,a+ii]=-1.
                maximum=torch.maximum(maximum,H.max());H.sub_(.5-1e-10).clamp_(min=0.)
                if power==2:
                    loss+=weight*.5*torch.sum(H*H);factor=H
                else:
                    loss+=weight*(.05*torch.sum(H*H)+.9/power/scale**(power-2)*torch.sum(H**power))
                    factor=.1*H+.9/scale**(power-2)*H**(power-1)
                g+=2*(factor@V)
            grad[a:a+1024]=g
        grad-=(grad*U).sum(1)[:,None]*U;grad/=norms[:,None]
        return loss,maximum,grad,U


def run(seconds):
    torch.set_num_threads(1);assert torch.cuda.is_available() and 'A4000' in torch.cuda.get_device_name(0)
    device=torch.device('cuda:0');fixed=27720
    # Independent finite differences for the higher-power objective.
    rng=np.random.default_rng(621);small=rng.normal(size=(17,7));base=rng.normal(size=(23,7));base/=np.linalg.norm(base,axis=1)[:,None]
    loss,_,gradient,_=evaluate(torch.tensor(small,dtype=torch.float64,device=device),torch.tensor(base,dtype=torch.float64,device=device),6,.2)
    def reference(V):
        U=V/np.linalg.norm(V,axis=1)[:,None];C=U@U.T;np.fill_diagonal(C,-1)
        H=np.maximum(C-(.5-1e-10),0.);D=np.maximum(U@base.T-(.5-1e-10),0.)
        return .05*np.sum(H*H)+.9/6/.2**4*np.sum(H**6)+.1*np.sum(D*D)+1.8/6/.2**4*np.sum(D**6)
    assert abs(float(loss)-reference(small))<1e-10
    error=0.;grad=gradient.cpu().numpy()
    for _ in range(30):
        i=int(rng.integers(17));j=int(rng.integers(7));A=small.copy();C=small.copy();A[i,j]+=1e-5;C[i,j]-=1e-5
        error=max(error,abs((reference(A)-reference(C))/2e-5-grad[i,j]))
    assert error<1e-7,error
    warm=np.load(DATA/'focused21_candidate.npy');assert warm.shape==(31561,21)
    B=torch.tensor(warm[:fixed],dtype=torch.float64,device=device)
    t=torch.tensor(warm[fixed:],dtype=torch.float64,device=device,requires_grad=True)
    started=time.monotonic();last=started;bestmax=1.;best=None;evals=0;stages=[]
    class Limit(Exception):pass
    for power in (2,6,2):
        stage_start=time.monotonic();stage_limit=seconds/3;scale=.005
        opt=torch.optim.LBFGS([t],lr=1.,max_iter=2000,max_eval=3000,history_size=40,
                              tolerance_grad=1e-12,tolerance_change=1e-16,line_search_fn='strong_wolfe')
        def closure():
            nonlocal last,bestmax,best,evals
            if time.monotonic()-stage_start>=stage_limit:raise Limit()
            l,m,g,U=evaluate(t,B,power,scale);t.grad=g;evals+=1;value=float(m)
            if value<bestmax:bestmax=value;best=U.cpu().numpy().copy()
            now=time.monotonic()
            if now-last>=30:
                print(json.dumps(dict(event='polish_progress',power=power,max_variable_dot=bestmax,
                                      evaluations=evals,seconds=now-started)),flush=True);last=now
            return l
        try:opt.step(closure)
        except Limit:pass
        stages.append(dict(power=power,seconds=time.monotonic()-stage_start,max_dot=bestmax))
        # Carry the least violated iterate into the next penalty stage.
        if best is not None:t=torch.tensor(best,dtype=torch.float64,device=device,requires_grad=True)
        if bestmax<=.5+1e-11:break
    if best is not None:np.save(DATA/'polish21_candidate.npy',np.concatenate([warm[:fixed],best]))
    report=dict(dimension=21,count=31561,seconds=time.monotonic()-started,evaluations=evals,
                max_variable_dot=bestmax,max_dot=max(.5,bestmax),stages=stages,gradient_error=error,certified=False,
                scope='Numerical exploration only; an exact subsequent check is required.')
    (DATA/'polish21_result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=240);a=p.parse_args();run(a.seconds)
