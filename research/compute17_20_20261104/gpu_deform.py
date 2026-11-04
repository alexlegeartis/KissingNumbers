#!/usr/bin/env python3
"""Bounded multistart deformation; float64 GPU gradients in memory-bounded chunks."""
import argparse,json,time
import numpy as np
import torch
from continuous import baseline
from orbits17 import DATA


def energy_gradient(t,threshold):
    with torch.no_grad():
        norms=torch.linalg.vector_norm(t,dim=1)
        U=t/norms[:,None]
        grad=torch.empty_like(U);energy=torch.zeros((),dtype=torch.float64,device=t.device)
        worst=torch.full((),-1.,dtype=torch.float64,device=t.device)
        block=1024 if U.shape[1]==20 else 2048
        for a in range(0,len(U),block):
            H=U[a:a+block]@U.T
            b=len(H);ii=torch.arange(b,device=t.device)
            H[ii,a+ii]=-1.
            worst=torch.maximum(worst,H.max())
            H.sub_(threshold).clamp_(min=0.)
            energy+=.5*torch.sum(H*H)
            grad[a:a+block]=2.*(H@U)
        grad-=(grad*U).sum(1)[:,None]*U
        grad/=norms[:,None]
        return energy,worst,grad,U


def check_gradient(device):
    rng=np.random.default_rng(171720);X=rng.normal(size=(37,17))
    t=torch.tensor(X,dtype=torch.float64,device=device)
    e,_,g,_=energy_gradient(t,.2);gradient=g.cpu().numpy()
    def loss(A):
        U=A/np.linalg.norm(A,axis=1)[:,None];G=U@U.T
        np.fill_diagonal(G,-1.);D=np.maximum(G-.2,0.)
        return .5*np.sum(D*D)
    assert abs(float(e)-loss(X))<1e-10
    worst=0.
    for _ in range(30):
        i,j=map(int,[rng.integers(len(X)),rng.integers(X.shape[1])]);eps=1e-5
        A=X.copy();A[i,j]+=eps;B=X.copy();B[i,j]-=eps
        numerical=(loss(A)-loss(B))/(2*eps)
        worst=max(worst,abs(numerical-gradient[i,j]))
    assert worst<1e-7,worst
    return worst


def run(dimension,seconds,seed,device):
    torch.set_num_threads(1)
    assert torch.cuda.is_available()
    info=torch.cuda.get_device_properties(device)
    assert 'A4000' in info.name, info.name
    cuda=torch.device(f'cuda:{device}')
    gradient_error=check_gradient(cuda)
    X=baseline(dimension)
    holes=np.load(DATA/f'holes{dimension}.npy')
    rng=np.random.default_rng(seed);started=time.monotonic()
    bestenergy=float('inf');bestmax=None;best=None;evals=0;restarts=0;last=started
    threshold=.5-1e-10
    print(json.dumps(dict(event='gpu_start',dimension=dimension,device=info.name,
                          dtype='float64',seconds=seconds,hole_directions=len(holes),
                          gradient_finite_difference_max_error=gradient_error)),flush=True)
    class Limit(Exception):pass
    while time.monotonic()-started<seconds:
        h=holes[restarts%len(holes)]
        V=np.concatenate([X,h[None,:]])
        perturb=[0.,.0001,.002,.01,.03,.06][restarts%6]
        V+=rng.normal(size=V.shape)*perturb
        t=torch.tensor(V,dtype=torch.float64,device=cuda,requires_grad=True)
        opt=torch.optim.LBFGS([t],lr=1.,max_iter=80,max_eval=120,tolerance_grad=1e-10,
                              tolerance_change=1e-14,history_size=20,line_search_fn='strong_wolfe')
        local_evals=0
        def closure():
            nonlocal bestenergy,bestmax,best,evals,local_evals,last
            if time.monotonic()-started>=seconds:raise Limit()
            with torch.no_grad():
                energy,worst,grad,U=energy_gradient(t,threshold)
                t.grad=grad
                value=float(energy);maximum=float(worst)
                evals+=1;local_evals+=1
                if value<bestenergy:
                    bestenergy=value;bestmax=maximum;best=U.cpu().numpy().copy()
                now=time.monotonic()
                if now-last>=30:
                    print(json.dumps(dict(event='gpu_progress',dimension=dimension,restarts=restarts,
                                          evaluations=evals,energy=bestenergy,max_dot=bestmax,
                                          seconds=now-started)),flush=True)
                    last=now
                return energy
        try:opt.step(closure)
        except Limit:break
        restarts+=1
        if bestmax is not None and bestmax<=.5+1e-11:break
    report=dict(dimension=dimension,count=len(X)+1,device=info.name,seconds=time.monotonic()-started,
                seed=seed,restarts=restarts,evaluations=evals,best_energy=bestenergy,max_dot=bestmax,
                gradient_finite_difference_max_error=gradient_error,
                certified=False,scope='Numerical exploration only; exact certification is required for a record.')
    if best is not None:np.save(DATA/f'gpu{dimension}_candidate.npy',best)
    (DATA/f'gpu{dimension}_result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dimension',type=int,choices=(17,20),required=True)
    p.add_argument('--seconds',type=float,default=180);p.add_argument('--seed',type=int,default=1)
    p.add_argument('--device',type=int,default=0);a=p.parse_args()
    run(a.dimension,a.seconds,a.seed,a.device)
