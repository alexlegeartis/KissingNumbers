#!/usr/bin/env python3
"""Deform only the auxiliary layer and one inserted point, keeping the exact base fixed."""
import argparse,json,time
import numpy as np
import torch
from baselines import DATA,baseline


def evaluate(t,B,threshold):
    with torch.no_grad():
        norms=torch.linalg.vector_norm(t,dim=1);U=t/norms[:,None]
        grad=torch.empty_like(U);loss=torch.zeros((),dtype=t.dtype,device=t.device)
        maximum=torch.full((),-1.,dtype=t.dtype,device=t.device)
        for a in range(0,len(U),1024):
            C=U[a:a+1024]@U.T;ii=torch.arange(len(C),device=t.device);C[ii,a+ii]=-1.
            maximum=torch.maximum(maximum,C.max());C.sub_(threshold).clamp_(min=0.)
            loss+=.5*torch.sum(C*C);g=2*(C@U)
            D=U[a:a+1024]@B.T
            maximum=torch.maximum(maximum,D.max());D.sub_(threshold).clamp_(min=0.)
            loss+=torch.sum(D*D);g+=2*(D@B)
            grad[a:a+1024]=g
        grad-=(grad*U).sum(1)[:,None]*U;grad/=norms[:,None]
        return loss,maximum,grad,U


def self_test(device):
    rng=np.random.default_rng(1921);X=rng.normal(size=(23,19));B=rng.normal(size=(29,19));B/=np.linalg.norm(B,axis=1)[:,None]
    l,_,g,_=evaluate(torch.tensor(X,dtype=torch.float64,device=device),torch.tensor(B,dtype=torch.float64,device=device),.2)
    grad=g.cpu().numpy()
    def energy(Y):
        U=Y/np.linalg.norm(Y,axis=1)[:,None];C=U@U.T;np.fill_diagonal(C,-1)
        return .5*np.maximum(C-.2,0).ravel()@np.maximum(C-.2,0).ravel()+np.sum(np.maximum(U@B.T-.2,0)**2)
    assert abs(float(l)-energy(X))<1e-10
    error=0.
    for _ in range(30):
        i=int(rng.integers(len(X)));j=int(rng.integers(19));A=X.copy();C=X.copy();A[i,j]+=1e-5;C[i,j]-=1e-5
        error=max(error,abs((energy(A)-energy(C))/2e-5-grad[i,j]))
    assert error<1e-7,error
    return error


def run(dimension,seconds,seed):
    torch.set_num_threads(1);assert torch.cuda.is_available()
    assert 'A4000' in torch.cuda.get_device_name(0)
    device=torch.device('cuda:0');error=self_test(device)
    X=baseline(dimension);fixed=10668 if dimension==19 else 27720
    base=X[:fixed];aux=X[fixed:];B=torch.tensor(base,dtype=torch.float64,device=device)
    direction_files=[DATA/f'holes{dimension}_structured.npy',DATA/f'holes{dimension}.npy']
    directions=np.concatenate([np.load(p) for p in direction_files if p.exists()])
    rng=np.random.default_rng(seed);started=time.monotonic();last=started
    best=float('inf');bestmax=None;bestU=None;restarts=0;evals=0
    print(json.dumps(dict(event='focused_start',dimension=dimension,fixed_points=fixed,
                          moving_points=len(aux)+1,seconds=seconds,gradient_error=error)),flush=True)
    class Limit(Exception):pass
    while time.monotonic()-started<seconds:
        if dimension==19 and restarts%2==0:
            extra=aux[int(rng.integers(768))].copy();extra[-3:]*=-1
        elif dimension==21 and restarts%2==0:
            z=np.load(DATA/'baseline21.npz');missing=np.flatnonzero(z['labels']<0)
            w=int(z['words'][rng.choice(missing)]);p=z['profiles'][rng.integers(21)]
            extra=p*np.asarray([1-2*((w>>i)&1) for i in range(21)])/np.sqrt(164.)
        else:extra=directions[restarts%len(directions)]
        V=np.concatenate([aux,extra[None,:]])
        V+=rng.normal(size=V.shape)*[0.,0.,.0001,.001,.005,.015][restarts%6]
        t=torch.tensor(V,dtype=torch.float64,device=device,requires_grad=True)
        opt=torch.optim.LBFGS([t],lr=1.,max_iter=100,max_eval=150,history_size=20,
                              tolerance_grad=1e-10,tolerance_change=1e-14,line_search_fn='strong_wolfe')
        def closure():
            nonlocal best,bestmax,bestU,evals,last
            if time.monotonic()-started>=seconds:raise Limit()
            l,m,g,U=evaluate(t,B,.5-1e-10);t.grad=g;evals+=1
            value=float(l);maximum=float(m)
            if value<best:best=value;bestmax=maximum;bestU=U.cpu().numpy().copy()
            now=time.monotonic()
            if now-last>=30:
                print(json.dumps(dict(event='focused_progress',dimension=dimension,restarts=restarts,
                                      evaluations=evals,energy=best,max_variable_dot=bestmax,seconds=now-started)),flush=True);last=now
            return l
        try:opt.step(closure)
        except Limit:break
        restarts+=1
        if bestmax is not None and bestmax<=.5+1e-11:break
    if bestU is not None:np.save(DATA/f'focused{dimension}_candidate.npy',np.concatenate([base,bestU]))
    report=dict(dimension=dimension,count=len(X)+1,fixed_points=fixed,restarts=restarts,evaluations=evals,
                seconds=time.monotonic()-started,best_energy=best,max_variable_dot=bestmax,
                max_dot=max(.5,bestmax) if bestmax is not None else None,gradient_error=error,
                certified=False,scope='Numerical candidates only; exact checking is required for any new bound.')
    (DATA/f'focused{dimension}_result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dimension',type=int,choices=(19,21),required=True)
    p.add_argument('--seconds',type=float,default=180);p.add_argument('--seed',type=int,default=1);a=p.parse_args()
    run(a.dimension,a.seconds,a.seed)
