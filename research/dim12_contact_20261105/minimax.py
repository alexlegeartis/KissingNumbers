"""Smooth minimax continuation, with explicitly measured tangent force residuals."""
from common import *
from scipy.optimize import minimize
import argparse
import time

class SoftMax:
    def __init__(self,shape,tau):
        self.shape=shape;self.tau=tau;self.i,self.j=np.triu_indices(shape[0],1);self.calls=0
    def evaluate(self,raw,details=False):
        raw=np.asarray(raw).reshape(self.shape);norm=np.linalg.norm(raw,axis=1,keepdims=True);x=raw/norm
        g=x@x.T;dots=g[self.i,self.j];peak=dots.max();w=np.exp((dots-peak)/self.tau);z=w.sum();w/=z
        matrix=np.zeros_like(g);matrix[self.i,self.j]=w;matrix[self.j,self.i]=w
        force=matrix@x;force-=(force*x).sum(axis=1,keepdims=True)*x
        if details:return dict(i=self.i,j=self.j,weight=w,force=force,peak=peak)
        self.calls+=1
        value=(peak-.5)+self.tau*np.log(z)
        return float(value),(force/norm).ravel()
    __call__=evaluate

def run(x,out,steps=2000,taus=(1e-4,3e-5,1e-5,3e-6,1e-6,3e-7,1e-7,3e-8,1e-8)):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);x=normalize(x);best=x.copy();bm=metrics(x)
    trace=[];start=time.monotonic();np.save(out/'initial.npy',x)
    for tau in taus:
        f=SoftMax(x.shape,tau)
        r=minimize(f,x.ravel(),jac=True,method='L-BFGS-B',
                   options=dict(maxiter=steps,maxfun=steps*3,maxcor=30,ftol=1e-16,gtol=1e-11,maxls=40))
        x=normalize(r.x.reshape(x.shape));m=metrics(x);d=f.evaluate(x,True)
        if m['max_cos']<bm['max_cos']:best=x.copy();bm=m
        row=dict(tau=tau,nit=int(r.nit),calls=f.calls,message=str(r.message),**m,
                 max_tangent_force=float(abs(d['force']).max()),force_norm=float(np.linalg.norm(d['force'])),
                 elapsed_seconds=time.monotonic()-start)
        trace.append(row);np.save(out/'best.npy',best);np.save(out/'final.npy',x)
        np.savez(out/'multipliers.npz',**d);save_json(out/'trace.json',trace)
        print(json.dumps(row),flush=True)
        if bm['max_cos']<.5:break
    save_json(out/'summary.json',dict(best=bm,trace=trace,source_sha256=sha(__file__)))
    return best,trace

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=OLD/'data'/'published842.npy')
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'minimax842');p.add_argument('--steps',type=int,default=2000)
    a=p.parse_args();reserve();run(np.load(a.input),a.out,a.steps)
