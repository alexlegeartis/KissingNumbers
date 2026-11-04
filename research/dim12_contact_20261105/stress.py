"""Fit nonnegative, approximately balanced weights on near-maximal pairs.

Residuals and complementarity gaps are reported explicitly. This numerical
diagnostic is not an exact stationary-point or rigidity certificate.
"""
from common import *
from scipy.sparse import coo_matrix
from scipy.optimize import minimize
import argparse
import time

def fit(x,band,steps=20000):
    x=normalize(x);n,d=x.shape;ii,jj=np.triu_indices(n,1);dots=(x@x.T)[ii,jj];peak=dots.max()
    selected=dots>peak-band;i,j,g=ii[selected],jj[selected],dots[selected];m=len(i)
    row=np.column_stack([i[:,None]*d+np.arange(d),j[:,None]*d+np.arange(d)]).ravel()
    col=np.repeat(np.arange(m),2*d)
    value=np.column_stack([x[j]-g[:,None]*x[i],x[i]-g[:,None]*x[j]]).ravel()
    a=coo_matrix((value,(row,col)),shape=(n*d,m)).tocsr()
    # y is m times the multiplier, to keep its typical entries near one.
    def fun(y):
        force=a@y;error=y.mean()-1.
        f=(force@force)/m + m*error*error
        grad=2*(a.T@force)/m+2*error
        return float(f),grad
    start=time.monotonic()
    opt=minimize(fun,np.ones(m),jac=True,method='L-BFGS-B',bounds=[(0,None)]*m,
                 options=dict(maxiter=steps,maxfun=steps*2,maxcor=20,ftol=1e-15,gtol=1e-10))
    lam=opt.x/opt.x.sum();force=(a@lam).reshape(n,d)
    score=np.bincount(np.r_[i,j],weights=np.r_[lam,lam],minlength=n)
    info=dict(band=band,edges=m,positive_weights=int((lam>1e-8).sum()),
              max_tangent_force=float(abs(force).max()),force_norm=float(np.linalg.norm(force)),
              weighted_complementarity_gap=float(lam@(peak-g)),
              max_support_gap=float((peak-g)[lam>1e-8].max()),
              score_quantiles=np.quantile(score,[0,.1,.25,.5,.75,.9,1]).tolist(),
              top_points=np.argsort(score)[-32:][::-1].tolist(),nit=int(opt.nit),
              message=str(opt.message),elapsed_seconds=time.monotonic()-start)
    return dict(i=i,j=j,weight=lam,force=force,score=score,coordinates=x),info

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=ROOT/'runs'/'minimax842'/'best.npy')
    p.add_argument('--out',type=Path,default=ROOT/'runs'/'stress842');a=p.parse_args();reserve();a.out.mkdir(parents=True,exist_ok=True)
    x=np.load(a.input);results=[]
    for k,band in enumerate([1e-5,3e-5,1e-4]):
        data,info=fit(x,band);np.savez(a.out/f'stress_{k}.npz',**data);results.append(info)
        print(json.dumps(info),flush=True);save_json(a.out/'summary.json',results)
