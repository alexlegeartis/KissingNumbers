"""Bounded augmented-Lagrangian sphere search after explicit symmetry release.
No record coordinates are used. Positive multipliers remember violated pairs.
"""
from pathlib import Path
import json,time,resource
import numpy as np
from scipy.optimize import minimize,check_grad
from search import full_screen,BudgetExpired
P=Path(__file__).resolve().parent

def objective(z,shape,multipliers,rho):
    Z=z.reshape(shape);norm=np.linalg.norm(Z,axis=1);X=Z/norm[:,None]
    g=X@X.T-.5
    W=np.maximum(multipliers+rho*g,0);np.fill_diagonal(W,0)
    loss=float(.5/rho*np.sum(W*W))
    grad=2*W@X;grad=(grad-X*np.sum(grad*X,axis=1)[:,None])/norm[:,None]
    return loss,grad.ravel()

def smooth_objective(z,shape,beta):
    Z=z.reshape(shape);norm=np.linalg.norm(Z,axis=1);X=Z/norm[:,None]
    A=X@X.T;np.fill_diagonal(A,-1)
    maximum=float(A.max());W=np.exp(beta*(A-maximum));np.fill_diagonal(W,0)
    total=float(W.sum());loss=maximum+np.log(total)/beta
    grad=2*(W/total)@X
    grad=(grad-X*np.sum(grad*X,axis=1)[:,None])/norm[:,None]
    return float(loss),grad.ravel()

def run(source,rng,seconds):
    start=time.monotonic();deadline=start+seconds
    X=np.load(P/source['file'])['vectors'];X=X+rng.normal(size=X.shape)*.001;X/=np.linalg.norm(X,axis=1)[:,None]
    shape=X.shape;lam=np.zeros((len(X),len(X)));rho=20.;best=(float('inf'),None,None);outerlog=[];calls=0
    for outer in range(60):
        if time.monotonic()>deadline:break
        localbest=None
        def fun(z):
            nonlocal localbest,calls
            if time.monotonic()>deadline:raise BudgetExpired
            value,grad=(smooth_objective(z,shape,min(500.,120.+outer*12)) if outer%3==2 else objective(z,shape,lam,rho));calls+=1
            if localbest is None or value<localbest[0]:localbest=(value,z.copy())
            return value,grad
        try:
            result=minimize(fun,X.ravel(),jac=True,method='L-BFGS-B',options={'maxiter':220,'maxls':25,'ftol':1e-12,'gtol':1e-8,'maxcor':12})
            X=result.x.reshape(shape);status=int(result.status)
        except BudgetExpired:
            if localbest is not None:X=localbest[1].reshape(shape)
            status='wall_budget'
        X/=np.linalg.norm(X,axis=1)[:,None]
        screen=full_screen(X)
        if screen['maximum_cosine']<best[0]:best=(screen['maximum_cosine'],X.copy(),screen)
        outerlog.append({'outer':outer,'rho':rho,'mode':'smooth_max' if outer%3==2 else 'augmented_lagrangian','status':status,'screen':screen,'seconds':time.monotonic()-start})
        print(json.dumps({'job':source['config']['name']+'_alm','outer':outer,'maximum_cosine':screen['maximum_cosine'],'violations':screen['violating_pairs_at_1e-10'],'seconds':time.monotonic()-start}),flush=True)
        if screen['numerically_admissible']:break
        lam=np.maximum(lam+rho*(X@X.T-.5),0);np.fill_diagonal(lam,0)
        rho=min(800.,rho*1.15)
        if outer and outer%7==0:
            # A controlled partial shake releases a different contact network.
            indices=rng.choice(len(X),max(1,len(X)//12),replace=False)
            X[indices]+=rng.normal(size=(len(indices),shape[1]))*.025
            X/=np.linalg.norm(X,axis=1)[:,None]
            lam[indices,:]*=.2;lam[:,indices]*=.2
    name=source['config']['name']+'_alm';filename=name+'.npz'
    np.savez_compressed(P/filename,vectors=best[1])
    report={'name':name,'file':filename,'source':source['file'],'method':'all point vectors free; augmented Lagrangian alternates smooth-max inner solves, positive pair multipliers and partial symmetry-breaking shakes','wall_budget_seconds':seconds,'seconds':time.monotonic()-start,'objective_evaluations':calls,'outer_iterations':outerlog,'best_screen':best[2],'status':'numerical point set only; exact_certificate false'}
    (P/(name+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    return report

def main():
    started=time.monotonic();rng=np.random.default_rng(2026110442)
    Z=rng.normal(size=(7,3));L=rng.random(size=(7,7));L=(L+L.T)/2;np.fill_diagonal(L,0)
    err=check_grad(lambda z:objective(z,(7,3),L,20)[0],lambda z:objective(z,(7,3),L,20)[1],Z.ravel())
    scale=max(1.,np.linalg.norm(objective(Z.ravel(),(7,3),L,20)[1]));relative=float(err/scale)
    if relative>2e-6:raise RuntimeError('Augmented gradient check failed')
    smooth_error=check_grad(lambda z:smooth_objective(z,(7,3),120)[0],lambda z:smooth_objective(z,(7,3),120)[1],Z.ravel())
    if smooth_error>2e-6:raise RuntimeError('Smooth maximum gradient check failed')
    releases=json.loads((P/'release_results.json').read_text())['runs']
    selected=[min((r for r in releases if r['best_screen']['dimension']==d),key=lambda r:r['best_screen']['maximum_cosine']) for d in (12,11)]
    reports=[run(r,rng,90 if r['best_screen']['dimension']==12 else 70) for r in selected]
    report={'random_seed':2026110442,'relative_gradient_error':relative,'smooth_max_gradient_error':float(smooth_error),'runs':reports,'seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'status':'Numerical only. No new SOTA without an admissible independently certified result.'}
    (P/'alm_results.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
