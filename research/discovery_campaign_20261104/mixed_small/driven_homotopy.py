"""Final bounded coupled-force homotopy and high-beta minimax on 842 points."""
from pathlib import Path
import json,time,resource
import numpy as np
from scipy.optimize import minimize,check_grad
from near_miss_search import objective,normalized,squared_violation
from search import full_screen,BudgetExpired
P=Path(__file__).resolve().parent

def forced(z,shape,target,beta,U,strength):
    value,grad,maximum=objective(z,shape,target,beta)
    if strength:
        X,norm=normalized(z,shape)
        value-=strength*float(np.sum(U*X))
        extra=-strength*U
        extra=(extra-X*np.sum(extra*X,axis=1)[:,None])/norm[:,None]
        grad+=extra.ravel()
    return value,grad,maximum

def main():
    start=time.monotonic();deadline=start+95;rng=np.random.default_rng(20261104843)
    X0=np.load(P/'published842_reconstructed.npz')['vectors'];shape=X0.shape
    A=X0@X0.T;np.fill_diagonal(A,-1);degrees=(A>.5008).sum(axis=1)
    best=(full_screen(X0)['maximum_cosine'],X0.copy());logs=[];calls=0
    test=rng.normal(size=(5,3));Utest=rng.normal(size=(5,3))
    error=check_grad(lambda z:forced(z,(5,3),.5,None,Utest,.3)[0],lambda z:forced(z,(5,3),.5,None,Utest,.3)[1],test.ravel())
    relative=float(error/max(1.,np.linalg.norm(forced(test.ravel(),(5,3),.5,None,Utest,.3)[1])))
    if relative>2e-6:raise RuntimeError('Driven gradient check failed')
    names=['high_beta_minimax','low_degree_cluster_drive','global_coupled_drive','hot_cluster_drive']
    for trial,name in enumerate(names):
        if time.monotonic()>deadline:break
        t=time.monotonic();limit=min(deadline,t+(25 if trial==0 else 23));X=best[1].copy();U=rng.normal(size=shape)
        U-=X*np.sum(U*X,axis=1)[:,None];U/=np.maximum(np.linalg.norm(U,axis=1)[:,None],1e-12)
        if trial==1:
            cluster=np.argsort(degrees)[:18];mask=np.zeros(len(X));mask[cluster]=1;U*=mask[:,None]
        elif trial==3:
            center=int(np.argmax(np.sum(np.maximum(A-.5008,0),axis=1)));cluster=np.unique(np.r_[center,np.argsort(A[center])[-23:]])
            mask=np.zeros(len(X));mask[cluster]=1;U*=mask[:,None]
        phases=[];trialbest=(float('inf'),None)
        schedule=([(.5,1e5,0.,200),(.5,4e5,0.,220),(.5,1.6e6,0.,350)] if trial==0 else [(.5035,None,.4 if trial==1 else .12,180),(.502,None,.1,180),(.5008,None,.015,180),(.5,None,0.,240),(.5,4e5,0.,300)])
        for target,beta,strength,nit in schedule:
            last=None
            def fun(z):
                nonlocal calls,last,best,trialbest
                if time.monotonic()>limit:raise BudgetExpired
                value,grad,maximum=forced(z,shape,target,beta,U,strength);calls+=1;last=z.copy()
                if maximum<trialbest[0]:trialbest=(maximum,normalized(z,shape)[0].copy())
                if maximum<best[0]:best=(maximum,normalized(z,shape)[0].copy())
                return value,grad
            try:
                r=minimize(fun,X.ravel(),jac=True,method='L-BFGS-B',options={'maxiter':nit,'maxls':25,'ftol':1e-15,'gtol':2e-12,'maxcor':16})
                X=normalized(r.x,shape)[0];status=int(r.status)
            except BudgetExpired:
                if last is not None:X=normalized(last,shape)[0]
                status='wall_budget'
            phases.append({'target':target,'beta':beta,'force_strength':strength,'status':status,'maximum_cosine':full_screen(X)['maximum_cosine']})
            if status=='wall_budget' or best[0]<=.5+1e-10:break
        file='near842_'+name+'.npz';np.savez_compressed(P/file,vectors=trialbest[1])
        row={'name':name,'file':file,'phases':phases,'best_screen':full_screen(trialbest[1]),'seconds':time.monotonic()-t};logs.append(row)
        print(json.dumps({'mode':name,'trial_best':trialbest[0],'global_best':best[0],'seconds':time.monotonic()-start}),flush=True)
        if best[0]<=.5+1e-10:break
    np.savez_compressed(P/'near842_driven_best.npz',vectors=best[1])
    result={'random_seed':20261104843,'source':'published842_reconstructed.npz','all_points_free':True,'gradient_relative_error':relative,'source_screen':full_screen(X0),'best_screen':full_screen(best[1]),'source_squared_violation':squared_violation(X0),'best_squared_violation':squared_violation(best[1]),'runs':logs,'objective_evaluations':calls,'wall_budget_seconds':95,'seconds':time.monotonic()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'status':'Numerical only; no exact/interval certificate; no new kissing-number claim.'}
    (P/'driven_results.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
