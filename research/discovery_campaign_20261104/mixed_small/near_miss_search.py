"""Competitive bounded search from the published 842-point numerical near-miss.
All 842 vectors move; cluster reseeds and homotopies are not append-to-841 tests.
"""
from pathlib import Path
import json,time,resource
import numpy as np
from scipy.optimize import minimize,check_grad
from search import full_screen,BudgetExpired
P=Path(__file__).resolve().parent

def normalized(z,shape):
    Z=z.reshape(shape);norms=np.linalg.norm(Z,axis=1);return Z/norms[:,None],norms

def objective(z,shape,target=.5,beta=None):
    X,norms=normalized(z,shape);A=X@X.T;np.fill_diagonal(A,-1)
    if beta is None:
        E=np.maximum(A-target,0);loss=float(5e3*np.sum(E*E));gradient=2e4*E@X
    else:
        maximum=float(A.max());W=np.exp(beta*(A-maximum));np.fill_diagonal(W,0)
        total=float(W.sum());loss=maximum+np.log(total)/beta;gradient=2*(W/total)@X
    gradient=(gradient-X*np.sum(gradient*X,axis=1)[:,None])/norms[:,None]
    return float(loss),gradient.ravel(),float(A.max())

def squared_violation(X):
    A=X@X.T;np.fill_diagonal(A,-1);E=np.maximum(A-.5,0)
    return float(.5*np.sum(E*E))

def seed_trial(base,mode,rng):
    X=base.copy();A=X@X.T;np.fill_diagonal(A,-1);degrees=(A>.5008).sum(axis=1);n,d=X.shape
    metadata={'mode':mode,'all_points_optimized':True}
    if mode=='direct_pressure':pass
    elif mode=='global_anisotropic_homotopy':
        Q,_=np.linalg.qr(rng.normal(size=(d,d)));scales=np.exp(rng.normal(size=d)*.012)
        X=X@(Q*scales)@Q.T;metadata['axis_scale_std']=.012
    elif mode=='low_degree_cluster':
        center=int(np.argmin(degrees));cluster=np.r_[center,np.argsort(A[center])[-7:]]
        X[cluster]+=rng.normal(size=(len(cluster),d))*.025;metadata.update(center=center,cluster=cluster.tolist(),sigma=.025)
    elif mode=='hot_contact_cluster':
        E=np.maximum(A-.5,0);center=int(np.argmax(np.sum(E*E,axis=1)));cluster=np.r_[center,np.argsort(A[center])[-15:]]
        X[cluster]+=rng.normal(size=(len(cluster),d))*.012;metadata.update(center=center,cluster=cluster.tolist(),sigma=.012)
    elif mode=='coupled_opposite_clusters':
        center=int(rng.integers(n));antipode=int(np.argmin(A[center]));cluster=np.unique(np.r_[center,antipode,np.argsort(A[center])[-5:],np.argsort(A[antipode])[-5:]])
        X[cluster]+=rng.normal(size=(len(cluster),d))*.02;metadata.update(center=center,opposite=antipode,cluster=cluster.tolist(),sigma=.02)
    elif mode=='global_symmetry_breaking':X+=rng.normal(size=X.shape)*.003;metadata['sigma']=.003
    else:raise ValueError(mode)
    X/=np.linalg.norm(X,axis=1)[:,None]
    return X,metadata

def main():
    start=time.monotonic();deadline=start+225.;rng=np.random.default_rng(20261104842)
    test=rng.normal(size=(7,3));errors=[]
    for target,beta in [(.5,None),(.5,12000.)]:
        err=check_grad(lambda z:objective(z,(7,3),target,beta)[0],lambda z:objective(z,(7,3),target,beta)[1],test.ravel())
        scale=max(1.,np.linalg.norm(objective(test.ravel(),(7,3),target,beta)[1]));relative=float(err/scale)
        if relative>2e-6:raise RuntimeError('Near-miss gradient check failed')
        errors.append({'target':target,'beta':beta,'relative_error':relative})
    initial=np.load(P/'published842_reconstructed.npz')['vectors'];shape=initial.shape
    initial_screen=full_screen(initial);best=(initial_screen['maximum_cosine'],initial.copy());trials=[];calls=0
    modes=['direct_pressure','global_anisotropic_homotopy','low_degree_cluster','hot_contact_cluster','coupled_opposite_clusters','global_symmetry_breaking']
    for number,mode in enumerate(modes):
        if time.monotonic()>deadline:break
        X,metadata=seed_trial(best[1],mode,rng);trial_start=time.monotonic();trial_deadline=min(deadline,trial_start+39.)
        before=full_screen(X);phases=[];trialbest=(before['maximum_cosine'],X.copy());aborted=False
        schedule=([(.5,None,450),(.5,15000.,350),(.5,60000.,450)] if mode=='direct_pressure' else [(.5022,None,160),(.5012,None,160),(.5005,None,220),(.5,None,320),(.5,30000.,420)])
        for target,beta,iterations in schedule:
            localbest=None
            def fun(z):
                nonlocal localbest,best,trialbest,calls
                if time.monotonic()>trial_deadline:raise BudgetExpired
                value,gradient,maximum=objective(z,shape,target,beta);calls+=1
                if localbest is None or value<localbest[0]:localbest=(value,z.copy())
                if maximum<trialbest[0]:trialbest=(maximum,normalized(z,shape)[0].copy())
                if maximum<best[0]:best=(maximum,normalized(z,shape)[0].copy())
                return value,gradient
            try:
                result=minimize(fun,X.ravel(),jac=True,method='L-BFGS-B',options={'maxiter':iterations,'maxls':30,'ftol':1e-15,'gtol':1e-11,'maxcor':16})
                X=normalized(result.x,shape)[0];status=int(result.status);nit=int(result.nit)
            except BudgetExpired:
                if localbest is not None:X=normalized(localbest[1],shape)[0]
                status='wall_budget';nit=None;aborted=True
            screen=full_screen(X);phases.append({'target':target,'smooth_beta':beta,'max_iterations':iterations,'actual_iterations':nit,'status':status,'screen':screen})
            if aborted or best[0]<=.5+1e-10:break
        screen=full_screen(trialbest[1]);filename=f'near842_trial{number}_{mode}.npz'
        np.savez_compressed(P/filename,vectors=trialbest[1])
        row={'trial':number,'mode':metadata,'file':filename,'initial_screen':before,'best_screen':screen,'phases':phases,'seconds':time.monotonic()-trial_start,'squared_violation_best':squared_violation(trialbest[1])}
        trials.append(row)
        print(json.dumps({'trial':number,'mode':mode,'before_maximum_cosine':before['maximum_cosine'],'best_maximum_cosine':screen['maximum_cosine'],'global_best_maximum_cosine':best[0],'seconds':time.monotonic()-start}),flush=True)
        np.savez_compressed(P/'near842_best.npz',vectors=best[1])
        if best[0]<=.5+1e-10:break
    best_screen=full_screen(best[1]);np.savez_compressed(P/'near842_best.npz',vectors=best[1])
    report={'random_seed':20261104842,'source':'published842_reconstructed.npz; source bytes and primary URL pinned in near_miss_source.json','all_points_free':True,'initial_screen':initial_screen,'initial_squared_violation':squared_violation(initial),'best_screen':best_screen,'best_squared_violation':squared_violation(best[1]),'improvement_in_maximum_cosine':initial_screen['maximum_cosine']-best_screen['maximum_cosine'],'trials':trials,'objective_evaluations':calls,'wall_budget_seconds':225,'seconds':time.monotonic()-start,'gradient_checks':errors,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'status':'Numerical near-miss search; not an exact or rigorous certificate; no new SOTA unless independently certified.'}
    (P/'near_miss_results.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
