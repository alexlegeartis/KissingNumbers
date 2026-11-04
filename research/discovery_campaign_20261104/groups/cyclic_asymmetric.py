"""Nonantipodal cyclic orbit unions with independently optimized seeds.
The optional single fixed pole breaks antipodal symmetry in odd dimensions.
"""
from pathlib import Path
import json,time,resource
import numpy as np
from scipy.optimize import minimize,check_grad
from cyclic_continuous import matrices
P=Path(__file__).resolve().parent

def objective(flat,k,d,G,pole):
    Z=flat.reshape(k,d);norm=np.linalg.norm(Z,axis=1);X=Z/norm[:,None]
    Y=np.einsum('tad,jd->jta',G,X);dots=np.einsum('ia,jta->ijt',X,Y)
    err=np.maximum(dots-.5,0);err[np.arange(k),np.arange(k),0]=0
    C=2*err
    grad=np.einsum('ijt,jta->ia',C,Y)+np.einsum('ijt,tad,ia->jd',C,G,X)
    loss=float(np.sum(err**2))
    if pole:
        e=np.maximum(X[:,0]-.5,0);loss+=len(G)*float(e@e);grad[:,0]+=2*len(G)*e
    grad=(grad-X*np.sum(grad*X,axis=1)[:,None])/norm[:,None]
    return loss,grad.ravel()

def run(d,q,k,freq,pole,restarts,seconds,rng):
    G=matrices(d,q,freq);start=time.perf_counter();best=None;trials=[]
    for trial in range(restarts):
        if trial and time.perf_counter()-start>seconds:break
        Z=rng.normal(size=(k,d))
        if best is not None and trial%4:
            Z=best[1]+rng.normal(size=(k,d))*[0,.02,.1,.4][trial%4]
        res=minimize(objective,Z.ravel(),args=(k,d,G,pole),jac=True,method='L-BFGS-B',options={'maxiter':3000,'ftol':1e-15,'gtol':3e-11,'maxls':35})
        X=res.x.reshape(k,d);X/=np.linalg.norm(X,axis=1)[:,None]
        V=np.einsum('tad,jd->jta',G,X).reshape(k*q,d)
        if pole:V=np.concatenate((V,np.eye(d)[:1]))
        A=V@V.T;np.fill_diagonal(A,-1);maximum=float(A.max())
        trials.append({'trial':trial,'loss':float(res.fun),'max_cosine_full':maximum,'iterations':int(res.nit),'status':int(res.status)})
        if best is None or res.fun<best[0]:best=(float(res.fun),X.copy(),maximum,V.copy())
        if maximum<=.5+1e-9:break
    tag=f'asym_d{d}_c{q}_k{k}_pole{int(pole)}_f'+','.join(map(str,freq))
    np.savez_compressed(P/(tag+'.npz'),seeds=best[1],vectors=best[3],generator=G[1])
    result={'dimension':d,'group_order':q,'seed_orbits':k,'fixed_pole':pole,'target_cardinality':k*q+int(pole),'frequencies':freq,'restart_limit':restarts,'restart_wall_budget_seconds':seconds,'wall_budget_checked_between_restarts':True,'optimizer_max_iterations':3000,'trials_completed':len(trials),'seconds':time.perf_counter()-start,'best_loss':best[0],'best_max_cosine':best[2],'admissible_numerically':best[2]<=.5+1e-9,'file':tag+'.npz','trials':trials,'certificate':'Numerical optimization only. Failed arrays are not kissing configurations.'}
    print(json.dumps({a:b for a,b in result.items() if a!='trials'}),flush=True)
    return result

def main():
    start=time.perf_counter();rng=np.random.default_rng(20261104)
    G=matrices(5,5,[1,2]);v=rng.normal(size=40)
    err=check_grad(lambda z:objective(z,8,5,G,True)[0],lambda z:objective(z,8,5,G,True)[1],v)
    assert err<1e-5
    jobs=[(5,5,8,[1,2],True,160,55),(5,7,6,[1,2],False,160,55),(6,5,15,[1,1,2],False,160,55),(6,7,11,[1,2,3],False,160,55),(7,7,18,[1,2,3],True,120,70)]
    results=[run(*j,rng) for j in jobs]
    report={'outcome':'No SOTA unless a fully screened admissible result explicitly says so.','random_seed':20261104,'gradient_finite_difference_error':float(err),'baselines':{'5':40,'6':72,'7':126},'results':results,'wall_seconds':time.perf_counter()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'Fully independent seed vectors with no antipodal constraint; reducible real cyclic representation with repeated irreducibles allowed; optional one-sided fixed pole. Bounded multi-start search, no nonexistence claim.'}
    (P/'asymmetric_results.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
