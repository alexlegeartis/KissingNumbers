"""Bounded continuous seed optimization in explicit real cyclic representations.
Numerical candidates are not certificates. No known record coordinates are read.
"""
from pathlib import Path
import itertools,json,time,resource,os
import numpy as np
from scipy.optimize import minimize,check_grad
P=Path(__file__).resolve().parent

def matrices(d,q,freq):
    ans=[]
    for t in range(q):
        G=np.eye(d);offset=d%2
        for j,f in enumerate(freq):
            a=2*np.pi*f*t/q;c=np.cos(a);s=np.sin(a);k=offset+2*j
            G[k:k+2,k:k+2]=[[c,-s],[s,c]]
        ans.append(G)
    return np.array(ans)

def objective(flat,k,d,G):
    Z=flat.reshape(k,d);norm=np.linalg.norm(Z,axis=1);X=Z/norm[:,None]
    Y=np.einsum('tad,jd->jta',G,X)
    dots=np.einsum('ia,jta->ijt',X,Y)
    error=np.maximum(abs(dots)-.5,0)
    error[np.arange(k),np.arange(k),0]=0
    C=2*error*np.sign(dots)
    grad=np.einsum('ijt,jta->ia',C,Y)+np.einsum('ijt,tad,ia->jd',C,G,X)
    grad=(grad-X*np.sum(grad*X,axis=1)[:,None])/norm[:,None]
    return float(np.sum(error**2)),grad.ravel()

def run(d,q,k,freq,restarts,seconds,rng):
    G=matrices(d,q,freq);start=time.perf_counter();best=None;logs=[]
    for trial in range(restarts):
        if trial and time.perf_counter()-start>seconds:break
        Z=rng.normal(size=(k,d))
        if best is not None and trial%3:Z=best[1]+rng.normal(size=(k,d))*(.03 if trial%3==1 else .2)
        res=minimize(objective,Z.ravel(),args=(k,d,G),jac=True,method='L-BFGS-B',options={'maxiter':1800,'ftol':1e-14,'gtol':1e-9,'maxls':30})
        X=res.x.reshape(k,d);X/=np.linalg.norm(X,axis=1)[:,None]
        Y=np.einsum('tad,jd->jta',G,X).reshape(k*q,d);V=np.concatenate((Y,-Y))
        A=V@V.T;np.fill_diagonal(A,-1);maximum=float(A.max())
        logs.append({'trial':trial,'loss':float(res.fun),'max_cosine_full':maximum,'iterations':int(res.nit),'status':int(res.status)})
        if best is None or res.fun<best[0]:best=(float(res.fun),X.copy(),maximum,V.copy())
        if maximum<=.5+1e-10:break
    tag=f'd{d}_c{q}_k{k}_f'+','.join(map(str,freq))
    np.savez_compressed(P/(tag+'.npz'),seeds=best[1],vectors=best[3],generators=G[1])
    result={'dimension':d,'group_order':q,'seed_orbits':k,'target_cardinality':2*k*q,'frequencies':freq,'restart_limit':restarts,'restart_wall_budget_seconds':seconds,'wall_budget_checked_between_restarts':True,'optimizer_max_iterations':1800,'trials_completed':len(logs),'seconds':time.perf_counter()-start,'best_loss':best[0],'best_max_cosine':best[2],'admissible_numerically':best[2]<=.5+1e-10,'file':tag+'.npz','trials':logs,'certificate':'Numerical optimization only. Inadmissible arrays are failed trial point sets, not kissing configurations.'}
    print(json.dumps({a:b for a,b in result.items() if a!='trials'}),flush=True)
    return result

def main():
    start=time.perf_counter();rng=np.random.default_rng(20261104)
    G=matrices(5,7,[1,2]);v=rng.normal(size=15)
    err=check_grad(lambda z:objective(z,3,5,G)[0],lambda z:objective(z,3,5,G)[1],v)
    assert err<1e-5
    jobs=[(5,7,3,[1,2],24,45),(5,5,4,[1,2],16,35),(6,19,2,[1,4,7],24,40),(6,19,2,[1,2,3],24,40),(7,8,8,[1,2,3],12,60),(8,11,11,[1,2,3,4],10,80)]
    results=[run(*j,rng) for j in jobs]
    report={'outcome':'No SOTA unless a fully screened admissible result explicitly says so. Failed point sets are saved for replay, not certificates.','random_seed':20261104,'gradient_finite_difference_error':float(err),'baselines':{'5':40,'6':72,'7':126,'8':240},'results':results,'wall_seconds':time.perf_counter()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'scope':'Continuous asymmetric seed vectors, independent amplitudes and phases in reducible cyclic representations; antipodality is retained. Bounded multi-start search, no nonexistence claim.'}
    (P/'continuous_results.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
