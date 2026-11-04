"""Polish the independently optimized 5D 40-point control; numerical only."""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parent

def main():
    start=time.perf_counter();data=np.load(P/'d5_c5_k4_f1,2.npz');G=data['generators'];Gs=np.array([np.linalg.matrix_power(G,t) for t in range(5)]);k,d=data['seeds'].shape
    def residual(z):
        X=z.reshape(k,d);X=X/np.linalg.norm(X,axis=1)[:,None]
        Y=np.einsum('tad,jd->jta',Gs,X);D=np.einsum('ia,jta->ijt',X,Y)
        R=np.maximum(abs(D)-.5,0);R[np.arange(k),np.arange(k),0]=0
        return R.ravel()
    result=least_squares(residual,data['seeds'].ravel(),max_nfev=1000,ftol=1e-15,xtol=1e-15,gtol=1e-15)
    X=result.x.reshape(k,d);X/=np.linalg.norm(X,axis=1)[:,None]
    Y=np.einsum('tad,jd->jta',Gs,X).reshape(20,5);V=np.r_[Y,-Y];A=V@V.T;np.fill_diagonal(A,-1)
    report={'dimension':5,'points':40,'baseline':40,'maximum_cosine':float(A.max()),'loss':float(result.fun@result.fun),'function_evaluations':result.nfev,'seconds':time.perf_counter()-start,'status':'Numerical reconstruction control only; no exact certificate and no new SOTA.'}
    np.savez_compressed(P/'d5_c5_40_polished_control.npz',seeds=X,vectors=V,generator=G)
    (P/'control_refinement.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
