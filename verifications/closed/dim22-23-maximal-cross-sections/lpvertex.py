"""Exact-ish test for a free point: maximise |z| over P = {z : |<z,xhat_i>| <= 1/2}.
A free point exists iff max|z| >= 1.  Iterated LP (objective <- current vertex) from
random starts; each fixed point is a vertex of P that locally maximises |z|."""
import numpy as np, time, sys
from scipy.optimize import linprog
from base import leech, lam22, lam23

def unitize(S,d):
    X=S.astype(np.float64); X/=np.linalg.norm(X,axis=1,keepdims=True)
    Vt=np.linalg.svd(X[:3000],full_matrices=False)[2][:d]
    Y=X@Vt.T; Y/=np.linalg.norm(Y,axis=1,keepdims=True)
    return Y

def maxnorm(Y, tries=40, seed=0, verbose=True):
    N,d=Y.shape
    A=np.vstack([Y,-Y]); b=np.full(2*N,0.5)
    rng=np.random.default_rng(seed); best=0.0; bz=None
    for t in range(tries):
        c=rng.normal(size=d); c/=np.linalg.norm(c)
        for it in range(30):
            r=linprog(-c, A_ub=A, b_ub=b, bounds=[(None,None)]*d, method='highs')
            if r.status!=0: break
            z=r.x; nz=np.linalg.norm(z)
            if nz<1e-12: break
            cn=z/nz
            if np.dot(cn,c)>1-1e-12: c=cn; break
            c=cn
        if nz>best: best,bz=nz,z.copy()
        if verbose: print("   try %2d  |z|=%.9f   best %.9f"%(t,nz,best), flush=True)
    return best,bz

if __name__=="__main__":
    A=leech()
    which=sys.argv[1] if len(sys.argv)>1 else "22"
    if which=="22": S,d,name=lam22(A)[0],22,"Lambda_22"
    elif which=="23": S,d,name=lam23(A)[0],23,"Lambda_23"
    else: S,d,name=A,24,"Leech"
    Y=unitize(S,d)
    t0=time.time()
    best,z=maxnorm(Y,tries=int(sys.argv[2]) if len(sys.argv)>2 else 25)
    print("%s: max |z| over P = %.9f   (free point iff >= 1)   [%.0fs]"%(name,best,time.time()-t0))
    print("   => min-max |cos| = %.9f"%(0.5/best))
    np.save("lpz_%s.npy"%name, z)
