"""(a) Is the configuration maximal?  Compute  m(S) = min_{|z|=1} max_i |<z, xhat_i>|.
A free point exists iff m(S) <= 1/2.   Smooth (log-sum-exp) minimisation, many restarts."""
import numpy as np, sys, time
from base import leech, lam23, lam22

def unitize(S, d):
    X = S.astype(np.float64); X /= np.linalg.norm(X,axis=1,keepdims=True)
    U,s,Vt = np.linalg.svd(X[:3000], full_matrices=False)
    B = Vt[:d].T
    Y = X @ B
    Y /= np.linalg.norm(Y,axis=1,keepdims=True)
    return Y

def minmax(Y, restarts=200, seed=0, iters=400, verbose=False):
    N,d = Y.shape
    rng = np.random.default_rng(seed)
    best = 1e9; bestz=None
    for r in range(restarts):
        z = rng.normal(size=d); z/=np.linalg.norm(z)
        beta = 40.0
        for it in range(iters):
            t = np.abs(Y@z)
            mx = t.max()
            w = np.exp(beta*(t-mx))
            w /= w.sum()
            g = (w*np.sign(Y@z)) @ Y          # gradient of soft-max of |<z,y>|
            g -= (g@z)*z
            ng = np.linalg.norm(g)
            if ng < 1e-14: break
            step = 0.5/ (1+it*0.05)
            z = z - step*g/ng
            z /= np.linalg.norm(z)
            if it%80==79: beta*=2.0
        m = np.abs(Y@z).max()
        if m < best: best, bestz = m, z.copy()
        if verbose and r%25==0: print("   restart %3d  best %.6f"%(r,best))
    return best, bestz

if __name__ == "__main__":
    A = leech()
    for name, S, d, R in (("Lambda_22", lam22(A)[0], 22, 200),
                          ("Lambda_23", lam23(A)[0], 23, 120),
                          ("Leech",     A,           24, 60)):
        Y = unitize(S,d)
        t0=time.time()
        m, z = minmax(Y, restarts=R, seed=7)
        print("%-10s N=%-6d d=%2d   min-max |cos| = %.9f   (target <= 0.5)   [%.0fs]"
              %(name, len(S), d, m, time.time()-t0))
        np.save("z_%s.npy"%name, z)
