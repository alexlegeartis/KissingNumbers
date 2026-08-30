"""Delsarte LP bounds needed for the LAYER lever in dims 22/23."""
import numpy as np
from scipy.optimize import linprog
def gegen(n, K, ts):
    G=np.zeros((K+1,len(ts))); G[0]=1.0
    if K>=1: G[1]=ts
    for k in range(1,K):
        G[k+1]=((2*k+n-2)*ts*G[k]-k*G[k-1])/(k+n-2)
    return G
def lp_bound(n, s, K=40, M=6000):
    ts=np.linspace(-1.0, s, M); G=gegen(n,K,ts)
    r=linprog(np.ones(K), A_ub=G[1:].T, b_ub=-np.ones(M), bounds=[(0,None)]*K, method='highs')
    return None if r.status!=0 else 1.0+r.fun

if __name__ == "__main__":
    print("reference kissing LP bounds:")
    for n in (21,22,23,24):
        print("   A(%d, 1/2) <= %.1f"%(n, lp_bound(n,0.5)))
    
    print("\nlayer lever: a dim-(n+1) code split over an axis at cosine-height h.")
    print("  equator  : code in R^n, cos<=1/2")
    print("  layer +-h: code in R^n, cos <= (1/2-h^2)/(1-h^2)  [internal]")
    print("  total    <= A(n,1/2) + 2 A(n, s(h)) + ... \n")
    for n,tot,name in ((21,49896,"dim22"),(22,93150,"dim23")):
        print("  --- %s (record %d), equator dim %d, A(%d,1/2) <= %.0f ---"%(name,tot,n,n,lp_bound(n,0.5)))
        for h2,lab in ((3/16,"h=sqrt3/4 (Lambda frame)"), (1/4,"h=1/2"), (1/8,"h=1/(2sqrt2)"),
                       (1/16,"h=1/4"), (3/8,"h=sqrt(3/8)"), (1/3,"h=1/sqrt3")):
            s=(0.5-h2)/(1-h2)
            b=lp_bound(n,s)
            print("     %-24s internal cos<=%.5f   A(%d,.)<=%9.1f   equator+2*layer <= %.0f"
                  %(lab,s,n,b, lp_bound(n,0.5)+2*b))
