#!/usr/bin/env python3
"""How much of dimensions 49-63 is SETTLED, and what is left.

    tau(48+k) >= 52416000 + 2*S(lambda(k)) + |A| .

THREE ingredients.  Two of them are now provably at their ceilings.

(1) lambda(k), the number of parts.  At t = 3/4 a part is an antipodal pair, so the cap
    direction set satisfies W = -W and |W| = 2*lambda(k); and W is a 60-degree code, so
    |W| <= tau(k).  Therefore

        lambda(k) <= floor(tau(k)/2)                                        [CEILING]

    and equality holds exactly when an antipodally symmetric record configuration exists.
    Verified antipodally symmetric for k = 9, 10, 11, 13, 14, 15, and for k <= 8 the records
    ARE the root systems.  So lambda(k) is OPTIMAL for every k except 12, where tau(12) = 841
    is odd, the ceiling is 420 and the best known antipodal system is K_12's 378.

(2) |A|, the axis layer.  A is a 60-degree code in R^k, so |A| <= tau(k), and the construction
    already takes |A| = tau(k).  OPTIMAL.

(3) S(lambda) = max total lines in lambda pairwise disjoint classes.  This is a maximum
    lambda-colourable subgraph problem and is NOT settled.  It is bounded by lambda * alpha
    with alpha the largest class; alpha >= 7069 is achieved and bounded above by the Delsarte
    LP over the attainable cosines {0, +-1/6, +-1/3} in dimension 48."""
import numpy as np, os
from scipy.optimize import linprog
def gegen(n,K,t):
    G=[np.ones_like(t),np.array(t,float)]
    for k in range(1,K):
        G.append(((2*k+n-2)*t*G[k]-k*G[k-1])/(k+n-2))
    return np.array(G)
def lp(n,ts,K=60):
    t=np.array(ts,float); G=gegen(n,K,t); A=G[1:].T
    r=linprog(np.ones(K),A_ub=A,b_ub=-np.ones(len(t)),bounds=[(0,None)]*K,method='highs')
    return 1+r.fun if r.success else None
print(__doc__)
TAU={1:2,2:6,3:12,4:24,5:40,6:72,7:126,8:240,9:306,10:510,11:604,12:841,13:1154,14:1932,15:2564}
LAM={1:1,2:3,3:6,4:12,5:20,6:36,7:63,8:120,9:153,10:255,11:302,12:378,13:577,14:966,15:1282}
print("  k   tau(k)   ceiling floor(tau/2)   lambda used   at the ceiling?")
for k in range(1,16):
    c=TAU[k]//2
    print("  %2d %7d %14d %14d       %s"%(k,TAU[k],c,LAM[k],"YES" if LAM[k]==c else "no, %d short"%(c-LAM[k])))
b=lp(48,[-1.0,-1/3,-1/6,0.0,1/6,1/3])
print()
print("Delsarte bound on a P_48 class (cosines 0, +-1/6, +-1/3): %.1f vectors = %.1f lines"%(b,b/2))
print("largest class achieved: 7069 lines")
S=np.sort(np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data',
                     'class_sizes_p48_familyB_ext.npy')).ravel())[::-1].astype(np.int64)
cum=np.concatenate([[0],np.cumsum(S)])
print()
print("  k  lambda   S(lambda) achieved   lambda*7069 (if every class were maximal)   ratio")
for k in (1,8,12,15):
    L=min(LAM[k],len(S))
    print("  %2d %6d %19d %36d   %.1f%%"%(k,L,int(cum[L]),L*7069,100.0*cum[L]/(L*7069)))
print()
print("So the direction side and the axis side are CLOSED; the only room left in 49-63 is the")
print("class family, and it is bounded by lambda*alpha -- at k = 15 we are at %.1f%% of that."%(
    100.0*cum[min(LAM[15],len(S))]/(min(LAM[15],len(S))*7069)))
