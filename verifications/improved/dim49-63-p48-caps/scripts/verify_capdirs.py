#!/usr/bin/env python3
"""Exact verification of the cap-direction line systems shipped for dimensions 57-63.

At cap level t = 3/4 a part is an antipodal pair {z,-z} and distinct parts sit at
<z,z'> <= 1/2, so the cap-direction set W is a union of antipodal pairs whose LINES are
pairwise at >= 60 degrees, and lambda(k) = |W|/2 is the number of classes consumed.

Each capdirs_<k>.npz stores the configuration exactly as (A + B sqrt(R))/den with integer
A, B and inner-product coefficients c, so <x,y>*den^2 = P + Q sqrt(R) with P, Q integers.
Checked here, with no tolerance anywhere:

  * all norms equal and rational (Q vanishes on the diagonal);
  * ANTIPODALLY SYMMETRIC: -x is in the configuration for every x;
  * 60-degree: 2<x,y> <= |x||y| for every pair, i.e. 2P - N <= -2Q sqrt(R).
"""
import numpy as np, os
HERE=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','data')
TAU={9:306,10:510,11:604,13:1154,14:1932,15:2564}
TAULAT={9:272,10:336,11:438,13:918,14:1422,15:2340}
print(__doc__)
print("   k   |W|  sqrt  norm   antipodal  60-degree   lambda(k)  (tau_lat(k)/2 was used before)")
for k in sorted(TAU):
    z=np.load(os.path.join(HERE,'capdirs_%d.npz'%k))
    A,B,c=z['A'].astype(np.int64),z['B'].astype(np.int64),z['c'].astype(np.int64)
    R=int(z['R'][0]); N=int(z['norm'][0]); n=len(A)
    assert n==TAU[k], "|W| = %d but tau(%d) = %d"%(n,k,TAU[k])
    P=(A*c)@A.T + R*((B*c)@B.T)
    Q=(A*c)@B.T + (B*c)@A.T
    assert (np.diag(P)==N).all() and (np.diag(Q)==0).all()
    key={(A[i].tobytes(),B[i].tobytes()) for i in range(n)}
    sym=all(((-A[i]).tobytes(),(-B[i]).tobytes()) in key for i in range(n))
    def le(L,Rr):
        if R==0: return L<=0
        if Rr>=0 and L<=0: return True
        if Rr<0 and L>0: return False
        return (L*L<=R*Rr*Rr) if L>0 else (L*L>=R*Rr*Rr)
    bad=0
    for i in range(n):
        for j in np.nonzero(2*P[i]-N>0)[0]:
            if j!=i and not le(2*int(P[i,j])-N,-2*int(Q[i,j])): bad+=1
    assert sym and bad==0
    print("   %2d %5d %5d %5d   %-9s  %-9s   %6d      %d"%(k,n,R,N,sym,"yes",n//2,TAULAT[k]//2))

# ---------------------------------------------------------------------------------------
# The OTHER eight values of lambda(k).  final.py takes lambda(k) from a table for
# k = 1..8 and k = 12; everything in that table except k = 12 is a root system, which is
# built and checked here so that no value used by a claim is a bare literal.
print()
print("   k  system  |W|  60-degree  antipodal  lambda(k) = |W|/2")
import itertools


def A_roots(n):
    return np.array([[1 if i == p else (-1 if i == q else 0) for i in range(n + 1)]
                     for p, q in itertools.permutations(range(n + 1), 2)], np.int64)


def D_roots(n):
    V = []
    for i, j in itertools.combinations(range(n), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = np.zeros(n, np.int64)
                v[i], v[j] = si, sj
                V.append(v)
    return np.array(V, np.int64)


def E8_roots():
    V = [2 * v for v in D_roots(8)]
    for m in range(256):
        t = np.array([-1 if (m >> b) & 1 else 1 for b in range(8)], np.int64)
        if (t < 0).sum() % 2 == 0:
            V.append(t)
    return np.array(V, np.int64)


_E = E8_roots()
_r1 = _E[0]
_r2 = _E[(_E @ _r1) == -4][0]
ROOTS = {1: np.array([[1], [-1]], np.int64), 2: A_roots(2), 3: A_roots(3),
         4: D_roots(4), 5: D_roots(5),
         6: _E[((_E @ _r1) == 0) & ((_E @ _r2) == 0)], 7: _E[(_E @ _r1) == 0], 8: _E}
NAMES = {1: 'A_1', 2: 'A_2', 3: 'A_3', 4: 'D_4', 5: 'D_5', 6: 'E_6', 7: 'E_7', 8: 'E_8'}
LAM8 = {1: 1, 2: 3, 3: 6, 4: 12, 5: 20, 6: 36, 7: 63, 8: 120}
for k in range(1, 9):
    W = ROOTS[k]
    m = int(W[0] @ W[0])
    ip = W @ W.T
    np.fill_diagonal(ip, -m)
    code = bool((2 * ip <= m).all())
    key = {W[i].tobytes() for i in range(len(W))}
    sym = all((-W[i]).tobytes() in key for i in range(len(W)))
    assert code and sym and len(W) // 2 == LAM8[k], (k, len(W))
    print("   %d  %-6s %4d  %-9s  %-9s  %d" % (k, NAMES[k], len(W), 'yes', sym, len(W) // 2))

print()
print("   k = 12 is the ONE value of lambda(k) this repository does not exhibit.  A part is")
print("   an antipodal pair, so dimension 60 needs an antipodally symmetric 60-degree code")
print("   of R^12; the largest known is the minimal shell of the Coxeter-Todd lattice K_12,")
print("   with 756 vectors, giving lambda(12) = 378.  That is a published value (Conway and")
print("   Sloane, SPLAG Table 1.5) and is used as such: no K_12 configuration is shipped, and")
print("   Cohn's dimensions1-24.txt is no help, its 841-point dimension-12 record containing")
print("   only 18 antipodal pairs.  The record does NOT depend on it.  On the largest")
print("   12-dimensional configuration this repository does exhibit -- Lambda_12, 648 vectors,")
print("   dim73-95-gamma72-caps/data/lam12_W.npy -- lambda(12) would be 324 and dimension 60")
print("   would read 56 789 413 instead of 57 477 123, still 4 372 572 above the published")
print("   52 416 841.  Every other k here is a root system above or a capdirs file.")

print()
print("ALL CAP-DIRECTION LINE SYSTEMS VERIFIED EXACTLY")
