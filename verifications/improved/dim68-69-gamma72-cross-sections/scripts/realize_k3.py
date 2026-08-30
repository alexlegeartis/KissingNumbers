#!/usr/bin/env python3
"""Dimension 69 recovered: the three-point Gram (4,-2,0) of Gamma_72 is REALIZABLE.

The k-point LP certifies, for a FIXED Gram G of minimal vectors u_1..u_k, a lower bound on
n[0..0] = #{v minimal : <u_i,v> = 0 for all i}; those v are pairwise at >= 60 degrees and
span a (72-k)-dimensional space, so tau(72-k) >= n[0..0].  The bound is valid for any ACTUAL
k-tuple with that Gram.

dim69-71/derive.py could only certify realizability THROUGH the LP -- a witnessing
(k-1)-point cell with strictly positive minimum -- and for every three-point Gram except the
pairwise-orthogonal one that cell has minimum 0.  Dimension 69 was therefore given up at
149 415 783, below the floor 331 737 984, and 527 597 644 was withdrawn as unusable.

Gamma_72 is available explicitly, so realizability does not need the LP: exhibit the triple.
Inner products of minimal vectors are x G y^T with G the catalogue Gram.

    python realize_k3.py
"""
import numpy as np, os
HERE=os.path.dirname(os.path.abspath(__file__))
D=os.path.join(HERE,'..','..','dim73-95-gamma72-caps','data')
if not os.path.exists(os.path.join(D, 'gamma72_pool.npy')):
    raise SystemExit(
        "This is the SEARCH that found a triple for each three-point Gram, not a verifier.  "
        "It needs a pool of Gamma_72 minimal vectors (488 MB), which this repository does not "
        "ship.  What it found IS shipped, and scripts/verify69.py checks it from scratch "
        "using nothing but the Gram.")
G=np.load(os.path.join(D,'gamma72_gram.npy')).astype(np.int64)
P=np.load(os.path.join(D,'gamma72_pool.npy'),mmap_mode='r')
print(__doc__)
X=np.asarray(P[:400000],dtype=np.int64)
nrm=np.einsum('ij,jk,ik->i',X,G,X)
assert set(nrm.tolist())=={8}
print("pool sample: %d vectors, all of norm 8 under the catalogue Gram"%len(X),flush=True)
XG=X@G
LP={(4,-2,0):527597644,(4,-1,0):500027708,(4,0,0):469123586,(3,0,0):409062220,
    (2,0,0):346628061,(1,0,0):245233175,(0,0,0):149415783}
rng=np.random.default_rng(0)
print()
print("  Gram (a,b,c)   realizable   LP certificate   beats the floor 331737984?")
best=None
for (a,b,c),cert in LP.items():
    trip=None
    for _ in range(400):
        i=int(rng.integers(len(X)))
        u1=X[i]; ip1=XG@u1
        J=np.nonzero(ip1==a)[0]; J=J[J!=i]
        if len(J)==0: continue
        for j in J[:40]:
            u2=X[j]; ip2=XG@u2
            K=np.nonzero((ip1==b)&(ip2==c))[0]
            K=K[(K!=i)&(K!=j)]
            if len(K):
                trip=np.array([u1,u2,X[int(K[0])]]); break
        if trip is not None: break
    if trip is None:
        print("   (%2d,%2d,%2d)      no          %-14d"%(a,b,c,cert)); continue
    M=trip@G@trip.T
    ok=bool((np.diag(M)==8).all() and M[0,1]==a and M[0,2]==b and M[1,2]==c and M[1,0]==a)
    print("   (%2d,%2d,%2d)      %-11s %-14d  %s"%(a,b,c,"YES" if ok else "MISMATCH",cert,
          "YES" if cert>331737984 else "no"))
    if ok and (best is None or cert>best[0]): best=(cert,trip,(a,b,c),M)
assert best is not None
cert,trip,g,M=best
np.save(os.path.join(HERE,'..','data','gamma72_triple_k3.npy'),trip)
print()
print("BEST REALIZABLE THREE-POINT GRAM: %s, exact Gram"%(str(g)))
for r in M.tolist(): print("     ",r)
print("   determinant %d"%int(round(np.linalg.det(M.astype(float)))))
print("   tau(69) >= %d      floor 331737984, factor %.2f"%(cert,cert/331737984))
print("   saved the witnessing triple to gamma72_triple_k3.npy")
