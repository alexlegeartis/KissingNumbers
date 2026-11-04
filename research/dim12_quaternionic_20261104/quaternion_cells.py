"""Search whole 128-point quaternionic cells before breaking their symmetry.

This deliberately finite family uses the six frames C_r,D_r, the four xor
relations c=a xor b xor delta, and all three placements of squared weights
(1/2,1/4,1/4). Pair compatibility is computed before any numerical optimization.
The test uses floating arithmetic with separated discrete scalar products; final
claims would still require an exact certificate. Greedy failure is not an upper bound.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools
import json
from pathlib import Path
import time
import numpy as np
from numba import njit
from search import DATA, ROOT, normalize, metrics


def qmul(a,b):
    return np.r_[a[0]*b[0]-a[1:]@b[1:],a[0]*b[1:]+b[0]*a[1:]+np.cross(a[1:],b[1:])]


def frames():
    omega=np.ones(4)/2
    tau=np.array([1,1,0,0])/np.sqrt(2)
    powers=[np.array([1.,0,0,0]),omega,qmul(omega,omega)]
    c=np.array([[qmul(p,e) for e in np.eye(4)] for p in powers])
    d=np.array([[qmul(tau,q) for q in f] for f in c])
    return np.concatenate([c,d])


@njit(cache=True)
def compatibility(cells,dots):
    n=len(cells)
    result=np.zeros((n,n),dtype=np.bool_)
    weights=np.ones((3,3))*.5
    for i in range(3): weights[i,i]=np.sqrt(.5)
    for i in range(n):
        h,fi,gi,ki,di=cells[i]
        for j in range(i):
            hh,fj,gj,kj,dj=cells[j]
            good=True
            for a in range(4):
                if not good: break
                for b in range(4):
                    if not good: break
                    c=a^b^di
                    for aa in range(4):
                        if not good: break
                        for bb in range(4):
                            cc=aa^bb^dj
                            maximum=(weights[h,0]*weights[hh,0]*dots[fi,fj,a,aa]
                                     +weights[h,1]*weights[hh,1]*dots[gi,gj,b,bb]
                                     +weights[h,2]*weights[hh,2]*dots[ki,kj,c,cc])
                            if maximum>.5000000001:
                                good=False
                                break
            if good: result[i,j]=True;result[j,i]=True
    return result


def materialize(cells,chosen,fs):
    core=[]
    axial=np.concatenate([fs[0],-fs[0],fs[1],-fs[1],fs[2],-fs[2]])
    for j in range(3):
        x=np.zeros((24,12));x[:,4*j:4*j+4]=axial;core.append(x)
    for index in chosen:
        h,f,g,k,delta=cells[index]
        w=np.ones(3)*.5;w[h]=np.sqrt(.5)
        for a,b in itertools.product(range(4),repeat=2):
            c=a^b^delta
            for signs in itertools.product([-1,1],repeat=3):
                core.append(np.concatenate([signs[0]*w[0]*fs[f,a],signs[1]*w[1]*fs[g,b],
                                             signs[2]*w[2]*fs[k,c]])[None,:])
    return np.vstack(core)


def main():
    start=time.monotonic()
    out=ROOT/'runs'/'quaternion_cells';out.mkdir(parents=True,exist_ok=True)
    fs=frames()
    cells=[]
    for h in range(3):
        for f,g,k in itertools.product(range(6),repeat=3):
            if [f,g,k][h]<3: continue
            for delta in range(4): cells.append([h,f,g,k,delta])
    cells=np.array(cells,dtype=np.int64)
    dots=abs(np.einsum('aik,bjk->abij',fs,fs))
    adj=compatibility(cells,dots)
    print('candidates',len(cells),'compatible_pairs',int(adj.sum()//2),flush=True)
    rng=np.random.default_rng(20261104)
    found={}
    best=[]
    for trial in range(30000):
        chosen=[]
        available=np.ones(len(cells),dtype=bool)
        while available.any():
            ids=np.flatnonzero(available)
            # Common-neighbor degree plus randomness samples different clique basins.
            if trial%3==0:
                vertex=int(rng.choice(ids))
            else:
                degree=adj[np.ix_(ids,ids)].sum(axis=1)
                score=degree+rng.gumbel(size=len(ids))*[1.,4.,12.][trial%3]
                vertex=int(ids[np.argmax(score)])
            chosen.append(vertex)
            available &= adj[vertex]
        if len(chosen)>len(best):
            best=chosen.copy();print('new_best_cells',len(best),'trial',trial,flush=True)
        if len(chosen)>=6:
            key=tuple(sorted(chosen))
            if key not in found:
                x=materialize(cells,chosen,fs)
                m=metrics(x)
                assert m['max_cos']<.500000001
                eigen=np.linalg.eigvalsh(x.T@x)
                heavy=np.bincount(cells[chosen,0],minlength=3).tolist()
                found[key]=dict(cells=list(key),heavy_counts=heavy,frame_eigenvalues=eigen.tolist(),**m)
                # Retain up to 256 seeds; cells plus this script reproduce all others.
                if len(found)<=256:
                    np.save(out/f'code_{len(found)-1:04d}.npy',x)
            if len(chosen)>=7:
                np.save(out/'BREAKTHROUGH_CANDIDATE.npy',x)
                break
        if trial%1000==0:
            print('trial',trial,'six_cell_codes',len(found),flush=True)
    report=dict(candidate_count=len(cells),compatible_pairs=int(adj.sum()//2),
                trials=trial+1,max_cells_found=len(best),
                status='finite numerical family search; no impossibility claim',
                cells=cells.tolist(),codes=list(found.values()),elapsed_seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('done','max_cells',len(best),'codes',len(found),'seconds',report['elapsed_seconds'],flush=True)


if __name__=='__main__': main()
