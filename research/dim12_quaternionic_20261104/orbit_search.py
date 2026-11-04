"""Explore a cyclic order-three quotient, then remove one point and release it.

Right multiplication by -omega on all three quaternion factors preserves the
canonical 840 code. Add an entire orbit (843 points), optimize its 281 orbit
generators, and then remove one high-stress point to obtain a new 842-point seed.
The symmetry is fully released for the final 12D search.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS']:
    os.environ[key]='1'
import argparse,hashlib,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from quaternion_cells import qmul
from kernels import NeighborList
from search import ROOT,DATA,normalize,metrics,relax
from continuation import optimize


def group():
    r=np.zeros((12,12))
    block=np.array([qmul(e,-np.ones(4)/2) for e in np.eye(4)])
    for j in range(3):r[j*4:j*4+4,j*4:j*4+4]=block
    assert np.max(abs(np.linalg.matrix_power(r,3)-np.eye(12)))<1e-14
    return [np.eye(12),r,r@r]


def generators(x,rotations):
    used=set();out=[]
    for i in range(len(x)):
        if i in used:continue
        orbit=[]
        for r in rotations:
            dots=x@(x[i]@r)
            j=int(np.argmax(dots));assert dots[j]>1-1e-12
            orbit.append(j)
        assert len(set(orbit))==3
        assert not any(j in used for j in orbit)
        used.update(orbit);out.append(x[i])
    assert len(used)==len(x)
    return np.array(out)


class OrbitHinge:
    def __init__(self,shape,rotations):
        self.shape=shape;self.rotations=rotations;self.neighbors=NeighborList(.5)
    def __call__(self,raw):
        raw=raw.reshape(self.shape);norms=np.linalg.norm(raw,axis=1,keepdims=True);z=raw/norms
        x=np.vstack([z@r for r in self.rotations])
        value,g=self.neighbors.evaluate(x)
        n=len(z);grad=sum(g[k*n:(k+1)*n]@r.T for k,r in enumerate(self.rotations))
        grad-=np.sum(grad*z,axis=1,keepdims=True)*z
        return value,(grad/norms).ravel()


def trial(spec):
    index,steps,rounds,out=spec;threadpool_limits(1)
    start=time.monotonic();out=Path(out);seed=84320261104+index;rng=np.random.default_rng(seed)
    rotations=group();base=generators(np.load(DATA/'canonical840.npy'),rotations)
    extra=np.zeros(12);extra[8:]=rng.normal(size=4);extra/=np.linalg.norm(extra)
    z=normalize(np.vstack([base,extra]))
    z+=rng.normal(scale=[.001,.01,.03,.08][index%4],size=z.shape);z=normalize(z)
    fun=OrbitHinge(z.shape,rotations);best=z.copy();bv=float('inf');trace=[]
    for r in range(rounds):
        opt=minimize(fun,z.ravel(),jac=True,method='L-BFGS-B',options=dict(
            maxiter=steps,maxcor=12,maxls=30,ftol=1e-16,gtol=1e-10))
        z=normalize(opt.x.reshape(z.shape))
        trace.append(dict(stage='cyclic843',round=r,objective=float(opt.fun),nit=opt.nit))
        if opt.fun<bv:bv=float(opt.fun);best=z.copy()
        z=best.copy()
        # Concentrate large changes in the last-factor orbits, with core motion allowed.
        concentration=np.sum(z[:,8:]**2,axis=1)
        chosen=rng.choice(np.argsort(concentration)[-20:],size=int(rng.choice([2,4,8])),replace=False)
        z[chosen]+=rng.normal(scale=float(rng.choice([.03,.10,.20])),size=(len(chosen),12))
        z=normalize(z)
    x=np.vstack([best@r for r in rotations]);sym=metrics(x)
    bad=np.maximum(x@x.T-.5,0);np.fill_diagonal(bad,0)
    scores=np.sum(bad*bad,axis=1)
    removed=int(np.argmax(scores))
    x=np.delete(x,removed,axis=0)
    x,opt=relax(x,steps*2);bm=metrics(x);trace.append(dict(stage='released842',**bm,**opt))
    np.save(out/f'orbit_{index:04d}.npy',x)
    mm=bm;minimax=x.copy()
    for power in [4,8,16,32]:
        x,opt=optimize(x,600,power=power);m=metrics(x)
        trace.append(dict(stage='minimax',power=power,**m,**opt))
        if m['max_cos']<mm['max_cos']:minimax=x.copy();mm=m
    np.save(out/f'orbit_{index:04d}_minimax.npy',minimax)
    result=dict(index=index,seed=seed,steps=steps,rounds=rounds,cyclic843=sym,
                removed_index=removed,best=bm,minimax=mm,trace=trace,
                elapsed_seconds=time.monotonic()-start)
    (out/f'orbit_{index:04d}.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=128)
    p.add_argument('--workers',type=int,default=64);p.add_argument('--steps',type=int,default=1800)
    p.add_argument('--rounds',type=int,default=4);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    cfg=vars(a).copy();cfg['out']=str(a.out.resolve())
    cfg['sources_sha256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['orbit_search.py','search.py','continuation.py','kernels.py','quaternion_cells.py']}
    (a.out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    results=[]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(trial,(i,a.steps,a.rounds,str(a.out.resolve()))) for i in range(a.trials)]):
            r=f.result();results.append(r)
            print(json.dumps({k:r[k] for k in ['index','cyclic843','best','minimax','elapsed_seconds']}),flush=True)
            (a.out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
