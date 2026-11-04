#!/usr/bin/env python3
"""Collectively deform N+1 points. Floating candidates are never claimed as certificates."""
import argparse
import itertools
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import linprog, minimize
from orbits17 import DATA


def baseline(dimension):
    if dimension == 17:
        z=np.load(DATA/'pool17.npz')
        E=np.pad(z['equator']/6.,((0,0),(0,1)))
        V=z['vertices'][z['seed']].astype(float)/6.
        V[:,-1]*=np.sqrt(2.)
        poles=np.zeros((2,17));poles[:,-1]=[np.sqrt(8.),-np.sqrt(8.)]
        return np.concatenate([E,V,poles])/np.sqrt(8.)
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'common'))
    from golay import golay24
    G=golay24(); B=G[(G[:,:4].sum(1)==0)&(G.sum(1)==8),4:]
    rows=[]
    for i,j in itertools.combinations(range(20),2):
        for si,sj in itertools.product((-2,2),repeat=2):
            v=np.zeros(20);v[i]=si;v[j]=sj;rows.append(v)
    for block in B:
        positions=np.flatnonzero(block)
        for bits in range(256):
            if bits.bit_count()%2:
                v=np.zeros(20)
                for j,p in enumerate(positions):v[p]=-1 if bits>>j&1 else 1
                rows.append(v)
    C=np.asarray(rows)/np.sqrt(8.)
    F=(1.-2.*G[G[:,0]==0,4:])/np.sqrt(20.)
    X=np.concatenate([C,F])
    assert X.shape==(19448,20)
    return X


def run(dimension,seed,seconds,restarts):
    start=time.monotonic();rng=np.random.default_rng(seed)
    X=baseline(dimension)
    bestnorm=0.;hole=None
    # LP ascent supplies an insertion direction near a hole, then all points move.
    for restart in range(restarts):
        v=rng.normal(size=dimension);v/=np.linalg.norm(v)
        for iteration in range(12):
            result=linprog(-v,A_ub=X,b_ub=np.full(len(X),.5),bounds=[(-2,2)]*dimension,
                           method='highs',options={'threads':1})
            if not result.success:break
            w=result.x; norm=np.linalg.norm(w);direction=w/norm
            if norm>bestnorm and float((X@w).max())<=.5+1e-9:
                bestnorm=norm;hole=direction.copy()
            if np.linalg.norm(direction-v)<1e-8:break
            v=direction
    assert hole is not None
    name=f'continuous{dimension}_s{seed}'
    print(json.dumps(dict(event='hole',dimension=dimension,seed=seed,
                          insertion_max_dot=float((X@hole).max()),seconds=time.monotonic()-start)),flush=True)
    V=np.concatenate([X,hole[None,:]])
    if seed%2:V+=rng.normal(size=V.shape)*1e-4
    shape=V.shape; evaluations=0;bestenergy=float('inf');best=None;bestmax=None;last=start
    threshold=.5-1e-10
    class TimeLimit(Exception):pass
    def objective(flat):
        nonlocal evaluations,bestenergy,best,bestmax,last
        if time.monotonic()-start>seconds:raise TimeLimit()
        V=flat.reshape(shape);norms=np.linalg.norm(V,axis=1);U=V/norms[:,None]
        H=U@U.T
        np.fill_diagonal(H,0.)
        maxdot=float(H.max())
        H-=threshold;np.maximum(H,0.,out=H)
        energy=.5*float(np.einsum('ij,ij->',H,H))
        grad=2.*(H@U)
        grad-=(grad*U).sum(1)[:,None]*U
        grad/=norms[:,None]
        evaluations+=1
        if energy<bestenergy:
            bestenergy=energy;best=U.copy();bestmax=maxdot
        now=time.monotonic()
        if now-last>30:
            print(json.dumps(dict(event='deformation',dimension=dimension,seed=seed,
                                  energy=bestenergy,max_dot=bestmax,evaluations=evaluations,
                                  seconds=now-start)),flush=True)
            last=now
        return energy,grad.ravel()
    status='TIME_LIMIT'
    try:
        result=minimize(objective,V.ravel(),jac=True,method='L-BFGS-B',
                        options=dict(maxiter=3000,maxfun=5000,ftol=1e-15,gtol=1e-10,maxcor=10))
        status=str(result.message)
    except TimeLimit:pass
    if best is not None:np.save(DATA/(name+'_candidate.npy'),best)
    report=dict(dimension=dimension,count=len(V),seed=seed,seconds=time.monotonic()-start,
                evaluations=evaluations,best_energy=bestenergy,max_dot=bestmax,status=status,
                certified=False,scope='Numerical exploration only; any record needs an exact subsequent certificate.')
    (DATA/(name+'_result.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dimension',type=int,choices=(17,20),required=True)
    p.add_argument('--seed',type=int,default=1);p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--hole-restarts',type=int,default=12);a=p.parse_args()
    run(a.dimension,a.seed,a.seconds,a.hole_restarts)
