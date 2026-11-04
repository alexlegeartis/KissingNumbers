"""Global degree-three necessary conditions for a dominating set of size k.

Moment m(a,b) = E[h(0)h(a)h(b)] / E[h(0)], averaged over verified
translations and graph automorphisms. All reductions preserve arbitrary H.
"""
import argparse,json,time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from model import *

Z=np.load(DATA/'moments.npz');O=Z['word_orbit'];T=Z['triple_orbit'];A=Z['pair_orbit']
TR=Z['triple_reps'];PR=Z['pair_reps'];F=Z['fourier']
nt=len(TR);na=len(PR);reps=np.array([np.flatnonzero(O==i)[0] for i in range(8)])
X=T[0,reps];sizes=np.bincount(O);I=np.eye(nt)


def constraints(k,independent=False):
    eq=[];eb=[];ge=[];gb=[]
    def equality(v,b):eq.append(v);eb.append(b)
    def lower(v,b):ge.append(v);gb.append(b)
    equality(I[T[0,0]],1)
    equality(np.bincount(X,weights=sizes,minlength=nt),k)
    for d in reps:
        equality(np.bincount(T[d],minlength=nt)-k*I[T[0,d]],0)
    for row in F:
        lower(np.bincount(X,weights=row,minlength=nt),1/k if k%2 else 0)
    for t,rep in enumerate(TR):
        d,e=divmod(int(rep),1024)
        pairs=[T[0,d],T[0,e],T[0,d^e]]
        for j in pairs:lower(I[j]-I[t],0)
        for j in range(3):lower(I[t]-I[pairs[j]]-I[pairs[(j+1)%3]],-1)
        lower(sum((I[j] for j in pairs),np.zeros(nt))-I[t],3-1024/k)
    for rep in PR:
        u,d=divmod(int(rep),1024);V=np.r_[u,N[u]]
        st=np.bincount(T[d,V],minlength=nt)
        s0=np.bincount(T[0,V],minlength=nt)
        sd=np.bincount(T[0,V^d],minlength=nt)
        xd=I[T[0,d]]
        lower(st-xd,0)
        lower(s0-st+xd,1)
        lower(sd-st+xd,1)
        lower(st-s0-sd-xd,1024/k-24)
    if independent:equality(I[T[0,P[0]]],0)
    bounds=[(0,1)]*nt
    if k%2:
        for j in X[1:]:bounds[j]=(0,(k-1)/k)
    return np.array(eq),np.array(eb),np.array(ge),np.array(gb),bounds


def algebra():
    path=DATA/'algebra.npz'
    if path.exists():return np.load(path)
    osize=np.bincount(A.ravel()).astype(np.int64)
    GT=np.zeros((nt,na,na),dtype=np.int64);GX=GT.copy()
    W=np.arange(1024)
    for i,rep in enumerate(PR):
        u,v=divmod(int(rep),1024)
        for result,ids in ((GT,T[u]),(GX,T[0,u^W])):
            counts=np.bincount(ids*na+A[:,v],minlength=nt*na).reshape(nt,na)
            result[:,i,:]=counts*osize[i]
    assert np.array_equal(GT,GT.transpose(0,2,1))
    assert np.array_equal(GX,GX.transpose(0,2,1))
    scale=np.sqrt(osize[:,None]*osize[None,:])
    np.savez_compressed(path,GT=GT,GX=GX,LT=GT/scale,LX=GX/scale,orbit_sizes=osize)
    return np.load(path)


def solve(k,independent,sdp):
    start=time.monotonic();E,e,G,g,bounds=constraints(k,independent)
    name=f'moments_k{k}_'+('ind' if independent else 'dom')+('_sdp' if sdp else '_lp')
    if not sdp:
        result=linprog(np.zeros(nt),A_ub=-G,b_ub=-g,A_eq=E,b_eq=e,bounds=bounds,method='highs')
        status=result.message;solution=result.x
        report=dict(run=name,status_code=int(result.status),status=status,seconds=time.monotonic()-start)
    else:
        import cvxpy as cp
        z=cp.Variable(nt);alg=algebra()
        lt=alg['LT'].reshape(nt,-1).T;lx=alg['LX'].reshape(nt,-1).T
        constraints_cp=[E@z==e,G@z>=g,z>=np.array([b[0] for b in bounds]),z<=np.array([b[1] for b in bounds])]
        constraints_cp+=[cp.reshape(lt@z,(na,na),order='C')>>0,
                         cp.reshape((lx-lt)@z,(na,na),order='C')>>0]
        prob=cp.Problem(cp.Minimize(0),constraints_cp)
        prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=150,verbose=False)
        solution=z.value
        report=dict(run=name,status=prob.status,solver_seconds=prob.solver_stats.solve_time,seconds=time.monotonic()-start)
    report.update(holes=k,independent=independent,variables=nt,equalities=len(E),inequalities=len(G),
                  numerical_only=True)
    if solution is not None:
        np.save(DATA/(name+'_moments.npy'),solution)
        report.update(equality_error=float(np.max(abs(E@solution-e))),minimum_linear_slack=float(np.min(G@solution-g)))
        if sdp:
            tmat=np.einsum('i,ijk->jk',solution,alg['LT'])
            xmat=np.einsum('i,ijk->jk',solution,alg['LX'])
            report['minimum_psd_eigenvalues']=[float(np.linalg.eigvalsh(tmat)[0]),float(np.linalg.eigvalsh(xmat-tmat)[0])]
    (DATA/(name+'_result.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--k',type=int,default=63);p.add_argument('--independent',action='store_true');p.add_argument('--sdp',action='store_true')
    a=p.parse_args();solve(a.k,a.independent,a.sdp)
