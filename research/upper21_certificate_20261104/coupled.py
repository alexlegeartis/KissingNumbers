"""Couple label pair moments to degree-three hole moments exactly."""
import argparse,json,time
import numpy as np
import cvxpy as cp
from model import DATA
from moments import constraints,algebra,T,reps,na,nt
from labels_sdp import prepare,four_cliques
from label_outer import linear


def coefficients(target):
    k=1024-target;z=prepare();nv=len(z['reps']);n=nv+nt
    El,el,Al,bl=linear(z)
    Eh,eh,Ah,bh,bounds=constraints(k,False)
    E=np.r_[np.c_[El,np.zeros((len(El),nt))],np.c_[np.zeros((len(Eh),nv)),Eh]]
    e=np.r_[el,eh]
    A=np.r_[np.c_[Al,np.zeros((len(Al),nt))],np.c_[np.zeros((len(Ah),nv)),k*Ah],
            np.c_[np.zeros((nt,nv)),k*np.eye(nt)],np.c_[np.zeros((nt,nv)),-k*np.eye(nt)]]
    b=np.r_[bl,k*bh,np.zeros(nt),-k*np.array([r[1] for r in bounds])]
    # Complement-conditioned clique inequality.
    A=np.r_[A,np.c_[target*z['cliques'],np.zeros((len(z['cliques']),nt))]]
    b=np.r_[b,np.full(len(z['cliques']),23*target-21504)]
    v=np.zeros(n);v[:nv]=z['row'];E=np.r_[E,v[None,:]];e=np.r_[e,target]
    for d in reps:
        v=np.zeros(n);v[:nv]=-target*np.bincount(z['orbits'][d].ravel(),minlength=nv)
        v[nv+T[0,d]]=21*k;E=np.r_[E,v[None,:]];e=np.r_[e,21*(2*k-1024)]
    for arr in (E,e,A,b):assert np.max(abs(arr-np.rint(arr)))<1e-7
    return tuple(np.rint(arr).astype(np.int64) for arr in (E,e,A,b))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--target',type=int,default=971);a=p.parse_args()
    start=time.monotonic();z=prepare();nv=len(z['reps']);alg=algebra();x=cp.Variable(nv+nt)
    E,e,A,b=coefficients(a.target)
    con=[E@x==e,A@x>=b]
    for B in z['blocks']:con.append(cp.reshape(B.reshape(441,nv)@x[:nv],(21,21),order='C')>>0)
    for B in (alg['LT'],alg['LX']-alg['LT']):
        con.append(cp.reshape(B.reshape(nt,-1).T@x[nv:],(na,na),order='C')>>0)
    prob=cp.Problem(cp.Minimize(0),con)
    prob.solve(solver='CLARABEL',max_iter=150,tol_feas=1e-9,tol_gap_abs=1e-9,tol_gap_rel=1e-9)
    name=f'coupled_target{a.target}'
    report=dict(run=name,target=a.target,status=prob.status,seconds=time.monotonic()-start,
                variables=nv+nt,equalities=len(E),inequalities=len(A),numerical_only=True)
    if x.value is not None:
        np.save(DATA/(name+'_moments.npy'),x.value)
        report.update(equality_error=float(np.max(abs(E@x.value-e))),minimum_linear_slack=float(np.min(A@x.value-b)))
    (DATA/(name+'_result.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
