"""Replace the sixteen algebraic points by normalized rational cavity vertices.

Rational vertex v is feasible against the fixed rational configuration if
A v<=16 and ||v||²>=32. Its radial contraction has norm²32. Two such contractions
are compatible exactly when v.w<=0 or 4(v.w)²<=||v||²||w||². Discovery uses floats;
only independently reconstructed rational vertices can certify an improvement.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import linprog,milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--starts',type=int,default=2048);ap.add_argument('--seed',type=int,default=20261113)
    ap.add_argument('--solve-seconds',type=int,default=90);a=ap.parse_args();start=time.monotonic()
    z=np.load(Path(__file__).with_name('normal26.npz'));Y=z['Y'];A=np.unique(np.vstack([z['B'],Y[:364]]),axis=0)
    assert np.array_equal(A,np.rint(A));A=np.rint(A).astype(np.int32)
    rng=np.random.default_rng(a.seed);vertices={};rows=[];best=0.
    for trial in range(a.starts):
        sigma=[0.,.03,.1,.3,1.][(trial//16)%5];c=Y[364+trial%16]/np.sqrt(32)+sigma*rng.normal(size=10)
        last=0.
        for it in range(30):
            res=linprog(-c,A_ub=A,b_ub=np.full(len(A),16.),bounds=[(None,None)]*10,method='highs',options={'threads':1})
            if not res.success:raise RuntimeError(res.message)
            p=res.x;r=float(p@p)
            if r-last<1e-10:break
            last=r;c=p
        if r>=32-1e-9:
            key=tuple(np.round(p,8));vertices[key]=p
        if r>best+1e-10 or trial%128==0:
            best=max(best,r);row=dict(trial=trial,norm2=r,best_norm2=best,vertices=len(vertices),seconds=time.monotonic()-start);rows.append(row);print(json.dumps(row),flush=True)
    P=np.asarray(list(vertices.values())).reshape((-1,10));r=(P*P).sum(1);n=len(P)
    G=P@P.T;bad=(G>0)&(4*G*G>np.outer(r,r)+1e-6);E=np.argwhere(np.triu(bad,1))
    np.savez_compressed(a.out.with_suffix('.npz'),vertices=P,constraints=A,edges=E,normal_indices=z['indices'])
    M=coo_matrix((np.ones(2*len(E)),(np.repeat(np.arange(len(E)),2),E.ravel())),shape=(len(E),n)).tocsr()
    report=dict(starts=a.starts,seed=a.seed,vertices=n,edges=len(E),best_norm2=best,progress=rows,scope='Floating discovery; exactification required before a lower bound')
    if n:
        res=milp(-np.ones(n),integrality=np.ones(n),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,1),options={'time_limit':a.solve_seconds,'threads':1,'mip_rel_gap':0})
        report.update(solver_status=int(res.status),message=res.message)
        if res.x is not None:
            sel=np.flatnonzero(res.x>.5);report.update(selected=len(sel),selected_indices=sel.tolist(),dual_bound=float(res.mip_dual_bound))
    report['seconds']=time.monotonic()-start;a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='progress'}),flush=True)


if __name__=='__main__':main()
