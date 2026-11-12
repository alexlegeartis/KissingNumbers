"""Optimize exact rational two-point split proposals and shared deletions.

The optimizer proposes an exchange. All returned combinatorial constraints
are checked independently before a delta is written.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import numpy as np,json,time,sys
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
HERE=Path(__file__).parent

def run(seconds=120,extra=False):
    suffix='_extra' if extra else ''
    p=json.loads((HERE/f'medium43_pool_f59{suffix}.json').read_text())
    o=json.loads((HERE/f'medium43_options{suffix}.json').read_text());f=np.load(HERE/f'medium43_pool_f59{suffix}.npz')
    V=p['vertices'];options=o['options'];B=f['selfbeta'];K=f['paircut'];n=len(options)
    blockers=sorted(set(b for _,_,bs in options for b in bs));index={b:i+n for i,b in enumerate(blockers)}
    rows=[];cols=[];data=[];upper=[]
    def add(items,rhs):
        row=len(upper)
        for col,val in items:rows.append(row);cols.append(col);data.append(val)
        upper.append(rhs)
    for i,(_,_,bs) in enumerate(options):
        for b in bs:add([(i,1),(index[b],-1)],0)
    conflicts=[]
    for i,(a,b,_) in enumerate(options):
        r=V[a]['row'];heads=[V[a]['head'],V[b]['head']]
        for j,(c,d,_) in enumerate(options[:i]):
            s=V[c]['row'];hh=[V[c]['head'],V[d]['head']]
            if r==s or any(B[r,s]>K[x,y] for x in heads for y in hh):
                add([(i,1),(j,1)],1);conflicts.append((i,j))
    A=coo_matrix((data,(rows,cols)),shape=(len(upper),n+len(blockers))).tocsc()
    objective=np.r_[np.full(n,-2.),np.ones(len(blockers))]
    print('MIP',A.shape,'nonzero',A.nnz,'conflicts',len(conflicts),flush=True)
    start=time.time()
    sol=milp(objective,integrality=np.ones(len(objective)),bounds=Bounds(0,1),
      constraints=LinearConstraint(A,np.full(len(upper),-np.inf),np.array(upper)),
      options=dict(time_limit=seconds,mip_rel_gap=0,threads=1))
    print('STATUS',sol.message,'objective',sol.fun,'bound',getattr(sol,'mip_dual_bound',None),'seconds',time.time()-start,flush=True)
    if sol.x is None:return
    chosen=np.flatnonzero(sol.x[:n]>.5).tolist();chosen_set=set(chosen)
    assert all(i not in chosen_set or j not in chosen_set for i,j in conflicts)
    newids=sorted(set(v for i in chosen for v in options[i][:2]));deleted=sorted(set(b for i in chosen for b in options[i][2]))
    assert len(newids)==2*len(chosen)
    for i in newids:
        assert set(V[i]['blockers'])<=set(deleted)
        for j in newids:
            if i!=j:assert B[V[i]['row'],V[j]['row']]<=K[V[i]['head'],V[j]['head']]
    result=dict(status='Exact combinatorial exchange; fresh geometric replay still required',gain=len(newids)-len(deleted),
      new_vertex_indices=newids,deleted_current_indices=deleted,option_indices=chosen,solver_message=sol.message,
      solver_upper_bound=-float(sol.mip_dual_bound),seconds=time.time()-start)
    (HERE/f'medium43_exchange{suffix}.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':run(float(sys.argv[1]) if len(sys.argv)>1 else 120,'--extra' in sys.argv)
