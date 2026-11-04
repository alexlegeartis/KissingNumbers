import os
import numpy as np,json
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
R='research/mixed_tier17_20261104/';rows=np.fromfile(R+'mixed17_equator_frontier.bin',np.int32).reshape(-1,6);n=len(rows);w=np.fromfile(R+'mixed17_switch_weights.bin',np.int32);E=np.fromfile(R+'mixed17_switch_edges.bin',np.int32).reshape(-1,2);assert len(w)==n+536
rr=np.repeat(np.arange(len(E)),2);cc=E.ravel();val=np.ones(len(rr));old=E[:,1]>=n;val[1::2][old]=-1;ub=np.where(old,0.,1.);A=coo_matrix((val,(rr,cc)),shape=(len(E),len(w))).tocsr();obj=np.r_[-w[:n],w[n:]].astype(float)
limit=float(os.environ.get('MIXED17_MILP_SECONDS','900'))
print(json.dumps({'candidates':n,'seed_atoms':536,'constraints':len(E),'time_limit_seconds':limit}),flush=True)
r=milp(obj,integrality=np.ones(len(w)),bounds=Bounds(0,1),constraints=LinearConstraint(A,-np.inf,ub),options={'time_limit':limit})
print(json.dumps({'status':r.message,'best_gain':None if r.fun is None else -r.fun,'upper_gain':None if r.mip_dual_bound is None else -r.mip_dual_bound}),flush=True)
if r.x is not None and r.fun<-.5:
 np.savez(R+'mixed17_equator_switch_candidate.npz',new_frontier_indices=np.flatnonzero(r.x[:n]>.5),removed_seed_indices=np.flatnonzero(r.x[n:]>.5),gain=-r.fun);print('POTENTIAL',flush=True)
