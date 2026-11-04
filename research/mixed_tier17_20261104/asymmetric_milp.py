import os
import numpy as np,json
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
R='research/mixed_tier17_20261104/';F=np.fromfile(R+'mixed17_equator_frontier.bin',np.int32).reshape(-1,6);n=2*len(F);w=np.fromfile(R+'mixed17_asymmetric_weights.bin',np.int32);E=np.fromfile(R+'mixed17_asymmetric_edges.bin',np.int32).reshape(-1,2);assert len(w)==n+568
rr=np.repeat(np.arange(len(E)),2);val=np.ones(len(rr));old=E[:,1]>=n;val[1::2][old]=-1;A=coo_matrix((val,(rr,E.ravel())),shape=(len(E),len(w))).tocsr();obj=np.r_[-w[:n],w[n:]].astype(float);print(json.dumps({'candidate_atoms':n,'old_atoms':568,'constraints':len(E)}),flush=True)
r=milp(obj,integrality=np.ones(len(w)),bounds=Bounds(0,1),constraints=LinearConstraint(A,-np.inf,np.where(old,0.,1.)),options={'time_limit':float(os.environ.get('MIXED17_MILP_SECONDS','900'))})
print(json.dumps({'status':r.message,'gain':None if r.fun is None else -r.fun,'upper_gain':None if r.mip_dual_bound is None else -r.mip_dual_bound}),flush=True)
if r.x is not None and r.fun<-.5:np.savez(R+'mixed17_asymmetric_candidate.npz',new_indices=np.flatnonzero(r.x[:n]>.5),removed_indices=np.flatnonzero(r.x[n:]>.5),gain=-r.fun)
