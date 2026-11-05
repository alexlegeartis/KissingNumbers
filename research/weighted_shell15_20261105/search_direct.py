import numpy as np,sys
from pathlib import Path
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
root=Path(__file__).parent;z=np.load(root/'source.npz');X=z['X'].astype(np.int16)*12;g=z['g'].astype(np.int16);C=np.load(root/'sparse_low.npz')['Q'].astype(np.int16);n=len(C);print('nonpure candidates',n,flush=True)
D=(C*g)@X.T;ii,jj=np.where(D>1152);ne=len(ii);G=(C*g)@C.T;i,j=np.where(np.triu(G>1152,1));e=len(i);print('source edges',ne,'candidate edges',e,flush=True)
rows=np.r_[np.repeat(np.arange(ne,dtype=np.int32),2),np.repeat(np.arange(ne,ne+e,dtype=np.int32),2),np.full(n+len(X),ne+e,np.int32)];cols=np.r_[np.column_stack((ii,n+jj)).ravel(),np.column_stack((i,j)).ravel(),np.arange(n+len(X))].astype(np.int32);vals=np.r_[np.tile([1.,-1.],ne),np.ones(2*e),-np.ones(n),np.ones(len(X))];M=coo_matrix((vals,(rows,cols)),shape=(ne+e+1,n+len(X))).tocsr()
if '--build-only' in sys.argv: raise SystemExit(0)
seconds=float(sys.argv[sys.argv.index('--seconds')+1]) if '--seconds' in sys.argv else 900
res=milp(np.r_[-np.ones(n),np.ones(len(X))],integrality=np.ones(n+len(X)),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,np.r_[np.zeros(ne),np.ones(e),-1.]),options={'time_limit':seconds});print('status',res.message,'gain',None if res.fun is None else -res.fun,'dual',res.mip_dual_bound,flush=True)
if res.fun is not None and res.fun<-.5:
 ids=np.flatnonzero(res.x[:n]>.5);rm=np.flatnonzero(res.x[n:]>.5);Y=np.r_[np.delete(X,rm,axis=0),C[ids]];np.savez_compressed('candidate_direct.npz',X=Y,g=g,den=12);print('POTENTIAL',len(Y),flush=True)
