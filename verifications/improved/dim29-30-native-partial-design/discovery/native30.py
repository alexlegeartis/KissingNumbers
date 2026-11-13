import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import sys,json,numpy as np
from pathlib import Path
from itertools import product
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import csr_matrix,coo_matrix,hstack,vstack
O=Path(__file__).parent;R=O.parent/'kn_latest';D=R/'verifications/improved/dim30-leech-sloane-section/data'
sys.path.insert(0,str(R/'research/reclaim_26_31_20261113/algebraic'))
from leech_sloane_module import model
import importlib.util
spec=importlib.util.spec_from_file_location('vm',D.parent/'verify_module.py');vm=importlib.util.module_from_spec(spec);spec.loader.exec_module(vm)
recipe=json.load(open(D/'construction.json'));parent=vm.reconstruct(recipe)[0]
blocks,A=model();ids=[i for i,b in enumerate(blocks) if max(b)<12];A=A[:,ids];blocks=[blocks[i] for i in ids]
un=[]
for b in blocks:
 for signs in product([-2,2],repeat=4):
  w=np.zeros(15,np.int16);w[list(b)]=signs;z=np.zeros(14,np.int16)
  for i in range(7):z[i]=w[2*i]+w[2*i+1]
  for i in range(6):z[8+i]=w[2*i]-w[2*i+1]
  un.append(z)
un=np.array(un);oldids=[i for i,b in enumerate(blocks) if list(b) in recipe['transversal_blocks']];print('oldids',len(oldids),flush=True)
keys={r.tobytes() for r in un.reshape(-1,16,14)[oldids].reshape(-1,14)}
section=np.delete(parent[parent[:,14]==0],14,axis=1).astype(np.int16);fixed=np.array([r for r in section if r.tobytes() not in keys]);print('fixed',len(fixed),'variable',len(un),flush=True)
mo=np.load(D/'module.npz');rot=mo['transformation_numerator'];den=int(mo['transformation_denominator']);un=(un@rot.T//den).astype(np.int16);fixed=(fixed@rot.T//den).astype(np.int16)
assert np.all(np.sum(un.astype(int)**2,axis=1)==32)
g=un.astype(int)@fixed.T.astype(int);bad=np.any(g.reshape(-1,16,len(fixed))>16,axis=(1,2));print('fixed incompatible blocks',sum(bad),flush=True)
x=np.load(D/'coordinates30_248084.npz')['points'].astype(np.int16)//2;E=np.r_[[0,2,4,5,6,10,11,23],np.arange(24,30)];out=np.setdiff1d(np.arange(30),E);outer=x[np.any(x[:,out],axis=1)];pr,counts=np.unique(outer[:,E],axis=0,return_counts=True)
forced=np.any(fixed.astype(np.float32)@pr.T.astype(np.float32)>16,axis=0)
cross=np.any((un.astype(np.float32)@pr.T.astype(np.float32)>16).reshape(-1,16,len(pr)),axis=1);cross[:,forced]=False;used=np.flatnonzero(cross.any(axis=0));cross=cross[:,used]
r,c=np.nonzero(cross);B=len(blocks);nr=len(r);mat=coo_matrix((np.tile([1,-1],nr),(np.repeat(np.arange(nr),2),np.c_[r,B+c].ravel())),shape=(nr,B+len(used))).tocsr();con=vstack([hstack([A,csr_matrix((280,len(used)))]),mat]);cost=np.r_[-16*np.ones(B),counts[used]];ub=np.ones(len(cost));ub[:B][bad]=0
res=milp(cost,integrality=np.ones(len(cost)),bounds=Bounds(np.zeros(len(cost)),ub),constraints=LinearConstraint(con,np.full(con.shape[0],-np.inf),np.r_[np.ones(280),np.zeros(nr)]),options={'time_limit':120,'threads':2})
print(res.message,res.fun,'bound',res.mip_dual_bound,flush=True)
if res.x is not None:
 sel=np.flatnonzero(res.x[:B]>.5);y=np.r_[fixed,un.reshape(-1,16,14)[sel].reshape(-1,14)];G=y.astype(int)@y.T.astype(int);np.fill_diagonal(G,-32);assert G.max()<=16
 np.savez_compressed(O/'native30module.npz',points=y,selected=sel)
 total=len(outer)+len(y)-int(counts[forced].sum())-int(counts[used[res.x[B:]>.5]].sum());print('TOTAL',total,'module',len(y),'blocks',len(sel),flush=True)
 (O/'native30.json').write_text(json.dumps({'count':total,'blocks':len(sel),'objective':res.fun,'bound':res.mip_dual_bound,'status':res.message}))
