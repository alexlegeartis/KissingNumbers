"""Weighted shared-blocker LNS with both cap and equator additions."""
import numpy as np,time,os
from scipy.sparse import csr_matrix,coo_matrix
from scipy.optimize import milp,Bounds,LinearConstraint
w=np.fromfile('research/mixed_tier17_root_sums_20261105/cache/weights.bin',np.uint8).astype(np.int16);E=np.fromfile('research/mixed_tier17_root_sums_20261105/cache/edges.bin',np.int32).reshape(-1,2);ng=int(np.load('research/mixed_tier17_root_sums_20261105/cache/nodes.npz')['ng']);n=len(w)-ng
cross=(E[:,0]<ng)&(E[:,1]>=ng);B=coo_matrix((np.ones(cross.sum(),np.int16),(E[cross,1]-ng,E[cross,0])),shape=(n,ng)).tocsr();C=E[(E[:,0]>=ng)&(E[:,1]>=ng)]-ng;A=coo_matrix((np.ones(2*len(C),np.uint8),(np.r_[C[:,0],C[:,1]],np.r_[C[:,1],C[:,0]])),shape=(n,n)).tocsr();cost=np.asarray(B@w[:ng]).ravel();low=np.flatnonzero(cost<=56);rng=np.random.default_rng(170922);best=-999;maxeligible=0;start=time.time();print('source',ng,'candidates',n,'eligible',len(low),'edges',len(E),flush=True)
for trial in range(int(os.getenv('ROOT_TRIALS','30000'))):
 seed=int(rng.choice(low));vict=np.zeros(ng,bool);vict[B.indices[B.indptr[seed]:B.indptr[seed+1]]]=1
 for step in range(35):
  hit=np.asarray(B@vict.astype(np.int16)).ravel();missing=np.diff(B.indptr)-hit;eligible=np.flatnonzero(missing==0)
  if len(eligible)>maxeligible:maxeligible=len(eligible);print('max eligible',maxeligible,'trial',trial,'step',step,flush=True)
  if step%4==0 and len(eligible):
   if len(eligible)>400:eligible=np.r_[seed,rng.choice(eligible[eligible!=seed],399,False)]
   elif seed not in eligible:eligible=np.r_[seed,eligible]
   Csub=A[eligible][:,eligible].tocoo();good=Csub.row<Csub.col;i=Csub.row[good];j=Csub.col[good];src=np.unique(B[eligible].indices);bi=B[eligible][:,src].tocoo();n1=len(eligible);nb=len(bi.data);ne=len(i);rows=np.r_[np.repeat(np.arange(nb),2),np.repeat(np.arange(nb,nb+ne),2)];cols=np.r_[np.c_[bi.row,n1+bi.col].ravel(),np.c_[i,j].ravel()];vals=np.r_[np.tile([1.,-1.],nb),np.ones(2*ne)];M=coo_matrix((vals,(rows,cols)),shape=(nb+ne,n1+len(src))).tocsr();lb=np.zeros(n1+len(src));lb[np.flatnonzero(eligible==seed)[0]]=1
   res=milp(np.r_[-w[ng+eligible].astype(float),w[src]],integrality=np.ones(n1+len(src)),bounds=Bounds(lb,np.ones(n1+len(src))),constraints=LinearConstraint(M,-np.inf,np.r_[np.zeros(nb),np.ones(ne)]),options={'time_limit':1})
   if res.fun is not None:
    picked=eligible[res.x[:n1]>.5];removed=src[res.x[n1:]>.5];gain=int(w[ng+picked].sum())-int(w[removed].sum())
    if gain>best:best=gain;print('best',gain,'added',len(picked),'removed points',int(w[removed].sum()),'eligible',len(eligible),'trial',trial,'step',step,'seconds',round(time.time()-start,1),flush=True)
    if gain>0:
     mask=np.r_[np.ones(ng,np.uint8),np.zeros(n,np.uint8)];mask[removed]=0;mask[ng+picked]=1;assert np.all(mask[E[:,0]]+mask[E[:,1]]<=1);assert int(mask@w.astype(np.int64))==5730+gain;mask.tofile('research/mixed_tier17_root_sums_20261105/cache/POTENTIAL.bin');print('POTENTIAL',5730+gain,'requires independent coordinate replay',flush=True);raise SystemExit
  opts=low[(missing[low]>0)&(missing[low]<=10)]
  if not len(opts):break
  if len(opts)>1500:opts=rng.choice(opts,1500,False)
  scores=hit[opts]/(missing[opts]+.5);top=scores.max();opts=opts[scores>=top*.8];k=int(rng.choice(opts));vict[B.indices[B.indptr[k]:B.indptr[k+1]]]=1
  if w[:ng][vict].sum()>500:break
 if trial%100==0:print('trial',trial,'best',best,'seconds',round(time.time()-start,1),flush=True)
print('DONE best',best,flush=True)
