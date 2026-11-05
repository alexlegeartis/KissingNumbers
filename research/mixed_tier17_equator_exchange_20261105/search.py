"""Jointly remove source equator/cap orbits and add transformed equator orbits."""
import numpy as np,time,os
from scipy.sparse import csr_matrix,coo_matrix
from scipy.optimize import milp,Bounds,LinearConstraint
z=np.load('research/mixed_tier17_equator_exchange_20261105/cache/pool.npz');X=z['X'].astype(np.int16);pts=z['points'].astype(np.int16);sizes=z['sizes'];assert np.all(sizes==4);Q=pts.reshape(-1,4,17);rep=Q[:,0,:16].astype(np.int32);weights=z['weights'];point_group=z['point_group'];indices=z['indices'];indptr=z['indptr'];n=len(Q);B=csr_matrix((np.ones(len(indices),np.int16),indices,indptr),shape=(n,len(weights)));cost=z['cost'];rng=np.random.default_rng(171024);H=np.fromfile('research/mixed_tier17_eval_subgroup_20261105/cache/mixed17_eval4_group.bin',np.uint16);K=[int(h) for h in H if not int(h)&6];sign=np.array([[1-2*(h>>j&1) for j in range(16)] for h in K],np.int32);cand=np.flatnonzero(cost<=56);print('candidates',n,'low-cost',len(cand),'source groups',len(weights),flush=True);best=-999;start=time.time()
for trial in range(int(os.getenv('EQUATOR_TRIALS','10000'))):
 seed=int(rng.choice(cand));vict=np.zeros(len(weights),bool);vict[B.indices[B.indptr[seed]:B.indptr[seed+1]]]=1
 for step in range(25):
  covered=np.asarray(B@vict.astype(np.int16)).ravel();missing=np.diff(B.indptr)-covered;eligible=np.flatnonzero(missing==0)
  if step%3==0 and len(eligible):
   if len(eligible)>350:eligible=np.r_[seed,rng.choice(eligible[eligible!=seed],349,False)]
   elif seed not in eligible:eligible=np.r_[seed,eligible]
   U=rep[eligible];D=np.zeros((len(U),len(U)),bool)
   for sg in sign:D|=(U*sg)@U.T>144
   np.fill_diagonal(D,False);i,j=np.where(np.triu(D,1));src=np.unique(B[eligible].indices);bm=B[eligible][:,src].tocoo();n1=len(eligible);nb=len(bm.data);ne=len(i);rows=np.r_[np.repeat(np.arange(nb),2),np.repeat(np.arange(nb,nb+ne),2)];cols=np.r_[np.c_[bm.row,n1+bm.col].ravel(),np.c_[i,j].ravel()];vals=np.r_[np.tile([1.,-1.],nb),np.ones(2*ne)];A=coo_matrix((vals,(rows,cols)),shape=(nb+ne,n1+len(src))).tocsr();lower=np.zeros(n1+len(src));lower[np.flatnonzero(eligible==seed)[0]]=1
   res=milp(np.r_[-4*np.ones(n1),weights[src]],integrality=np.ones(n1+len(src)),bounds=Bounds(lower,np.ones(n1+len(src))),constraints=LinearConstraint(A,-np.inf,np.r_[np.zeros(nb),np.ones(ne)]),options={'time_limit':1})
   if res.fun is not None:
    picked=eligible[res.x[:n1]>.5];removed=src[res.x[n1:]>.5];gain=4*len(picked)-int(weights[removed].sum())
    if gain>best:best=gain;print('best',gain,'added groups',len(picked),'removed points',int(weights[removed].sum()),'eligible',len(eligible),'step',step,'trial',trial,'secs',round(time.time()-start,1),flush=True)
    if gain>0:
     keep=~np.isin(point_group,removed);Y=np.r_[X[keep],Q[picked].reshape(-1,17)];np.savez('research/mixed_tier17_equator_exchange_20261105/cache/POTENTIAL.npz',X=Y,g=np.r_[np.ones(16,dtype=np.int16),2],norm=288);print('POTENTIAL',len(Y),'requires independent full-coordinate check',flush=True);raise SystemExit
  options=cand[(missing[cand]>0)&(missing[cand]<=8)]
  if not len(options):break
  if len(options)>1000:options=rng.choice(options,1000,False)
  scores=covered[options]/(missing[options]+.5);top=scores.max();options=options[scores>=top*.8];j=int(rng.choice(options));vict[B.indices[B.indptr[j]:B.indptr[j+1]]]=1
  if weights[vict].sum()>500:break
 if trial%100==0:print('trial',trial,'best',best,'secs',round(time.time()-start,1),flush=True)
print('DONE best',best,flush=True)
