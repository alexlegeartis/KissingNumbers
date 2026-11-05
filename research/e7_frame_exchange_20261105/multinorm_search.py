"""Search mixed pointwise exchanges in norm-3, norm-8 and norm-18 shells."""
from pathlib import Path
import numpy as np,time,os
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack,hstack,csr_matrix
p=Path(__file__).resolve().parent/'cache';z=np.load(p/'multinorm9.npz');Q=z['Q'].astype(np.int16);B=z['B'];norm=z['norm'].astype(np.int16);X=z['X'].astype(np.int16);w=np.load(p/'candidates14.npz');U=w['Q'].astype(np.int16);UB=w['B'];old=w['old'].astype(np.int16);sel=(UB!=65535).sum(1)==4;U=U[sel];UB=UB[sel,:4];V=np.r_[U,-U];D=V@X.T;bad=D>4;assert np.all(bad.sum(1)==4);VB=np.full((len(V),9),65535,np.uint16)
for i in range(len(V)):VB[i,:4]=np.flatnonzero(bad[i])
Q=np.r_[Q,V];B=np.r_[B,VB];norm=np.r_[norm,np.full(len(V),8,np.int16)];cost=(B!=65535).sum(1);M=len(Q);N=len(X);dot=Q@Q.T;threshold=np.floor(np.sqrt(norm[:,None].astype(float)*norm[None,:])/2).astype(np.int16);i,j=np.where(np.triu(dot>threshold,1));a,b=np.where(B!=65535);ne=len(a);e=len(i);rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)];cols=np.r_[np.c_[a,M+B[a,b]].ravel(),np.c_[i,j].ravel()];vals=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)];A=coo_matrix((vals,(rows,cols)),shape=(ne+e,M+N)).tocsr();upper=np.r_[np.zeros(ne),np.ones(e)]
# The exact multinorm verifier proves this local inequality for all old points.
cheap=(norm==3)|(cost==4);affected=np.unique(B[cheap]);affected=affected[affected!=65535];rr=[];cc=[];vv=[]
for k,v in enumerate(affected):
 members=np.flatnonzero(cheap&(B==v).any(axis=1));rr.extend([k]*(len(members)+1));cc.extend(members.tolist()+[M+int(v)]);vv.extend(np.where(norm[members]==3,4.,3.).tolist()+[-12.])
cut=coo_matrix((vv,(rr,cc)),shape=(len(affected),M+N)).tocsr();A=vstack([A,cut],format='csr');upper=np.r_[upper,np.zeros(len(affected))]
# Each nine-conflict candidate hits six points untouched by cheap classes;
# the exact verifier bounds the compatible high candidates per such point by 12.
high=np.flatnonzero(cost==9);free=np.unique(B[high]);free=np.array([v for v in free if v!=65535 and v not in set(affected)],np.int32);assert len(free)==288
rr=[];cc=[];vv=[]
for k,v in enumerate(free):
 members=high[(B[high]==v).any(axis=1)];assert len(members)==72
 rr.extend([k]*(len(members)+1));cc.extend(members.tolist()+[M+int(v)]);vv.extend([1.]*len(members)+[-12.])
cut=coo_matrix((vv,(rr,cc)),shape=(len(free),M+N)).tocsr();A=vstack([A,cut],format='csr');upper=np.r_[upper,np.zeros(len(free))]
c=np.r_[-np.ones(M),np.ones(N)];gain=coo_matrix(c.reshape(1,-1)).tocsr();A=vstack([A,gain],format='csr');upper=np.r_[upper,-1.];high=np.flatnonzero(cost==9);assert len(high)==3456
# The cheap-only incidence bound requires at least one nine-conflict point.
row=coo_matrix((-np.ones(len(high)),(np.zeros(len(high)),high)),shape=(1,M+N)).tocsr();A=vstack([A,row],format='csr');upper=np.r_[upper,-4.]
# Any high point loses at least one cheap point over the 220 cheap-source set.
cheapids=np.flatnonzero(cheap);cols=np.r_[cheapids,M+affected];data=np.r_[np.ones(len(cheapids)),-np.ones(len(affected))];row=coo_matrix((data,(np.zeros(len(cols)),cols)),shape=(1,M+N)).tocsr();A=vstack([A,row],format='csr');upper=np.r_[upper,-1.]
# Aggregate exact high/cheap witness cuts by source point.
witness_source=np.array([v for v in affected if np.any((B[high]==v).any(axis=1))],np.int32);Z=len(witness_source);where=np.full(N,-1,np.int32);where[witness_source]=np.arange(Z);A=hstack([A,csr_matrix((A.shape[0],Z))],format='csr');rr=[];cc=[];vv=[]
for k,v in enumerate(witness_source):
 members=np.flatnonzero(cheap&(B==v).any(axis=1));rr.extend([k]*(len(members)+2));cc.extend(members.tolist()+[M+int(v),M+N+k]);vv.extend(np.where(norm[members]==3,4.,3.).tolist()+[-12.,3.])
cut=coo_matrix((vv,(rr,cc)),shape=(Z,M+N+Z)).tocsr();A=vstack([A,cut],format='csr');upper=np.r_[upper,np.zeros(Z)]
linkrows=[];linkcols=[]
for k,v in enumerate(witness_source):
 members=high[(B[high]==v).any(axis=1)];linkrows.extend(members.tolist());linkcols.extend([k]*len(members))
rr=np.repeat(np.arange(len(linkrows)),2);cc=np.c_[linkrows,M+N+np.array(linkcols)].ravel();vv=np.tile([1.,-1.],len(linkrows));link=coo_matrix((vv,(rr,cc)),shape=(len(linkrows),M+N+Z)).tocsr();A=vstack([A,link],format='csr');upper=np.r_[upper,np.zeros(len(linkrows))];c=np.r_[c,np.zeros(Z)]
print('candidates',M,'source edges',ne,'candidate edges',e,'capacity cuts',len(affected),'witness cuts',Z,flush=True);start=time.time();r=milp(c,integrality=np.r_[np.ones(M+N),np.zeros(Z)],bounds=Bounds(0,1),constraints=LinearConstraint(A,-np.inf,upper),options={'time_limit':int(os.environ.get('E7_SEARCH_SECONDS','300')),'presolve':True});print('result',r.fun,'dual',getattr(r,'mip_dual_bound',None),'status',r.message,'seconds',round(time.time()-start,1),flush=True)
if r.x is not None:
 pick=np.flatnonzero(r.x[:M]>.5);remove=np.flatnonzero(r.x[M:M+N]>.5);P=np.r_[np.delete(X,remove,axis=0),Q[pick]].astype(np.int64);n=np.r_[np.full(N-len(remove),8),norm[pick]].astype(np.int64);assert len(P)>N and np.all((P*P).sum(1)==n)
 # Verify every pair with integer arithmetic, including unequal norms.
 maxblocks=[]
 for st in range(0,len(P),128):
  d=P[st:st+128]@P.T;bad=(d>0)&(4*d*d>n[st:st+128,None]*n[None,:]);
  for k in range(len(d)):bad[k,st+k]=False
  assert not np.any(bad),'invalid candidate'
 np.savez_compressed(p/'multinorm_POTENTIAL.npz',points=P,norm=n);print('EXACT',len(P),'points, all pair inequalities checked',flush=True)
