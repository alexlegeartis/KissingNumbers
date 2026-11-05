import numpy as np,time,os
from pathlib import Path
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack,hstack,load_npz,csr_matrix
cache=Path(__file__).resolve().parent/'cache';z=np.load(cache/'candidates14.npz');Q=z['Q'].astype(np.int16);B=z['B'];old=z['old'].astype(np.int16);L=len(Q);O=len(old);M=2*L;N=2*O;counts=(B!=65535).sum(1);th=4
# Each line conflict gives one blocker for each chosen candidate sign.
a,b=np.where(B!=65535);bl=B[a,b].astype(int);dots=(Q[a]*old[bl]).sum(axis=1);assert np.all(abs(dots)>th);source=np.where(dots>0,bl,bl+O);ca=np.r_[a,a+L];cb=np.r_[source,np.where(source<O,source+O,source-O)];ne=len(ca)
# All candidate-to-candidate point conflicts, including independently selected signs.
all_i=[];all_j=[]
for i in range(0,L,256):
 D=Q[i:i+256]@Q.T
 for k in range(len(D)):D[k,:i+k+1]=0
 a1,b1=np.where(D>th);a2,b2=np.where(D< -th);a1+=i;a2+=i
 all_i.extend([a1,a1+L,a2,a2+L]);all_j.extend([b1,b1+L,b2+L,b2])
i=np.concatenate(all_i);j=np.concatenate(all_j);e=len(i);print('point candidates',M,'source',N,'source edges',ne,'candidate edges',e,flush=True)
rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)];cols=np.r_[np.c_[ca,M+cb].ravel(),np.c_[i,j].ravel()];vals=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)];A=coo_matrix((vals,(rows,cols)),shape=(ne+e,M+N)).tocsr();upper=np.r_[np.zeros(ne),np.ones(e)]
# Valid capacity cuts from exact 32-vertex checks, for both signs.
rows=[];cols=[];vals=[];k=0
for oldpt in range(N):
 for cost,cap in ((4,4),(8,8)):
  members=ca[(cb==oldpt)&(counts[ca%L]==cost)]
  if not len(members):continue
  assert len(members)==(32 if cost==4 else 64)
  rows.extend([k]*(len(members)+1));cols.extend(members.tolist()+[M+oldpt]);vals.extend([1.]*len(members)+[-float(cap)]);k+=1
cut=coo_matrix((vals,(rows,cols)),shape=(k,M+N)).tocsr();A=vstack([A,cut],format='csr');upper=np.r_[upper,np.zeros(k)]
# Seek strictly more points directly. Each addition costs one variable;
# each deleted baseline point costs one variable.
c=np.r_[-np.ones(M),np.ones(N)];gainrow=coo_matrix(c.reshape(1,-1)).tocsr();A=vstack([A,gainrow],format='csr');upper=np.r_[upper,-1.]
# Exact verifier also shows at most six 14-conflict candidates can never improve
# the code when all remaining candidates are from the 4/8 classes.
high=np.flatnonzero(counts==14);high=np.r_[high,high+L];assert len(high)==3328
row=coo_matrix((-np.ones(len(high)),(np.zeros(len(high)),high)),shape=(1,M+N)).tocsr();A=vstack([A,row],format='csr');upper=np.r_[upper,-7.]
# Any 14-conflict point incurs a deficit of at least six cheap points when seven or more are selected.
cheap=np.flatnonzero(counts<14);cheap=np.r_[cheap,cheap+L];cols=np.r_[cheap,M+np.arange(N)];vals=np.r_[np.ones(len(cheap)),-np.ones(N)];row=coo_matrix((vals,(np.zeros(len(cols)),cols)),shape=(1,M+N)).tocsr();A=vstack([A,row],format='csr');upper=np.r_[upper,-6.]
# Link each 14-conflict point to its exact deficient old-point witnesses.
# Auxiliary z variables encode whether any selected high point witnesses a
# given old point; they may be continuous because selected high points are binary.
W=load_npz(cache/'witness14.npz').tocoo();assert W.shape==(len(high),N)
witness_source=np.unique(W.col);Z=len(witness_source);where=np.full(N,-1,np.int32);where[witness_source]=np.arange(Z)
A=hstack([A,csr_matrix((A.shape[0],Z))],format='csr');rows=[];cols=[];vals=[]
for j,oldpt in enumerate(witness_source):
 members=ca[(cb==oldpt)&(counts[ca%L]==8)];assert len(members)==64
 rows.extend([j]*(len(members)+2));cols.extend(members.tolist()+[M+int(oldpt),M+N+j]);vals.extend([1.]*len(members)+[-8.,1.])
cut=coo_matrix((vals,(rows,cols)),shape=(Z,M+N+Z)).tocsr();A=vstack([A,cut],format='csr');upper=np.r_[upper,np.zeros(Z)]
rr=np.repeat(np.arange(W.nnz),2);cc=np.c_[high[W.row],M+N+where[W.col]].ravel();vv=np.tile([1.,-1.],W.nnz);link=coo_matrix((vv,(rr,cc)),shape=(W.nnz,M+N+Z)).tocsr();A=vstack([A,link],format='csr');upper=np.r_[upper,np.zeros(W.nnz)]
c=np.r_[c,np.zeros(Z)]
print('matrix',A.shape,'capacity cuts',k,'witness variables',Z,'witness links',W.nnz,flush=True);start=time.time();r=milp(c,integrality=np.r_[np.ones(M+N),np.zeros(Z)],bounds=Bounds(0,1),constraints=LinearConstraint(A,-np.inf,upper),options={'time_limit':int(os.environ.get('E7_SEARCH_SECONDS','300')),'presolve':True});print('result',r.fun,'dual',getattr(r,'mip_dual_bound',None),'status',r.message,'time',round(time.time()-start,1),flush=True)
if r.x is not None:
 pick=np.flatnonzero(r.x[:M]>.5);remove=np.flatnonzero(r.x[M:M+N]>.5);P=np.r_[Q,-Q];X=np.r_[old,-old];Y=np.r_[np.delete(X,remove,axis=0),P[pick]];D=Y@Y.T;np.fill_diagonal(D,-9);assert len(set(map(tuple,Y)))==len(Y) and D.max()<=4 and np.all((Y*Y).sum(1)==8);np.save(cache/'point14_POTENTIAL.npy',Y);print('EXACT',len(Y),'max dot',D.max(),flush=True)
