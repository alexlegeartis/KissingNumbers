import numpy as np,collections,sys
from pathlib import Path
from scipy.sparse import coo_matrix
from scipy.optimize import milp,Bounds,LinearConstraint
root=Path(__file__).parent;z=np.load(root/'source.npz');X=z['X'].astype(np.int16)*12;g=z['g'].astype(np.int16);z=np.load(root/'pool.npz');v=np.load(root/'sparse_low.npz');Q=np.r_[z['Q'],v['Q']];B=np.r_[z['B'],v['blockers']];_,idx=np.unique(Q,axis=0,return_index=True);Q=Q[idx];B=B[idx]; B=np.pad(B,((0,0),(0,2)),constant_values=65535);n0=len(Q)
z=np.load(root/'radical50_low.npz');R=z['Q'];RB=z['blockers'];assert len(R)==1024
Q=np.r_[Q,R];B=np.r_[B,RB];israd=np.arange(len(Q))>=n0
cost=(B!=65535).sum(1);rng=np.random.default_rng(150739)
def category(V):
 a=V[:,:14];t=V[:,14];out=np.full(len(V),6,np.int8)
 e=(t==0)&(np.count_nonzero(a,axis=1)==4)&np.isin(a,[-24,0,24]).all(1);h=(np.abs(a)==12).all(1)&(np.abs(t)==12);l=(np.abs(t)==24);axes=(t==0)&(np.count_nonzero(a,axis=1)==1)&(np.abs(a)==48).any(1)
 out[e]=0;out[h&(t==12)]=1;out[h&(t==-12)]=2;out[l&(t==24)]=3;out[l&(t==-24)]=4;out[axes]=5
 return out
kind=category(Q);source_kind=category(X);assert np.bincount(source_kind,minlength=7).tolist()==[1456,512,512,28,28,28,0];assert (kind==6).sum()==6144
trials=int(sys.argv[sys.argv.index('--trials')+1]) if '--trials' in sys.argv else 1000
seconds=float(sys.argv[sys.argv.index('--seconds')+1]) if '--seconds' in sys.argv else 3
radseeds=np.flatnonzero(israd);qseeds=np.flatnonzero((kind==6)&(~israd)&(cost<=8));seeds=np.r_[radseeds,qseeds];low_seeds=seeds;assert len(radseeds)==1024 and len(qseeds)>0
print('pool',len(Q),'nonpure seeds',len(seeds),'low-blocker seeds',len(low_seeds),'trials',trials,'seconds',seconds,flush=True)
for trial in range(trials):
 seed=int(rng.choice(low_seeds if trial%2==0 else seeds));vict=set(map(int,B[seed][B[seed]!=65535]));last=-999
 for step in range(90):
  allow=np.zeros(65536,bool);allow[65535]=1;allow[list(vict)]=1;missing=np.count_nonzero(~allow[B],axis=1);eligible=np.flatnonzero(missing==0);eligible=eligible[eligible!=seed];V=Q[eligible];d=(V*g)@Q[seed];eligible=eligible[d<=np.where(israd[seed],np.where(israd[eligible],25,169),np.where(israd[eligible],169,1152))]
  if trial<2 and step%10==0:print('progress',trial,'step',step,'victims',len(vict),'eligible',len(eligible),flush=True)
  if step%8==0:
   if len(eligible)>600:
    sp=eligible[kind[eligible]==6];pu=eligible[kind[eligible]!=6];sp=rng.choice(sp,min(len(sp),350),False);pu=rng.choice(pu,min(len(pu),600-len(sp)),False);eligible=np.r_[sp,pu]
   ids=np.r_[seed,eligible];V=Q[ids];G=(V*g)@V.T;threshold=np.where(israd[ids,None],np.where(israd[ids][None,:],25,169),np.where(israd[ids][None,:],169,1152));i,j=np.where(np.triu(G>threshold,1));rm=np.array(sorted(vict));i1,j1=np.where(np.isin(B[ids],rm));lookup=np.full(65536,-1,np.int32);lookup[rm]=np.arange(len(rm));jj=lookup[B[ids][i1,j1]];n=len(ids);ne=len(i1);e=len(i);rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)];cols=np.r_[np.column_stack((i1,n+jj)).ravel(),np.column_stack((i,j)).ravel()];vals=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)];cutrows=[];cutcols=[];cutvals=[]
   for c in range(6):
    u=np.flatnonzero(kind[ids]==c);v=np.flatnonzero(source_kind[rm]==c);cutrows.extend([ne+e+c]*(len(u)+len(v)));cutcols.extend(np.r_[u,n+v]);cutvals.extend(np.r_[np.ones(len(u)),-np.ones(len(v))])
   rows=np.r_[rows,np.array(cutrows,int)];cols=np.r_[cols,np.array(cutcols,int)];vals=np.r_[vals,np.array(cutvals,float)];M=coo_matrix((vals,(rows,cols)),shape=(ne+e+6,n+len(rm))).tocsr();lb=np.zeros(n+len(rm));lb[0]=1;res=milp(np.r_[-np.ones(n),np.ones(len(rm))],integrality=np.ones(n+len(rm)),bounds=Bounds(lb,np.ones(n+len(rm))),constraints=LinearConstraint(M,-np.inf,np.r_[np.zeros(ne),np.ones(e),np.zeros(6)]),options={'time_limit':seconds,'presolve':False});last=0 if res.fun is None else round(-res.fun)
   if last>0:
    pick=ids[res.x[:n]>.5];rem=rm[res.x[n:]>.5];keep=np.ones(len(X),bool);keep[rem]=False;out=np.r_[X[keep],Q[pick[~israd[pick]]]];radical=Q[pick[israd[pick]]];np.savez_compressed('search15_hybrid_lns_potential.npz',normal=out,radical=radical,g=g);print('POTENTIAL',len(out)+len(radical),'trial',trial,flush=True);raise SystemExit
  opts=np.flatnonzero((missing>0)&(missing<=6));
  if not len(opts):break
  extras=np.where(~allow[B[opts]],B[opts],65535);extras.sort(axis=1);groups,inverse=np.unique(extras,axis=0,return_inverse=True);counts=np.bincount(inverse,weights=np.where(kind[opts]==6,8.0,1.0));sizes=np.count_nonzero(groups!=65535,axis=1);score=counts/(sizes+.25);top=np.flatnonzero(score>=score.max()*(.65+.35*rng.random()));choice=int(rng.choice(top));vict.update(map(int,groups[choice][groups[choice]!=65535]));
  if len(vict)>200:break
 if trial%10==0:print('trial',trial,'victims',len(vict),'score',last,'step',step,flush=True)
