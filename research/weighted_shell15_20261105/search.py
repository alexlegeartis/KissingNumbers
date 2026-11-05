import numpy as np,collections,sys
from pathlib import Path
from scipy.sparse import coo_matrix
from scipy.optimize import milp,Bounds,LinearConstraint
root=Path(__file__).parent;z=np.load(root/'source.npz');X=z['X'].astype(np.int16)*12;g=z['g'].astype(np.int16);z=np.load(root/'pool.npz');Q=z['Q'];B=z['B'];cost=(B!=65535).sum(1);rng=np.random.default_rng(159913)
def category(V):
 a=V[:,:14];t=V[:,14];out=np.full(len(V),6,np.int8)
 e=(t==0)&(np.count_nonzero(a,axis=1)==4)&np.isin(a,[-24,0,24]).all(1);h=(np.abs(a)==12).all(1)&(np.abs(t)==12);l=(np.count_nonzero(a,axis=1)==2)&np.isin(a,[-24,0,24]).all(1)&(np.abs(t)==24);axes=(t==0)&(np.count_nonzero(a,axis=1)==1)&(np.abs(a)==48).any(1)
 out[e]=0;out[h&(t==12)]=1;out[h&(t==-12)]=2;out[l&(t==24)]=3;out[l&(t==-24)]=4;out[axes]=5
 return out
kind=category(Q);source_kind=category(X);assert np.bincount(source_kind,minlength=7).tolist()==[1456,512,512,28,28,28,0];assert (kind==6).sum()==246
trials=int(sys.argv[sys.argv.index('--trials')+1]) if '--trials' in sys.argv else 1000
seconds=float(sys.argv[sys.argv.index('--seconds')+1]) if '--seconds' in sys.argv else 3
print('pool',len(Q),'nonpure seeds',int((kind==6).sum()),'trials',trials,'seconds',seconds,flush=True)
for trial in range(trials):
 seed=int(rng.choice(np.flatnonzero((cost>=2)&(cost<=10)&(kind==6))));vict=set(map(int,B[seed][B[seed]!=65535]));last=-999
 for step in range(90):
  allow=np.zeros(65536,bool);allow[65535]=1;allow[list(vict)]=1;missing=np.count_nonzero(~allow[B],axis=1);eligible=np.flatnonzero(missing==0);eligible=eligible[eligible!=seed];V=Q[eligible];d=(V*g)@Q[seed];eligible=eligible[d<=1152]
  if step%4==0 and len(eligible)<2500:
   ids=np.r_[seed,eligible];V=Q[ids];G=(V*g)@V.T;i,j=np.where(np.triu(G>1152,1));rm=np.array(sorted(vict));i1,j1=np.where(np.isin(B[ids],rm));lookup=np.full(65536,-1,np.int32);lookup[rm]=np.arange(len(rm));jj=lookup[B[ids][i1,j1]];n=len(ids);ne=len(i1);e=len(i);rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)];cols=np.r_[np.column_stack((i1,n+jj)).ravel(),np.column_stack((i,j)).ravel()];vals=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)];cutrows=[];cutcols=[];cutvals=[]
   for c in range(6):
    u=np.flatnonzero(kind[ids]==c);v=np.flatnonzero(source_kind[rm]==c);cutrows.extend([ne+e+c]*(len(u)+len(v)));cutcols.extend(np.r_[u,n+v]);cutvals.extend(np.r_[np.ones(len(u)),-np.ones(len(v))])
   rows=np.r_[rows,np.array(cutrows,int)];cols=np.r_[cols,np.array(cutcols,int)];vals=np.r_[vals,np.array(cutvals,float)];M=coo_matrix((vals,(rows,cols)),shape=(ne+e+6,n+len(rm))).tocsr();lb=np.zeros(n+len(rm));lb[0]=1;res=milp(np.r_[-np.ones(n),np.ones(len(rm))],integrality=np.ones(n+len(rm)),bounds=Bounds(lb,np.ones(n+len(rm))),constraints=LinearConstraint(M,-np.inf,np.r_[np.zeros(ne),np.ones(e),np.zeros(6)]),options={'time_limit':seconds});last=0 if res.fun is None else round(-res.fun)
   if last>0:
    pick=ids[res.x[:n]>.5];rem=rm[res.x[n:]>.5];keep=np.ones(len(X),bool);keep[rem]=False;out=np.r_[X[keep],Q[pick]];np.savez_compressed('candidate.npz',X=out,g=g,den=12);print('POTENTIAL',len(out),'trial',trial,flush=True);raise SystemExit
  opts=np.flatnonzero((missing>0)&(missing<=4));
  if not len(opts):break
  extras=np.where(~allow[B[opts]],B[opts],65535);extras.sort(axis=1);groups,counts=np.unique(extras,axis=0,return_counts=True);sizes=np.count_nonzero(groups!=65535,axis=1);score=counts/(sizes+.25);top=np.flatnonzero(score>=score.max()*(.65+.35*rng.random()));choice=int(rng.choice(top));vict.update(map(int,groups[choice][groups[choice]!=65535]));
  if len(vict)>150:break
 if trial%10==0:print('trial',trial,'victims',len(vict),'score',last,'step',step,flush=True)
