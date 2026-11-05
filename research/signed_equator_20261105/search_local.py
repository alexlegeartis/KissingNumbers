import numpy as np,itertools,collections,sys
from pathlib import Path
root=Path(__file__).parent
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
z=np.load(root/'source.npz');X0=np.c_[z['rational'],z['radical']].astype(np.int64);X=5*X0;A=X[:,:13];R=X[:,13:];words=np.arange(4096,dtype=np.uint16);bits=np.arange(12,dtype=np.uint16);h=1-2*((words[:,None]>>bits)&1).astype(np.int16);H=np.zeros((8192,26),np.int64);H[:4096,:12]=5*h;H[4096:,:12]=5*h;H[:4096,12]=10;H[4096:,12]=-10
E=[]
for loc in itertools.combinations(range(12),4):
 for signs in itertools.product((-10,10),repeat=4):
  v=np.zeros(26,np.int64);v[list(loc)]=signs;E.append(v)
E=np.array(E);T=np.r_[5*np.load(root/'extras.npy').astype(np.int64),np.load(root/'extras_den5.npy').astype(np.int64)];Q=np.r_[H,E,T];a=Q[:,:13];b=Q[:,13:]
def bad(a,b,c,d):
 p=a@c.T+3*b@d.T-200;q=a@d.T+b@c.T
 return ((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q))
def kind(V):
 a=V[:,:13];b=V[:,13:];zero=(b==0).all(1);e=zero&(a[:,12]==0)&(np.count_nonzero(a[:,:12],axis=1)==4)&np.isin(a[:,:12],[-10,0,10]).all(1);h=zero&(np.abs(a[:,:12])==5).all(1)&(np.abs(a[:,12])==10);axes=zero&(np.count_nonzero(a,axis=1)==1)&(np.abs(a)==20).any(1);out=np.full(len(V),4,np.int8);out[e]=0;out[h&(a[:,12]==10)]=1;out[h&(a[:,12]==-10)]=2;out[axes]=3;return out
K=kind(Q);S=kind(X);S[S==4]=3;assert np.bincount(S,minlength=4).tolist()==[816,144,144,50]
B=bad(a,b,A,R);cost=B.sum(1);seeds=np.flatnonzero((K==4)&(cost<=10));print('pool',len(Q),'nonpure seeds',len(seeds),'costs',collections.Counter(cost[seeds]),flush=True);rng=np.random.default_rng(131115)
trials=int(sys.argv[sys.argv.index('--trials')+1]) if '--trials' in sys.argv else 1000
seconds=float(sys.argv[sys.argv.index('--seconds')+1]) if '--seconds' in sys.argv else 3
for trial in range(trials):
 seed=int(rng.choice(seeds));vict=B[seed].copy();last=-999
 for step in range(100):
  missing=B[:,~vict].sum(1);eligible=np.flatnonzero(missing==0);eligible=eligible[eligible!=seed];eligible=eligible[~bad(a[eligible],b[eligible],a[seed:seed+1],b[seed:seed+1]).ravel()]
  if step%4==0 and len(eligible)<2200:
   ids=np.r_[seed,eligible];G=bad(a[ids],b[ids],a[ids],b[ids]);i,j=np.where(np.triu(G,1));rm=np.flatnonzero(vict);ii,jj=np.where(B[ids][:,rm]);n=len(ids);ne=len(ii);e=len(i);rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)];cols=np.r_[np.column_stack((ii,n+jj)).ravel(),np.column_stack((i,j)).ravel()];vals=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)];cr=[];cc=[];cv=[]
   for k in range(4):
    u=np.flatnonzero(K[ids]==k);v=np.flatnonzero(S[rm]==k);cr.extend([ne+e+k]*(len(u)+len(v)));cc.extend(np.r_[u,n+v]);cv.extend(np.r_[np.ones(len(u)),-np.ones(len(v))])
   rows=np.r_[rows,np.array(cr,int)];cols=np.r_[cols,np.array(cc,int)];vals=np.r_[vals,np.array(cv,float)];M=coo_matrix((vals,(rows,cols)),shape=(ne+e+4,n+len(rm))).tocsr();lo=np.zeros(n+len(rm));lo[0]=1;res=milp(np.r_[-np.ones(n),np.ones(len(rm))],integrality=np.ones(n+len(rm)),bounds=Bounds(lo,np.ones(n+len(rm))),constraints=LinearConstraint(M,-np.inf,np.r_[np.zeros(ne),np.ones(e),np.zeros(4)]),options={'time_limit':seconds});last=0 if res.fun is None else round(-res.fun)
   if last>0:
    pick=ids[res.x[:n]>.5];rem=rm[res.x[n:]>.5];keep=np.ones(len(X),bool);keep[rem]=False;Y=np.r_[X[keep],Q[pick]];np.savez_compressed('candidate_local_den5.npz',rational=Y[:,:13],radical=Y[:,13:],den=5);print('POTENTIAL',len(Y),'trial',trial,flush=True);raise SystemExit
  opts=np.flatnonzero((missing>0)&(missing<=4)&(cost<=20));
  if not len(opts):break
  groups=collections.Counter(tuple(np.flatnonzero(B[t]&~vict)) for t in opts);key=max(groups,key=lambda k:groups[k]/(len(k)+.25)*(0.8+.4*rng.random()));vict[list(key)]=True
  if vict.sum()>160:break
 if trial%10==0:print('trial',trial,'victims',vict.sum(),'score',last,'step',step,flush=True)
