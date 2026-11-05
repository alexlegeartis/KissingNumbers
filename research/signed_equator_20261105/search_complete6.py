import numpy as np,itertools,collections,sys
from pathlib import Path
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
z=np.load(Path(__file__).with_name('source.npz'));X0=np.c_[z['rational'],z['radical']].astype(np.int64);X=5*X0;A=X[:,:13];R=X[:,13:];A0=X0[:,:13];R0=X0[:,13:];words=np.arange(4096,dtype=np.uint16);bits=np.arange(12,dtype=np.uint16);head=(1-2*((words[:,None]>>bits)&1).astype(np.int16)).astype(np.int16);H=np.zeros((8192,13),np.int16);H[:4096,:12]=head;H[4096:,:12]=head;H[:4096,12]=2;H[4096:,12]=-2
E=[];owners={};groups={};owner_inc=[]
for loc in itertools.combinations(range(12),4):
 for signs in itertools.product((-2,2),repeat=4):
  v=np.zeros(13,np.int16);v[list(loc)]=signs;idx=len(E);E.append(v)
  for triple in itertools.combinations(loc,3):
   key=(triple,tuple(int(v[j]) for j in triple));extra=next(j for j in loc if j not in triple);o=(key,extra)
   if o not in owners:owners[o]=len(owners);groups.setdefault(key,[]).append(owners[o])
   owner_inc.append((idx,owners[o]))
E=np.array(E,np.int16);assert len(E)==7920 and len(owners)==15840 and len(groups)==1760;T=np.load(Path(__file__).with_name('extras_complete_den5.npy')).astype(np.int64);N=len(H)+len(E)+len(T);D=len(X);P=len(owners);VA=np.r_[5*H,5*E,T[:,:13]].astype(np.int64);VB=np.r_[np.zeros((len(H)+len(E),13),np.int64),T[:,13:]]
p=VA@A.T+3*VB@R.T-200;q=VA@R.T+VB@A.T;bad=((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q));ii,jj=np.where(bad);ne=len(ii);print('candidate costs',collections.Counter(bad.sum(1)).most_common(12),'source edges',ne,flush=True)
rows=[np.repeat(np.arange(ne,dtype=np.int32),2)];cols=[np.column_stack((ii,N+jj)).ravel().astype(np.int32)];vals=[np.tile(np.array([1.,-1.]),ne)];ubs=[np.zeros(ne)];row=ne
# Every pair of head words at Hamming distance at most three lies in a two-center radius-one ball.
cliques=[]
for base in (0,4096):
 for bit in range(12):
  for w in range(4096):
   u=w^(1<<bit)
   if w>u:continue
   vs={w,u}
   for k in range(12):vs.add(w^(1<<k));vs.add(u^(1<<k))
   assert len(vs)==24;cliques.append([base+x for x in sorted(vs)])
cliques=np.array(cliques,np.int32);nc=len(cliques);rows.append(np.repeat(np.arange(row,row+nc,dtype=np.int32),24));cols.append(cliques.ravel());vals.append(np.ones(nc*24));ubs.append(np.ones(nc));row+=nc
# Owners identify the fourth support coordinate of each signed triple.
inc=np.array(owner_inc,np.int32);ni=len(inc);rows.append(np.repeat(np.arange(row,row+ni,dtype=np.int32),2));cols.append(np.column_stack((len(H)+inc[:,0],N+D+inc[:,1])).ravel());vals.append(np.tile([1.,-1.],ni));ubs.append(np.zeros(ni));row+=ni
rr=[];cc=[]
for ids in groups.values():rr.extend([row]*len(ids));cc.extend(N+D+np.array(ids));row+=1
rows.append(np.array(rr,np.int32));cols.append(np.array(cc,np.int32));vals.append(np.ones(len(rr)));ubs.append(np.ones(len(groups)))
# Extra low-blocker shell vectors use direct exact Q(sqrt3) pair comparisons.
t0=len(H)+len(E);p=T[:,:13]@VA.T+3*T[:,13:]@VB.T-200;q=T[:,:13]@VB.T+T[:,13:]@VA.T;bad=((q==0)&(p>0))|((q>0)&((p>=0)|((p<0)&(3*q*q>p*p))))|((q<0)&(p>0)&(p*p>3*q*q));i,j=np.where(bad);i=i+t0;use=(j<t0)|(j>i);i=i[use];j=j[use];nt=len(i);rows.append(np.repeat(np.arange(row,row+nt,dtype=np.int32),2));cols.append(np.column_stack((i,j)).ravel());vals.append(np.ones(nt*2));ubs.append(np.ones(nt));row+=nt;print('extra',len(T),'extra edges',nt,flush=True)
# Published A(12,4)=144; signed-triple counting gives equator <=816.
for lo,hi,ub in [(0,4096,144),(4096,8192,144),(8192,8192+len(E),816)]:
 rows.append(np.full(hi-lo,row,np.int32));cols.append(np.arange(lo,hi,dtype=np.int32));vals.append(np.ones(hi-lo));ubs.append(np.array([ub],float));row+=1
# Leave-parity packing bound: ≤136 for each signed coordinate, ≤20 for each signed pair.
for i in range(12):
 for sign in (-2,2):
  ids=np.flatnonzero(E[:,i]==sign)+len(H);rows.append(np.full(len(ids),row,np.int32));cols.append(ids);vals.append(np.ones(len(ids)));ubs.append(np.array([136.0]));row+=1
for i in range(12):
 for j in range(i+1,12):
  for si in (-2,2):
   for sj in (-2,2):
    ids=np.flatnonzero((E[:,i]==si)&(E[:,j]==sj))+len(H);rows.append(np.full(len(ids),row,np.int32));cols.append(ids);vals.append(np.ones(len(ids)));ubs.append(np.array([20.0]));row+=1
# Exact source-saturation cuts: E=816, each head sign=144, and source outsiders/axes <=50.
z=(R0==0).all(1);oldE=z&(A0[:,12]==0)&(np.count_nonzero(A0[:,:12],axis=1)==4)&np.isin(A0[:,:12],[-2,0,2]).all(1);oldH=z&(np.abs(A0[:,:12])==1).all(1)&(np.abs(A0[:,12])==2);oldP=oldH&(A0[:,12]==2);oldM=oldH&(A0[:,12]==-2);oldO=~(oldE|oldH);assert (oldE.sum(),oldP.sum(),oldM.sum(),oldO.sum())==(816,144,144,50)
for newids,oldids in [(np.arange(8192,8192+len(E)),np.flatnonzero(oldE)),(np.arange(4096),np.flatnonzero(oldP)),(np.arange(4096,8192),np.flatnonzero(oldM)),(np.arange(8192+len(E),8192+len(E)+26),np.flatnonzero(oldO))]:
 rows.append(np.full(len(newids)+len(oldids),row,np.int32));cols.append(np.r_[newids,N+oldids].astype(np.int32));vals.append(np.r_[np.ones(len(newids)),-np.ones(len(oldids))]);ubs.append(np.array([0.]));row+=1
# The target is a strict improvement, and pure-family saturation makes extras mandatory.
rows.append(np.full(N+D,row,np.int32));cols.append(np.arange(N+D,dtype=np.int32));vals.append(np.r_[-np.ones(N),np.ones(D)]);ubs.append(np.array([-1.]));row+=1
ids=np.arange(len(H)+len(E)+26,N,dtype=np.int32);rows.append(np.full(len(ids),row,np.int32));cols.append(ids);vals.append(-np.ones(len(ids)));ubs.append(np.array([-1.]));row+=1
M=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(row,N+D+P)).tocsr();obj=np.r_[-np.ones(N),np.ones(D),np.zeros(P)];integrality=np.r_[np.ones(N+D),np.zeros(P)];print('model vars',len(obj),'rows',row,'nnz',M.nnz,'head',len(H),'equator',len(E),flush=True)
if '--build-only' in sys.argv: raise SystemExit(0)
seconds=float(sys.argv[sys.argv.index('--seconds')+1]) if '--seconds' in sys.argv else 3600
res=milp(obj,integrality=integrality,bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,np.concatenate(ubs)),options={'time_limit':seconds});print('status',res.message,'gain',None if res.fun is None else -res.fun,'dual',res.mip_dual_bound,flush=True)
if res.fun is not None and res.fun<-.5:
 ids=np.flatnonzero(res.x[:N]>.5);rem=np.flatnonzero(res.x[N:N+D]>.5);Y=np.r_[np.delete(X,rem,axis=0),np.c_[VA[ids],VB[ids]]];np.savez_compressed('candidate_complete6.npz',rational=Y[:,:13],radical=Y[:,13:],den=5);print('POTENTIAL',len(Y),flush=True)
