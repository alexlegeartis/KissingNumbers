"""Point-level shared-deletion neighborhoods of an arbitrary image-union graph."""
import sys,collections,numpy as np
from scipy.sparse import coo_matrix,csr_matrix
from scipy.optimize import milp,Bounds,LinearConstraint
label=sys.argv[1];seed0=int(sys.argv[2]);rng=np.random.default_rng(seed0);base=f'search{label}_multiimage_'
o=np.fromfile(base+'origin.bin',np.uint8).astype(bool);E=np.fromfile(base+'edges.bin',np.int32).reshape(-1,2);sid=np.flatnonzero(o);qid=np.flatnonzero(~o);ns=len(sid);nq=len(qid);sm=np.full(len(o),-1,np.int32);sm[sid]=np.arange(ns);qm=np.full(len(o),-1,np.int32);qm[qid]=np.arange(nq)
a,b=E[:,0],E[:,1];cross=o[a]^o[b];ss=np.where(o[a[cross]],a[cross],b[cross]);qq=np.where(o[a[cross]],b[cross],a[cross]);B=coo_matrix((np.ones(len(ss),np.int32),(qm[qq],sm[ss])),shape=(nq,ns)).tocsr();bad=E[~o[a]&~o[b]];qedges=np.c_[qm[bad[:,0]],qm[bad[:,1]]];del E,a,b,ss,qq,cross,bad
adj=coo_matrix((np.ones(2*len(qedges),np.int8),(np.r_[qedges[:,0],qedges[:,1]],np.r_[qedges[:,1],qedges[:,0]])),shape=(nq,nq)).tocsr();cost=np.diff(B.indptr);mincost=int(cost.min());low=np.flatnonzero(cost<=mincost+4);wide=np.flatnonzero(cost<=mincost+20);print(label,'source',ns,'candidates',nq,'candidate edges',len(qedges),'mincost',mincost,'seeds',len(low),len(wide),flush=True)
for trial in range(1000):
 seed=int(rng.choice(low if trial%2==0 else wide));vict=set(map(int,B.indices[B.indptr[seed]:B.indptr[seed+1]]));last=None
 neighbor=adj.indices[adj.indptr[seed]:adj.indptr[seed+1]];compatible=np.ones(nq,bool);compatible[neighbor]=False;compatible[seed]=False
 for step in range(120):
  allow=np.zeros(ns,np.int32);allow[list(vict)]=1;miss=cost-np.asarray(B@allow).ravel();eligible=np.flatnonzero((miss==0)&compatible)
  if step%6==0:
   if len(eligible)>700:eligible=rng.choice(eligible,700,replace=False)
   ids=np.r_[seed,eligible];rm=np.array(sorted(vict),np.int32);nb=len(ids);nr=len(rm);lookup=np.full(nq,-1,np.int32);lookup[ids]=np.arange(nb);mask=(lookup[qedges[:,0]]>=0)&(lookup[qedges[:,1]]>=0);ce=lookup[qedges[mask]];bi=B[ids][:,rm].tocoo();ec=len(ce);bc=len(bi.data)
   rows=np.r_[np.repeat(np.arange(bc),2),np.repeat(np.arange(bc,bc+ec),2)];cols=np.r_[np.column_stack((bi.row,nb+bi.col)).ravel(),ce.ravel()];vals=np.r_[np.tile([1.,-1.],bc),np.ones(2*ec)];M=coo_matrix((vals,(rows,cols)),shape=(bc+ec,nb+nr)).tocsr();lb=np.zeros(nb+nr);lb[0]=1
   result=milp(np.r_[-np.ones(nb),np.ones(nr)],integrality=np.ones(nb+nr),bounds=Bounds(lb,np.ones(nb+nr)),constraints=LinearConstraint(M,-np.inf,np.r_[np.zeros(bc),np.ones(ec)]),options={'time_limit':2,'presolve':False});last=None if result.fun is None else round(-result.fun)
   if last is not None and last>0:
    pick=ids[result.x[:nb]>.5];removed=rm[result.x[nb:]>.5];selected=o.copy();selected[sid[removed]]=0;selected[qid[pick]]=1;selected.astype(np.uint8).tofile(base+'lns_potential.bin');print('POTENTIAL',int(selected.sum()),'gain',last,'trial',trial,'step',step,flush=True);raise SystemExit
  opts=np.flatnonzero((miss>0)&(miss<=40)&compatible)
  if not len(opts):break
  if len(opts)>2500:opts=rng.choice(opts,2500,replace=False)
  groups=collections.Counter()
  for j in opts:
   arr=B.indices[B.indptr[j]:B.indptr[j+1]];extra=tuple(int(x) for x in arr if not allow[x]);groups[extra]+=1
  if not groups:break
  scores=[(count/(len(key)+.25),key) for key,count in groups.items()];top=max(v for v,k in scores);short=[k for v,k in scores if v>=top*(.7+.3*rng.random())];vict.update(short[int(rng.integers(len(short)))])
  if len(vict)>400:break
 if trial%5==0:print('trial',trial,'victims',len(vict),'eligible',len(eligible),'gain',last,'step',step,flush=True)
