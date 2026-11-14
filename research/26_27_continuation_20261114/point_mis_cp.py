"""Optimize individual points in a mixed-field exchange graph, with an incumbent."""
import argparse,time,json,sys
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from ortools.sat.python import cp_model
from radical_radii import bsign,qsign

def main():
 ap=argparse.ArgumentParser();ap.add_argument('graph',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--max-cost',type=int,default=4);ap.add_argument('--seconds',type=float,default=1200);ap.add_argument('--seed',type=int,default=353);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();z=np.load(a.graph/'nodes.npz');A=z['A'];B=z['B'];rad=z['radicands'];den=z['denominators'];commonA=z['common_A'];commonB=z['common_B'];baseline=int(z['baseline_count']);origids=z['old_indices']
 with open(a.graph/'graph.bin','rb') as f:
  n,m=map(int,np.fromfile(f,'<u4',2));ne=int(np.fromfile(f,'<u8',1)[0]);initial=np.fromfile(f,'<u4',m);ptr=np.fromfile(f,'<u8',n+1);idx=np.fromfile(f,'<u4',ne)
 assert np.array_equal(initial,np.arange(m));G=csr_matrix((np.ones(ne,np.int8),idx.astype(np.int32),ptr.astype(np.int64)),shape=(n,n));cost=np.asarray(G[m:,:m].getnnz(axis=1)).ravel();new=m+np.flatnonzero(cost<=a.max_cost);cross=G[new,:m];old=np.unique(cross.indices);untouched=np.setdiff1d(np.arange(m),old);ids=np.r_[old,new];GG=G[ids][:,ids].tocsr();mm=len(old);nn=len(ids);assert GG[:mm,:mm].nnz==0;model=cp_model.CpModel();x=[model.new_bool_var(f'p{i}') for i in range(nn)];ii,jj=GG.nonzero();keep=ii<jj;ii=ii[keep];jj=jj[keep]
 for i,j in zip(ii,jj):model.add_at_most_one(x[int(i)],x[int(j)])
 cliques=0
 for v in range(mm):
  neighbors=GG.indices[GG.indptr[v]:GG.indptr[v+1]]
  if len(neighbors)<2:continue
  local=GG[neighbors][:,neighbors].toarray().astype(bool);left=np.ones(len(neighbors),bool)
  while left.any():
   allowed=left.copy();chosen=[]
   while allowed.any():
    opts=np.flatnonzero(allowed);w=int(opts[np.argmax(local[opts][:,allowed].sum(1))]);chosen.append(w);left[w]=False;allowed[w]=False;allowed&=local[w]
   if len(chosen)>1:model.add_at_most_one([x[v]]+[x[int(neighbors[j])] for j in chosen]);cliques+=1
 model.add(sum(x)>=mm);model.maximize(sum(x))
 for i,v in enumerate(x):model.add_hint(v,i<mm)
 class Improve(cp_model.CpSolverSolutionCallback):
  def __init__(self):super().__init__();self.best=-1
  def on_solution_callback(self):
   if self.objective_value>self.best:self.best=self.objective_value;print('INCUMBENT',int(self.best),'target',mm+1,'seconds',self.wall_time,flush=True)
   if self.objective_value>=mm+1:self.stop_search()
 solver=cp_model.CpSolver();solver.parameters.num_search_workers=1;solver.parameters.max_time_in_seconds=a.seconds;solver.parameters.random_seed=a.seed;solver.parameters.symmetry_level=0;solver.parameters.linearization_level=2;solver.parameters.log_search_progress=True;print('MODEL',json.dumps(dict(nodes=nn,old=mm,new=len(new),edges=len(ii),cliques=cliques,max_cost=a.max_cost,seconds=time.monotonic()-start)),flush=True);st=solver.solve(model,Improve());report=dict(status=solver.status_name(st),nodes=nn,old=mm,new=len(new),edges=len(ii),max_cost=a.max_cost,seconds=time.monotonic()-start)
 if st in (cp_model.FEASIBLE,cp_model.OPTIMAL):
  selected=np.flatnonzero([solver.boolean_value(v) for v in x]);report.update(selected=len(selected),bound=float(solver.best_objective_bound))
  if len(selected)>mm:
   chosen=ids[selected];fixedA=np.r_[commonA,A[untouched]];fixedB=np.r_[commonB,B[untouched]];RA=np.r_[fixedA,A[chosen]];RB=np.r_[fixedB,B[chosen]];RR=np.r_[np.where(np.any(fixedB,axis=1),3,1),rad[chosen]];DD=np.r_[np.ones(len(fixedA),np.int16),den[chosen]];assert len(RA)==baseline+len(selected)-mm;assert np.all(np.sum(RA.astype(np.int64)**2+RR[:,None]*RB.astype(np.int64)**2,axis=1)==128*DD.astype(np.int64)**2) and not np.any(np.sum(RA.astype(np.int64)*RB.astype(np.int64),axis=1));A64=RA.astype(np.int64);B64=RB.astype(np.int64);addedlocal=np.flatnonzero(chosen>=m)
   for loc in addedlocal:
    j=len(fixedA)+loc;ca=A64[j];cb=B64[j];rr=int(RR[j]);dd=int(DD[j]);p=A64@ca;q=B64@ca;c=A64@cb;e=B64@cb;limit=64*DD.astype(np.int64)*dd;same=RR==rr;bad=np.empty(len(RA),bool);bad[same]=qsign(p[same]-limit[same]+rr*e[same],q[same]+c[same],rr)>0;bad[~same]=bsign(p[~same]-limit[~same],q[~same],c[~same],e[~same],RR[~same],rr)>0;bad[j]=False;assert not bad.any(),('direct changed geometry',j,np.flatnonzero(bad)[:8])
   removed=origids[old[np.flatnonzero([not solver.boolean_value(x[i]) for i in range(mm)])]];np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicands=RR,denominators=DD,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added_A=A[chosen[addedlocal]],added_B=B[chosen[addedlocal]],added_radicands=rad[chosen[addedlocal]],added_denominators=den[chosen[addedlocal]],removed_indices=removed);report.update(count=len(RA),gain=len(selected)-mm,added=len(addedlocal),removed=len(removed));print('POSITIVE_EXACT_CHANGED_PAIRS',json.dumps(report),flush=True)
 (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
