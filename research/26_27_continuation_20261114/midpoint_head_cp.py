"""Joint midpoint, head, quadratic-orbit and Leech-parent exchanges."""
import argparse,itertools,time,json,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from ortools.sat.python import cp_model
from image_search import above

def main():
 ap=argparse.ArgumentParser();ap.add_argument('dim',type=int);ap.add_argument('pools',type=Path,nargs='+');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=1800);ap.add_argument('--source-cuts',type=int,default=24);ap.add_argument('--max-old-deletions',type=int,default=32);ap.add_argument('--head-deletions',type=int,default=2);ap.add_argument('--seed',type=int,default=337);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();r=(Path(__file__).resolve().parent/'data');d=a.dim;basecount=201180 if d==26 else 209862;z=np.load(r/f'PackingStar/25D-31D/new_bounds_20261111/D{d}_{basecount}_coordinates.npz')
 if d==26:O=2*np.r_[z['rational_points'],z['quadratic_A']];OB=2*np.r_[np.zeros_like(z['rational_points']),z['quadratic_B']]
 else:O=2*z['points'];OB=np.zeros_like(O)
 L=np.load(r/'leech.npy');lk={v.tobytes():i for i,v in enumerate(L)};canon=np.array([i if v[np.flatnonzero(v)[0]]>0 else lk[(-v).tobytes()] for i,v in enumerate(L)]);g=np.load(r/f'parentgraph{d}/graph.npz');sf=g['source_for_leech'];T=[]
 for i,j in itertools.combinations(range(d-24),2):
  for s in [-1,1]:
   t=np.zeros(d-24,np.int16);t[[i,j]]=[4,4*s];T.append(t)
 T=np.array(T);nt=len(T);tl={v.tobytes():i for i,v in enumerate(T)};bodies={};sourcecolor={};parents={}
 for h,t,p in zip(g['bodies'],g['tails'],g['parents']):
  if h[np.flatnonzero(h)[0]]<0:continue
  t=t.copy()
  if t[np.flatnonzero(t)[0]]<0:t=-t
  key=h.tobytes();c=tl[t.tobytes()];bodies[key]=[h.copy(),{c,c^1}];sourcecolor[key]=c;parents[key]=np.sort(canon[p])
 for path in a.pools:
  z=np.load(path);A=z['candidate_A'] if 'candidate_A' in z else z['candidates'];B=z['candidate_B'] if 'candidate_B' in z else np.zeros_like(A);den=z['denominators'] if 'denominators' in z else np.full(len(A),int(z['denominator']) if 'denominator' in z else 1)
  for i in np.flatnonzero(~np.any(B,axis=1)):
   if np.any(A[i]%int(den[i])):continue
   v=(A[i]//int(den[i])).astype(np.int16)
   if np.sum(v[:24].astype(np.int64)**2)!=96:continue
   t=v[24:].copy()
   if not np.any(t):continue
   if t[np.flatnonzero(t)[0]]<0:t=-t
   if t.tobytes() not in tl:continue
   h=v[:24].copy()
   if h[np.flatnonzero(h)[0]]<0:h=-h
   key=h.tobytes();c=tl[t.tobytes()]
   if key not in bodies:bodies[key]=[h,set()]
   bodies[key][1].update([c,c^1])
   if 'parent_indices' in z:parents[key]=np.sort(canon[z['parent_indices'][i]])
 unknown=[key for key in bodies if key not in parents];LF=L.astype(np.float32)
 for at in range(0,len(unknown),64):
  keys=unknown[at:at+64];HH=np.array([bodies[k][0] for k in keys],np.float32);D=HH@LF.T
  for i,key in enumerate(keys):
   pp=np.flatnonzero(D[i]==48)
   if len(pp)!=2:continue
   assert np.array_equal(L[pp].sum(0),bodies[key][0]);parents[key]=np.sort(canon[pp])
 vals=[];options=[]
 for key,(h,co) in bodies.items():
  if key not in parents:continue
  group=[]
  for c in sorted(co):group.append(len(vals));vals.append((h,c,parents[key],sourcecolor.get(key,-1)==c))
  options.append(group)
 H=np.array([v[0] for v in vals],np.int16);colors=np.array([v[1] for v in vals],np.int16);P=np.array([v[2] for v in vals],np.int32);old=np.array([v[3] for v in vals]);origids=np.flatnonzero(old);assert len(origids)==len(sourcecolor);core=np.sum(O[:,:24].astype(np.int64)**2+3*OB[:,:24].astype(np.int64)**2,axis=1);headids=np.flatnonzero(core==64);headgroups={}
 for i in headids:
  hh=O[i,:24].copy();tt=O[i,24:].copy()
  if hh[np.flatnonzero(hh)[0]]<0:hh=-hh
  if tt[np.flatnonzero(tt)[0]]<0:tt=-tt
  headgroups.setdefault(np.r_[hh,tt].tobytes(),[]).append(int(i))
 headrows=list(headgroups.values());assert all(len(v)==4 for v in headrows);HA=np.array([O[v[0],:24] for v in headrows],np.int16);HT=np.array([O[v[0],24:] for v in headrows],np.int16);quadgroups={}
 for i in np.flatnonzero(np.any(OB,axis=1)):
  bb=OB[i,:24].copy();tt=O[i,24:].copy()
  if bb[np.flatnonzero(bb)[0]]<0:bb=-bb
  if tt[np.flatnonzero(tt)[0]]<0:tt=-tt
  quadgroups.setdefault(np.r_[bb,tt].tobytes(),[]).append(int(i))
 quadrows=list(quadgroups.values());assert all(len(v)==4 for v in quadrows);QB=np.array([OB[v[0],:24] for v in quadrows],np.int16).reshape(-1,24);QT=np.array([O[v[0],24:] for v in quadrows],np.int16).reshape(-1,d-24);QP=np.array([canon[lk[v.tobytes()]] for v in QB],np.int32);pure=np.flatnonzero(core==0);vertices=np.unique(np.r_[P.ravel(),QP]);where={int(v):i for i,v in enumerate(vertices)};incident=[[] for _ in vertices]
 for i,ps in enumerate(P):
  for p in ps:incident[where[int(p)]].append(i)
 model=cp_model.CpModel();x=[model.new_bool_var(f'e{i}') for i in range(len(H))];heads=[model.new_bool_var(f'h{i}') for i in headrows];quads=[model.new_bool_var(f'q{i}') for i in quadrows];used=[model.new_bool_var(f'r{i}') for i in vertices];double=[model.new_bool_var(f'd{i}') for i in vertices];qowner={int(p):j for j,p in enumerate(QP)};origdeg=np.bincount(np.ravel([[where[int(p)] for p in P[i]] for i in origids]),minlength=len(vertices));hintvars=[];hintvals=[]
 def hint(v,b):hintvars.append(v);hintvals.append(int(b))
 for j,(p,es) in enumerate(zip(vertices,incident)):
  q=quads[qowner[int(p)]] if int(p) in qowner else 0;deg=sum(x[i] for i in es);model.add(deg+2*q<=2*used[j]);model.add(used[j]<=deg+q);model.add(double[j]==deg+2*q-used[j]);hint(used[j],sf[p]<0);hint(double[j],int(origdeg[j])+2*int(int(p) in qowner)-int(sf[p]<0))
  for i in es:model.add(x[i]<=used[j])
 for group in options:
  if len(group)>1:model.add_at_most_one([x[i] for i in group])
 for i,v in enumerate(x):hint(v,bool(old[i]))
 for v in heads+quads:hint(v,1)
 if a.max_old_deletions:model.add(sum(1-x[int(i)] for i in origids)<=a.max_old_deletions)
 if a.head_deletions>=0:model.add(sum(heads)>=len(heads)-a.head_deletions)
 HF=H.astype(np.float32);headbad=np.abs(HF@HA.astype(np.float32).T)+np.abs(T[colors].astype(np.float32)@HT.astype(np.float32).T)>64
 for i,j in np.argwhere(headbad):model.add_at_most_one(x[int(i)],heads[int(j)])
 if len(quads):
  qb=above(np.abs(T[colors].astype(np.float32)@QT.astype(np.float32).T),np.abs(HF@QB.astype(np.float32).T))
  for i,j in np.argwhere(qb):model.add_at_most_one(x[int(i)],quads[int(j)])
 score=(basecount-len(L)-len(pure))//2;assert int(np.sum(origdeg+2*np.isin(vertices,QP)-(sf[vertices]<0)))+2*len(heads)==score;model.maximize(sum(double)+2*sum(heads));cost=np.abs(T.astype(np.int64)@T.T);cuts=set()
 for at in range(0,len(H),128):
  D=np.abs(HF[at:at+128]@HF[origids].T);bad=D+cost[colors[at:at+128,None],colors[origids][None,:]]>64
  for j in range(len(D)):
   ii=at+j;ids=np.flatnonzero(bad[j])
   if len(ids)>a.source_cuts:ids=ids[np.argsort(-D[j,ids],kind='stable')[:a.source_cuts]]
   for k in ids:
    jj=int(origids[k]);pair=(min(ii,jj),max(ii,jj))
    if pair[0]!=pair[1] and pair not in cuts:cuts.add(pair);model.add_at_most_one(x[pair[0]],x[pair[1]])
 for v,b in zip(hintvars,hintvals):model.add_hint(v,b)
 print('MODEL',json.dumps(dict(dimension=d,bodies=len(H),parent_lines=len(vertices),heads=len(heads),quads=len(quads),original_bodies=len(origids),source_cuts=len(cuts),target_score=score+1,seconds=time.monotonic()-start)),flush=True);rounds=[]
 class Improve(cp_model.CpSolverSolutionCallback):
  def __init__(self):super().__init__();self.best=-1
  def on_solution_callback(self):
   if self.objective_value>self.best:self.best=self.objective_value;print('MASTER_INCUMBENT',int(self.best),'target',score+1,'bound',self.best_objective_bound,'seconds',self.wall_time,flush=True)
   if self.objective_value>=score+1:self.stop_search()
 for rd in range(100):
  remain=a.seconds-(time.monotonic()-start)
  if remain<=0:break
  solver=cp_model.CpSolver();solver.parameters.num_search_workers=1;solver.parameters.random_seed=a.seed+rd;solver.parameters.max_time_in_seconds=remain;solver.parameters.symmetry_level=0;solver.parameters.linearization_level=2;st=solver.solve(model,Improve());rec=dict(round=rd,status=solver.status_name(st),cuts=len(cuts),seconds=time.monotonic()-start);rounds.append(rec)
  if st not in (cp_model.FEASIBLE,cp_model.OPTIMAL):print('ROUND',json.dumps(rec),flush=True);break
  rec.update(score=float(solver.objective_value),bound=float(solver.best_objective_bound))
  if solver.objective_value<score+1:print('ROUND_NO_POSITIVE',json.dumps(rec),flush=True);break
  chosen=np.flatnonzero([solver.boolean_value(v) for v in x]);newcuts=[]
  for at in range(0,len(chosen),128):
   ids=chosen[at:at+128];D=np.abs(HF[ids]@HF[chosen].T);bad=D+cost[colors[ids,None],colors[chosen][None,:]]>64;ii,jj=np.where(bad);aa=ids[ii];bb=chosen[jj]
   for u,v in zip(aa[aa<bb],bb[aa<bb]):
    pair=(int(u),int(v))
    if pair not in cuts:newcuts.append(pair);cuts.add(pair);model.add_at_most_one(x[pair[0]],x[pair[1]])
  rec.update(selected_bodies=len(chosen),new_conflicts=len(newcuts),heads=int(sum(solver.boolean_value(v) for v in heads)),quads=int(sum(solver.boolean_value(v) for v in quads)));print('ROUND',json.dumps(rec),flush=True)
  if not newcuts:
   rp=vertices[np.flatnonzero([solver.boolean_value(v) for v in used])];anti=np.array([lk[(-v).tobytes()] for v in L]);keep=np.ones(len(L),bool);keep[rp]=False;keep[anti[rp]]=False;LC=np.c_[2*L[keep],np.zeros((int(keep.sum()),d-24),np.int16)];fixed=list(pure)
   for j,v in enumerate(heads):
    if solver.boolean_value(v):fixed.extend(headrows[j])
   for j,v in enumerate(quads):
    if solver.boolean_value(v):fixed.extend(quadrows[j])
   C=np.array([np.r_[s*H[j],t*T[colors[j]]] for j in chosen for s,t in itertools.product([-1,1],repeat=2)],np.int16);RA=np.r_[LC,O[fixed],C];RB=np.r_[np.zeros_like(LC),OB[fixed],np.zeros_like(C)];assert len(RA)>basecount;assert np.all(np.sum(RA.astype(np.int64)**2+3*RB.astype(np.int64)**2,axis=1)==128) and not np.any(np.sum(RA.astype(np.int64)*RB.astype(np.int64),axis=1));sourcekeys={np.r_[v,b].tobytes() for v,b in zip(O,OB)};keys=[np.r_[v,b].tobytes() for v,b in zip(RA,RB)];assert len(set(keys))==len(RA);added=np.flatnonzero([k not in sourcekeys for k in keys]);newkeys=set(keys);removed=np.flatnonzero([np.r_[v,b].tobytes() not in newkeys for v,b in zip(O,OB)]);A64=RA.astype(np.int64);B64=RB.astype(np.int64)
   for j in added:
    fail=above(A64@A64[j]+3*B64@B64[j],A64@B64[j]+B64@A64[j]);fail[j]=False;assert not fail.any(),('changed geometry',j,np.flatnonzero(fail)[:8])
   np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added_A=RA[added],added_B=RB[added],removed_indices=removed);rec.update(count=len(RA),gain=len(RA)-basecount,added=len(added),removed=len(removed));print('POSITIVE_EXACT_CHANGED_PAIRS',json.dumps(rec),flush=True);break
  model.clear_hints()
  for v in hintvars:model.add_hint(v,solver.value(v))
 (a.out/'report.json').write_text(json.dumps(dict(dimension=d,rounds=rounds,seconds=time.monotonic()-start),indent=2)+'\n');print('FINISHED',json.dumps(rounds[-1] if rounds else {}),flush=True)
if __name__=='__main__':main()
