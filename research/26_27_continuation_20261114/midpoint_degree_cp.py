"""Optimize midpoint edges and the released Leech parents directly, allowing paths."""
import argparse,itertools,time,json,sys
from pathlib import Path
import numpy as np
from ortools.sat.python import cp_model
from image_search import above

def main():
 ap=argparse.ArgumentParser();ap.add_argument('pack',type=Path);ap.add_argument('priced',type=Path,nargs='+');ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=1200);ap.add_argument('--seed',type=int,default=317);a=ap.parse_args();a.out.mkdir(exist_ok=False);start=time.monotonic();z=np.load(a.pack/'data.npz');d=int(z['dimension']);T=z['tail_palette'];nt=len(T);L=z['leech'];F=z['fixed_A'];FB=z['fixed_B'];blocked=z['blocked'];flat=z['state'][z['initial']];ci=flat//nt;phase=flat%nt;OH=z['H'].reshape(-1,4,24)[ci];OC=z['offset'].reshape(-1,4)[ci]^phase[:,None];OP=z['parents'][ci];pool={}
 def put(h,c,parents,original=False):
  h=h.copy()
  if h[np.flatnonzero(h)[0]]<0:h=-h
  key=np.r_[h,np.int16(c)].tobytes();dots=np.abs(L[parents].astype(np.int64)@h.astype(np.int64));pp=np.sort(parents[dots==48]);assert len(pp)==2
  if key in pool:
   assert np.array_equal(pool[key][2],pp)
   if original:pool[key][3]=True
  else:pool[key]=[h,int(c),pp,bool(original)]
 for hs,cs,ps in zip(OH,OC,OP):
  for h,c in zip(hs,cs):put(h,c,ps,True)
 original_count=len(pool);assert original_count==4*len(ci)
 for path in a.priced:
  zz=np.load(path)
  for hs,off,ph,ps in zip(zz['H'],zz['offset'],zz['phase'],zz['parents']):
   for h,c in zip(hs,off^int(ph)):put(h,c,ps)
 vals=list(pool.values());H=np.array([v[0] for v in vals],np.int16);colors=np.array([v[1] for v in vals],np.int16);P=np.array([v[2] for v in vals],np.int32);old=np.array([v[3] for v in vals]);vertices=np.unique(P);where={int(v):i for i,v in enumerate(vertices)};VP=np.array([[where[int(x)] for x in ps] for ps in P]);incident=[[] for _ in vertices]
 for i,(u,v) in enumerate(VP):incident[u].append(i);incident[v].append(i)
 original_parents=set(map(int,OP.ravel()));base=int(z['base']);baseline=201180 if d==26 else 209862;target=(baseline-base)//2+1;assert (baseline-base)%2==0;assert sum(int(v) in original_parents for v in vertices)==target-1;model=cp_model.CpModel();x=[model.new_bool_var(f'edge{i}') for i in range(len(H))];used=[model.new_bool_var(f'removed{i}') for i in range(len(vertices))];double=[model.new_bool_var(f'double{i}') for i in range(len(vertices))]
 for j,es in enumerate(incident):
  deg=sum(x[i] for i in es);model.add(deg==used[j]+double[j]);model.add(double[j]<=used[j]);model.add_hint(used[j],int(vertices[j]) in original_parents);model.add_hint(double[j],int(vertices[j]) in original_parents)
  for i in es:model.add(x[i]<=used[j])
 for i in range(len(H)):model.add_hint(x[i],bool(old[i]))
 HF=H.astype(np.float32);cost=np.abs(T.astype(np.int64)@T.T);edges=0
 for at in range(0,len(H),128):
  D=np.abs(HF[at:at+128]@HF.T);bad=D+cost[colors[at:at+128,None],colors[None,:]]>64;ii,jj=np.where(bad);ii+=at;up=ii<jj;ii=ii[up];jj=jj[up];assert not np.any(old[ii]&old[jj]);edges+=len(ii)
  for i,j in zip(ii,jj):model.add_at_most_one(x[int(i)],x[int(j)])
 model.add(sum(double)>=target);solver=cp_model.CpSolver();solver.parameters.num_search_workers=1;solver.parameters.random_seed=a.seed;solver.parameters.max_time_in_seconds=a.seconds;solver.parameters.symmetry_level=0;solver.parameters.linearization_level=2;solver.parameters.log_search_progress=True;print('MODEL',json.dumps(dict(dimension=d,bodies=len(H),parent_lines=len(vertices),original_bodies=original_count,target_double=target,pair_conflicts=edges,seconds=time.monotonic()-start)),flush=True);st=solver.solve(model);report=dict(dimension=d,status=solver.status_name(st),body_variables=len(H),parent_lines=len(vertices),conflicts=edges,seconds=time.monotonic()-start)
 if st in (cp_model.FEASIBLE,cp_model.OPTIMAL):
  chosen=np.flatnonzero([solver.boolean_value(v) for v in x]);removedlines=vertices[np.flatnonzero([solver.boolean_value(v) for v in used])];lk={v.tobytes():i for i,v in enumerate(L)};anti=np.array([lk[(-v).tobytes()] for v in L]);keep=~blocked.copy();keep[removedlines]=False;keep[anti[removedlines]]=False;LC=np.c_[2*L[keep],np.zeros((int(keep.sum()),d-24),np.int16)];C=np.array([np.r_[s*H[j],t*T[colors[j]]] for j in chosen for s,t in itertools.product([-1,1],repeat=2)],np.int16);RA=np.r_[LC,F,C];RB=np.r_[np.zeros_like(LC),FB,np.zeros_like(C)];assert len(RA)>baseline and len(RA)==base+4*len(chosen)-2*len(removedlines);assert np.all(np.sum(RA.astype(np.int64)**2+3*RB.astype(np.int64)**2,axis=1)==128) and not np.any(np.sum(RA.astype(np.int64)*RB.astype(np.int64),axis=1));src=np.load(f'{Path(__file__).resolve().parent}/data/PackingStar/25D-31D/new_bounds_20261111/D{d}_{baseline}_coordinates.npz')
  if d==26:SA=2*np.r_[src['rational_points'],src['quadratic_A']];SB=2*np.r_[np.zeros_like(src['rational_points']),src['quadratic_B']]
  else:SA=2*src['points'];SB=np.zeros_like(SA)
  sourcekeys={np.r_[v,b].tobytes() for v,b in zip(SA,SB)};keys=[np.r_[v,b].tobytes() for v,b in zip(RA,RB)];assert len(set(keys))==len(RA);added=np.flatnonzero([k not in sourcekeys for k in keys]);newkeys=set(keys);deleted=np.flatnonzero([np.r_[v,b].tobytes() not in newkeys for v,b in zip(SA,SB)]);A64=RA.astype(np.int64);B64=RB.astype(np.int64)
  for j in added:
   fail=above(A64@A64[j]+3*B64@B64[j],A64@B64[j]+B64@A64[j]);fail[j]=False;assert not fail.any(),('changed geometry',j,np.flatnonzero(fail)[:8])
  np.savez_compressed(a.out/'candidate.npz',rational=RA,radical=RB,radicand=3,squared_norm=128,inner_product_limit=64,coordinate_denominator=2);np.savez_compressed(a.out/'exchange.npz',added_A=RA[added],added_B=RB[added],removed_indices=deleted);report.update(count=len(RA),gain=len(RA)-baseline,added=len(added),removed=len(deleted));print('POSITIVE_EXACT_CHANGED_PAIRS',json.dumps(report),flush=True)
 (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('FINISHED',json.dumps(report),flush=True)
if __name__=='__main__':main()
