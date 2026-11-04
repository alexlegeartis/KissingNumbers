"""Coupled short-dual point search against all three retained 40D shells.

All existing points may be deleted; shared deletions are charged once. A
blocker cutoff bounds the finite neighbourhood, never the kissing number.
"""
from pathlib import Path
from math import isqrt
import argparse,hashlib,json,sys,time
import numpy as np
from scipy.sparse import csr_matrix,coo_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.optimize import milp,Bounds,LinearConstraint
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--blockers',type=int,default=64)
 ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seconds',type=float,default=120)
 args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 P=ROOT/'verifications/improved/dim40-affine-shells/source'
 R=ROOT/'verifications/improved/dim41-46-p48-heads/partial43'
 Q=np.array(json.loads((P/'gram.json').read_text()),np.int64)
 data=np.load(R/'body_complement.npz');X=np.load(R/'tiny_bodies.npz')['H'].astype(np.int64)
 if not np.array_equal(Q,data['q']):raise ValueError('Shared exact Gram')
 qi3=np.rint(3*np.linalg.inv(Q)).astype(np.int64)
 if not np.array_equal(Q@qi3,3*np.eye(40,dtype=np.int64)):raise ValueError('Exact inverse')
 Y=(X@data['b'].T)@qi3
 if np.any(Y%288):raise ValueError('Integral short-dual coefficients')
 Y//=288
 if len(Y)!=30240 or not np.all(np.einsum('ij,ij->i',Y@Q,Y)==30):raise ValueError('Complete tiny source norms')
 z=np.load(P/'auxiliary.npz');A=np.load(P/'core.npz')['coefficients']
 B=z['norm36'];C=z['norm30'];known={row.astype(np.int64).tobytes() for row in C}
 Y=np.array([row for row in Y if row.tobytes() not in known],np.int64)
 if len(Y)!=29106:raise ValueError('Baseline tiny membership')
 initial_count=len(Y);L=Y@Q;active=np.arange(len(Y));blockers=[[] for _ in Y]
 old_count=len(A)+len(B)+len(C);checks=pruned=0;certificates=[]
 rng=np.random.default_rng(20261104)
 for name,old,norm,offset in [('medium',B,36,len(A)),('small',C,30,len(A)+len(B)),('core',A,54,0)]:
  bound=int(abs(L).max(initial=0))*int(np.abs(old.astype(np.int64)).sum(axis=1).max(initial=0))
  if bound>=2**24:raise ValueError(('Exact binary32 bound',bound))
  threshold=isqrt(30*norm)//2;order=rng.permutation(len(old));scanned=0
  for at in range(0,len(old),256):
   if not len(active):break
   ids=order[at:at+256];dots=L[active].astype(np.float32)@old[ids].astype(np.float32).T
   bad=dots>threshold;checks+=len(active)*len(ids);scanned+=len(ids)
   counts=np.array([len(blockers[int(i)]) for i in active])+bad.sum(axis=1)
   good=counts<=args.blockers
   for row,i in enumerate(active):
    if good[row]:blockers[int(i)].extend((offset+ids[np.flatnonzero(bad[row])]).tolist())
    else:blockers[int(i)]=[]
   pruned+=int((~good).sum());active=active[good]
  certificates.append({'layer':name,'norm':norm,'threshold':threshold,'source_points_scanned':scanned,
    'survivors':len(active),'arithmetic_absolute_bound':bound})
  print(json.dumps(certificates[-1]),flush=True)
 selected=[blockers[int(i)] for i in active];ptr=np.r_[0,np.cumsum(list(map(len,selected)))]
 edges=np.array([j for row in selected for j in row],np.int64)
 graph=csr_matrix((np.ones(len(edges),np.int8),edges,ptr),shape=(len(active),old_count))
 matching=maximum_bipartite_matching(graph,perm_type='column')
 if len(set(map(int,matching[matching>=0])))!=int((matching>=0).sum()):raise ValueError('Matching injection')
 for i,j in enumerate(matching):
  if j>=0 and int(j) not in selected[i]:raise ValueError('Exact blocker membership')
 report={'dimension':40,'baseline':old_count,'generated_short_dual':30240,'new_candidates':initial_count,
  'blocker_cutoff':args.blockers,'restricted_candidates':len(active),'pruned':pruned,
  'layers':certificates,'dot_products':checks,'matching_size':int((matching>=0).sum()),
  'matching_deficiency':int((matching<0).sum()),'seconds':time.monotonic()-start,'new_sota':False,
  'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
   [P/'gram.json',P/'auxiliary.npz',P/'core.npz',R/'tiny_bodies.npz',R/'body_complement.npz']},
  'scope':'All baseline layers mutable. Complete exact blockers for admitted candidates. A saturating injection excludes gain only from this finite cutoff pool; no global bound.'}
 np.savez_compressed(args.out/'pool.npz',candidates=Y[active],candidate_source_indices=active,
  blocker_indptr=ptr,blocker_indices=edges,matching=matching,gram=Q)
 if np.any(matching<0):
  # The matching deficit ignores new/new conflicts. Solve the coupled model;
  # a deficit alone is never reported as a geometric improvement.
  used=sorted(set(map(int,edges)));variables={j:len(active)+k for k,j in enumerate(used)}
  rr=[];cc=[];vv=[];upper=[]
  for i,row in enumerate(selected):
   for j in row:
    t=len(upper);rr.extend((t,t));cc.extend((i,variables[j]));vv.extend((1,-1));upper.append(0)
  Z=Y[active];G=(Z@Q)@Z.T
  for i,j in zip(*np.where(np.triu(G>15,1))):
   t=len(upper);rr.extend((t,t));cc.extend((int(i),int(j)));vv.extend((1,1));upper.append(1)
  mat=coo_matrix((vv,(rr,cc)),shape=(len(upper),len(active)+len(used))).tocsc()
  result=milp(np.r_[-np.ones(len(active)),np.ones(len(used))],integrality=np.ones(len(active)+len(used)),
   bounds=Bounds(0,1),constraints=LinearConstraint(mat,-np.inf,upper),options={'time_limit':args.seconds,'mip_rel_gap':0})
  additions=[] if result.x is None else np.flatnonzero(result.x[:len(active)]>.5).tolist()
  deletions=[] if result.x is None else [used[i] for i in np.flatnonzero(result.x[len(active):]>.5)]
  report['optimization']={'status':int(result.status),'gain':len(additions)-len(deletions),'added':additions,'removed':deletions,
   'scope':'Numerical proposal only; independent exact final geometry required'}
  np.savez_compressed(args.out/'proposal.npz',added=Z[additions],removed_indices=np.array(deletions,np.int64),gram=Q)
 report['seconds']=time.monotonic()-start
 (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
