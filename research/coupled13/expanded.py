"""Expanded coupled search with mixed Q(sqrt3) heads; no antipodal constraints."""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import json,time
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
import experiment as ex
P=Path(__file__).resolve().parent


def pool():
 A=ex.candidates().tolist();B=[[0]*13 for _ in A]
 # A ±3, A ±2, B ±1 at three distinct coordinates: norm16.
 for i in range(13):
  for j in range(13):
   if i==j:continue
   for k in range(13):
    if k in (i,j):continue
    for x,y,z in product((-1,1),repeat=3):
     a=[0]*13;b=[0]*13;a[i]=3*x;a[j]=2*y;b[k]=z;A.append(a);B.append(b)
 # Radical four-sign cap bodies at both heights, norm16.
 for supp in combinations(range(12),4):
  for signs in product((-1,1),repeat=4):
   for h in (-2,2):
    a=[0]*12+[h];b=[0]*13
    for i,s in zip(supp,signs):b[i]=s
    A.append(a);B.append(b)
 # Deterministic sampled integer norm16 shape3,2,1,1,1, not a full enumeration.
 rng=np.random.default_rng(20261104)
 for trial in range(20000):
  ids=rng.choice(13,size=5,replace=False);signs=rng.choice((-1,1),size=5);a=[0]*13
  for i,s,w in zip(ids,signs,(3,2,1,1,1)):a[int(i)]=int(s*w)
  A.append(a);B.append([0]*13)
 return np.array(A,dtype=np.int16),np.array(B,dtype=np.int16)


def filter_pool(C,D,A,B,maxblock=8):
 norms=np.sum(C.astype(np.int64)**2+3*D.astype(np.int64)**2,axis=1);keep=[];blocks=[];hist=Counter()
 for start in range(0,len(C),1024):
  X,Y=C[start:start+1024],D[start:start+1024]
  p=X@A.T+3*(Y@B.T);q=X@B.T+Y@A.T
  dot=p+np.sqrt(3)*q;bad=dot>2*np.sqrt(norms[start:start+len(X),None])+1e-9;counts=bad.sum(axis=1)
  for val,num in zip(*np.unique(counts,return_counts=True)):hist[int(val)]+=int(num)
  for local in np.flatnonzero(counts<=maxblock):
   index=start+int(local);block=np.flatnonzero(bad[local]).tolist();duplicate=False
   for j in block:
    pp,qq=int(p[local,j]),int(q[local,j])
    if ex.engine.sign_sqrt3(pp,qq)>0 and ex.engine.sign_sqrt3(pp*pp+3*qq*qq-int(norms[index])*16,2*pp*qq)==0:duplicate=True;break
   if not duplicate:keep.append(index);blocks.append(block)
 return C[keep],D[keep],norms[keep],blocks,hist


def solve(C,D,norms,blocks,A,B,maxc=3200,seconds=20):
 # Prefer sparse blockers, with a per-blocker cap to avoid drowning mixed
 # candidates in one cheap but unprofitable family. This is a bounded subset.
 order=sorted(range(len(C)),key=lambda i:(len(blocks[i]),tuple(blocks[i]),i))[:maxc]
 C=C[order].astype(np.int64);D=D[order].astype(np.int64);norms=norms[order];blocks=[blocks[i] for i in order];n=len(C)
 p=C@C.T+3*(D@D.T);q=C@D.T+D@C.T;dot=p+np.sqrt(3)*q
 conflict=dot>np.sqrt(norms[:,None]*norms[None,:])/2+1e-9
 rr=[];cc=[];vv=[];ub=[];r=0;old=sorted(set(j for b in blocks for j in b));oldids={j:n+i for i,j in enumerate(old)}
 for i,b in enumerate(blocks):
  for j in b:rr.extend([r,r]);cc.extend([i,oldids[j]]);vv.extend([1,-1]);ub.append(0);r+=1
 for i,j in zip(*np.where(np.triu(conflict,1))):rr.extend([r,r]);cc.extend([int(i),int(j)]);vv.extend([1,1]);ub.append(1);r+=1
 M=coo_matrix((vv,(rr,cc)),shape=(r,n+len(old))).tocsr()
 objective=np.r_[-np.ones(n),np.ones(len(old))]
 result=milp(objective,integrality=np.ones(n+len(old)),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,np.array(ub)),options={'time_limit':seconds,'mip_rel_gap':0})
 chosen=[] if result.x is None else np.flatnonzero(result.x[:n]>.5).tolist();removed=[] if result.x is None else [old[i] for i in np.flatnonzero(result.x[n:]>.5)]
 report={'candidate_subset':n,'mixed_radical_candidates_in_subset':int(np.any(D!=0,axis=1).sum()),'candidate_conflicts':int(np.triu(conflict,1).sum()),'old_points_in_model':len(old),'solver_status':int(result.status),'solver_message':result.message,'gain':len(chosen)-len(removed),'selected':len(chosen),'removed':len(removed),'removed_equator':sum(i<816 for i in removed),'removed_caps_or_poles':sum(i>=816 for i in removed),'finite_numerical_model_gain_bound':-float(result.mip_dual_bound) if getattr(result,'mip_dual_bound',None) is not None else None,'antipodality_constraints':False,'certificate_scope':'numerically constructed bounded conflict model; no rigorous negative bound or exhaustive coverage'}
 if len(chosen)>len(removed):
  retain=[i for i in range(len(A)) if i not in set(removed)];AA=np.concatenate([A[retain],C[chosen]]);BB=np.concatenate([B[retain],D[chosen]])
  fresh=ex.Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':AA.tolist(),'radical_coordinates':BB.tolist()});report['exact_improvement']=ex.engine.verify_artifact(fresh).to_dict();np.savez_compressed(P/'expanded_improvement.npz',A=AA,B=BB)
 return report


def main():
 C,D=pool();print('expanded candidate grammar',len(C),flush=True);source=ex.parse_source();runs=[]
 for left,right in ((0,5),(2,3),(4,5)):
  start=time.time();A,B=ex.seed(left,right,source);X,Y,norms,blocks,hist=filter_pool(C,D,A,B)
  print('seed',left,right,'eligible',len(X),'radical',int(np.any(Y!=0,axis=1).sum()),'hist',sorted(hist.items())[:10],flush=True)
  result=solve(X,Y,norms,blocks,A,B);result.update(left=left,right=right,eligible=len(X),candidate_grammar_size=len(C),blocker_histogram=dict(hist),seconds=time.time()-start)
  runs.append(result);(P/'expanded_results.json').write_text(json.dumps({'outcome':'no SOTA without exact1155+ witness','runs':runs},indent=2));print('result',result,flush=True)
if __name__=='__main__':main()
