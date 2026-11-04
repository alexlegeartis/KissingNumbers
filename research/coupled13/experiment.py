"""Primary-formula free-zone seeds and joint point-removal/replacement search."""
from pathlib import Path
from itertools import combinations, product
from collections import Counter
import sys,re,json,time,hashlib
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from discovery.core import Artifact,RunContext
from discovery.engines import signed_layers as engine
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix
P=Path(__file__).resolve().parent


def parse_source():
 A=[];B=[]
 for line in (P/'ze99.md').read_text().splitlines():
  tok=line.strip(' `').split(',')
  if len(tok)!=13:continue
  a=[];b=[]
  for t in tok:
   t=t.strip()
   if re.fullmatch(r'-?\d+',t):a.append(int(t));b.append(0)
   elif re.fullmatch(r'-?\d+\*sqrt\(3\)',t):a.append(0);b.append(int(t.split('*')[0]))
   else:break
  if len(a)==13:A.append(a);B.append(b)
 if len(A)!=1154:raise ValueError(('source count',len(A)))
 return np.array(A,dtype=np.int64),np.array(B,dtype=np.int64)


def zone(variant):
 result=[]
 for support in ((0,1,2,3),(0,1,4,5),(2,3,4,5)):
  for signs in product((-1,1),repeat=4):
   row=[0]*6
   for i,s in zip(support,signs):row[i]=s
   if support[0]==0:
    pair=tuple(row[:2]);swap=None
    if variant==1 and pair==(1,1):swap=(2,5)
    if variant==2 and row[0]==1:swap=(2,5)
    if variant==3:
     if pair==(1,-1):swap=(2,5)
     elif pair==(1,1):swap=(2,4)
    if variant==4 and row[0]!=row[1]:swap=(2,5)
    if variant==5:
     if row[0]!=row[1]:swap=(2,5)
     elif pair==(1,1):swap=(2,4)
    if swap:row[swap[0]],row[swap[1]]=row[swap[1]],row[swap[0]]
   result.append(row)
 V=np.array(result,dtype=np.int64);G=V@V.T;np.fill_diagonal(G,-4)
 if len(set(map(tuple,V)))!=48 or G.max()>2:raise ValueError('zone formula')
 return V


def equator(left,right):
 rows=[]
 for color in range(5):
  edges=[(5,color),((color+1)%5,(color-1)%5),((color+2)%5,(color-2)%5)]
  for edgea in edges:
   for edgeb in edges:
    support=list(edgea)+[i+6 for i in edgeb]
    for signs in product((-1,1),repeat=4):
     row=[0]*12
     for i,s in zip(support,signs):row[i]=s
     rows.append(row)
 rows.extend(np.pad(zone(left),((0,0),(0,6))).tolist())
 rows.extend(np.pad(zone(right),((0,0),(6,0))).tolist())
 E=np.array(rows,dtype=np.int64);G=E@E.T;np.fill_diagonal(G,-4)
 if len(E)!=816 or len(set(map(tuple,E)))!=816 or G.max()>2:raise ValueError('equator formula')
 return E


def seed(left,right,source):
 A0,B0=source;keep=A0[:,12]!=0
 A=np.concatenate([np.c_[2*equator(left,right),np.zeros(816,dtype=int)],A0[keep]])
 B=np.concatenate([np.zeros((816,13),dtype=int),B0[keep]])
 return A,B


def candidates():
 rows=[]
 # Include all signed supports up to five across every coordinate, hence both
 # height-zero and nonzero directions, independently for both orientations.
 for k in range(1,6):
  for supp in combinations(range(13),k):
   for signs in product((-1,1),repeat=k):
    row=[0]*13
    for i,s in zip(supp,signs):row[i]=s
    rows.append(row)
 # Body supports with height +/-2 couple cap and equator blocker sets.
 for k in range(1,6):
  for supp in combinations(range(12),k):
   for signs in product((-1,1),repeat=k):
    for h in (-2,2):
     row=[0]*12+[h]
     for i,s in zip(supp,signs):row[i]=s
     rows.append(row)
 return np.array(rows,dtype=np.int16)


def filter_pool(C,A,B,maxblock=5):
 norms=np.sum(C*C,axis=1).astype(np.int64);inds=[];blocks=[];hist=Counter()
 for start in range(0,len(C),1024):
  X=C[start:start+1024]
  G=X@A.T+np.sqrt(3)*(X@B.T)
  bad=G>2*np.sqrt(norms[start:start+len(X),None])+1e-9
  counts=bad.sum(axis=1)
  for value,num in zip(*np.unique(counts,return_counts=True)):hist[int(value)]+=int(num)
  for local in np.flatnonzero(counts<=maxblock):
   index=start+int(local);block=np.flatnonzero(bad[local]).tolist()
   # All candidates must be new normalized points. Equality has cosine1.
   duplicate=False
   for j in block:
    p=int(C[index]@A[j]);q=int(C[index]@B[j])
    if engine.sign_sqrt3(p,q)>0 and engine.sign_sqrt3(p*p+3*q*q-int(norms[index])*16,2*p*q)==0:duplicate=True;break
   if not duplicate:inds.append(index);blocks.append(block)
 return C[inds],norms[inds],blocks,hist


def solve(C,norms,blocks,A,B,seconds=10,maxc=2000):
 if len(C)>maxc:
  order=sorted(range(len(C)),key=lambda i:(len(blocks[i]),tuple(blocks[i]),i))[:maxc]
  C=C[order];norms=norms[order];blocks=[blocks[i] for i in order]
 # Candidate mutual compatibility is rational: positive dot needs4dot²<=normproduct.
 G=C.astype(np.int64)@C.astype(np.int64).T
 conflict=(G>0)&(4*G*G>norms[:,None]*norms[None,:])
 rr=[];cc=[];vv=[];ub=[];r=0;n=len(C);old=sorted(set(j for b in blocks for j in b));oldids={j:n+i for i,j in enumerate(old)}
 for i,b in enumerate(blocks):
  for j in b:rr.extend([r,r]);cc.extend([i,oldids[j]]);vv.extend([1,-1]);ub.append(0);r+=1
 for i,j in zip(*np.where(np.triu(conflict,1))):rr.extend([r,r]);cc.extend([int(i),int(j)]);vv.extend([1,1]);ub.append(1);r+=1
 if n==0:return {'candidates':0,'gain':0},None
 objective=np.r_[-np.ones(n),np.ones(len(old))]
 M=coo_matrix((vv,(rr,cc)),shape=(r,n+len(old))).tocsr()
 result=milp(objective,integrality=np.ones(n+len(old)),bounds=Bounds(0,1),constraints=LinearConstraint(M,-np.inf,np.array(ub)),options={'time_limit':seconds,'mip_rel_gap':0.0})
 if result.x is None:chosen=[];removed=[]
 else:
  chosen=np.flatnonzero(result.x[:n]>.5).tolist();removed=[old[i] for i in np.flatnonzero(result.x[n:]>.5)]
 gain=len(chosen)-len(removed)
 report={'candidates':n,'candidate_conflicts':int(np.triu(conflict,1).sum()),'old_points_in_model':len(old),'solver_status':int(result.status),'solver_message':result.message,'gain':gain,'selected':len(chosen),'removed':len(removed),'removed_equator':sum(i<816 for i in removed),'removed_caps_or_poles':sum(i>=816 for i in removed),'finite_model_gain_upper_bound':-float(result.mip_dual_bound) if getattr(result,'mip_dual_bound',None) is not None else None,'antipodality_constraints':False}
 if gain>0:
  keep=[i for i in range(len(A)) if i not in set(removed)];AA=np.concatenate([A[keep],C[chosen]]);BB=np.concatenate([B[keep],np.zeros((len(chosen),13),dtype=int)])
  return report,(AA,BB)
 return report,None


def main():
 source=parse_source();np.savez_compressed(P/'source1154.npz',A=source[0],B=source[1])
 sourceartifact=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':source[0].tolist(),'radical_coordinates':source[1].tolist()})
 receipt=engine.verify_artifact(sourceartifact).to_dict();(P/'source_exact.json').write_text(json.dumps(receipt,indent=2))
 C=candidates();print('candidate grammar',len(C),flush=True);runs=[]
 pairs=[(0,2),(0,3),(0,4),(0,5),(2,2),(2,3),(3,4),(4,5)]
 for left,right in pairs:
  start=time.time();A,B=seed(left,right,source)
  artifact=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':A.tolist(),'radical_coordinates':B.tolist()})
  report=engine.verify_artifact(artifact).to_dict()
  np.savez_compressed(P/f'seed_{left}_{right}.npz',A=A,B=B)
  X,norms,blocks,hist=filter_pool(C,A,B,maxblock=5)
  print('seed',left,right,'filtered',len(X),'hist',sorted(hist.items())[:8],flush=True)
  result,witness=solve(X,norms,blocks,A,B,seconds=12,maxc=2400)
  result.update(left=left,right=right,seed_cardinality=1154,seed_exact_receipt=report,candidate_grammar_size=len(C),candidate_blocker_histogram=dict(hist),eligible_candidates=len(X),seconds=time.time()-start)
  if witness is not None:
   aa,bb=witness;fresh=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':aa.tolist(),'radical_coordinates':bb.tolist()});result['improvement_exact_receipt']=engine.verify_artifact(fresh).to_dict();np.savez_compressed(P/f'improvement_{left}_{right}.npz',A=aa,B=bb)
  runs.append(result);(P/'results.json').write_text(json.dumps({'outcome':'new SOTA only if an exact1155+ receipt exists','runs':runs},indent=2));print('result',result['gain'],result['solver_status'],result['seconds'],flush=True)
if __name__=='__main__':main()
