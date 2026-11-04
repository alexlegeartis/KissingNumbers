"""Targeted exact compatibility search on bounded joint removal sets."""
from itertools import combinations
from pathlib import Path
from functools import lru_cache
import json,time,random
import numpy as np
import expanded as xp
P=Path(__file__).resolve().parent


def find_trade(C,D,norms,blocks,seconds=20):
 start=time.monotonic();lookup={}
 for i,b in enumerate(blocks):lookup.setdefault(tuple(b),[]).append(i)
 unique=list(lookup);sets=list(unique)
 rng=random.Random(20261104)
 for trial in range(12000):
  a,b=rng.sample(unique,2);joint=tuple(sorted(set(a)|set(b)))
  if len(joint)<=12:sets.append(joint)
 sets=list(dict.fromkeys(sets));checked=0;eligible_count=0;max_unlocked=0;nodes=0
 @lru_cache(maxsize=300000)
 def conflict(i,j):
  a=int(C[i]@C[j]+3*(D[i]@D[j]));b=int(C[i]@D[j]+D[i]@C[j])
  return not xp.ex.engine.compatible_dot(a,b,int(norms[i])*int(norms[j]))
 for removed in sets:
  if time.monotonic()-start>seconds:break
  k=len(removed);eligible=[]
  for count in range(k+1):
   for subset in combinations(removed,count):eligible.extend(lookup.get(subset,()))
  checked+=1;max_unlocked=max(max_unlocked,len(eligible))
  if len(eligible)<=k:continue
  eligible_count+=1
  # Find k+1 compatible insertions. Prune only with an exact remaining count.
  target=k+1
  def dfs(chosen,remaining):
   nonlocal nodes
   nodes+=1
   if len(chosen)>=target:return chosen
   if len(chosen)+len(remaining)<target or nodes>1500000:return None
   while remaining:
    v=remaining.pop()
    nexts=[u for u in remaining if not conflict(v,u)]
    answer=dfs(chosen+[v],nexts)
    if answer is not None:return answer
    if len(chosen)+len(remaining)<target:break
   return None
  answer=dfs([],eligible[:])
  if answer is not None:return {'removed':list(removed),'selected':answer,'gain':len(answer)-k,'tested_removal_sets':checked,'max_unlocked':max_unlocked,'search_nodes':nodes,'seconds':time.monotonic()-start},answer,list(removed)
  if nodes>1500000:break
 return {'gain':0,'tested_removal_sets':checked,'total_proposed_removal_sets':len(sets),'nontrivial_sets':eligible_count,'max_unlocked':max_unlocked,'search_nodes':nodes,'seconds':time.monotonic()-start,'scope':'single candidate blocker sets and bounded sampled pair unions, max12 old removals; exact compatibility but no general exclusion'},None,None


def main():
 C,D=xp.pool();source=xp.ex.parse_source();runs=[]
 for left,right in ((0,5),(2,3),(4,5)):
  A,B=xp.ex.seed(left,right,source);X,Y,norms,blocks,hist=xp.filter_pool(C,D,A,B,maxblock=8)
  result,chosen,removed=find_trade(X,Y,norms,blocks,seconds=22)
  result.update(left=left,right=right,eligible=len(X),candidate_pool=len(C),antipodality_constraints=False)
  if chosen is not None:
   keep=[i for i in range(len(A)) if i not in set(removed)];AA=np.concatenate([A[keep],X[chosen]]);BB=np.concatenate([B[keep],Y[chosen]])
   fresh=xp.ex.Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':AA.tolist(),'radical_coordinates':BB.tolist()});result['exact_improvement_receipt']=xp.ex.engine.verify_artifact(fresh).to_dict();np.savez_compressed(P/'small_trade_improvement.npz',A=AA,B=BB)
  runs.append(result);(P/'small_trade_results.json').write_text(json.dumps({'runs':runs},indent=2));print(result,flush=True)
if __name__=='__main__':main()
