"""Exact binary-torsion and low-residue-coset search for unresolved Z4 modules."""
from pathlib import Path
from itertools import combinations
import json,time,hashlib
import numpy as np
OUT=Path(__file__).resolve().parent

def wt(w):return sum(min(x,4-x) for x in w)
def basis(rows):
 pivots={}
 for i,row in enumerate(rows):
  label=1<<i
  while row:
   p=row.bit_length()-1
   if p in pivots:row^=pivots[p][0];label^=pivots[p][1]
   else:pivots[p]=(row,label);break
 return pivots

def reduce(row,pivots):
 label=0
 for p in sorted(pivots,reverse=True):
  if row>>p&1:row^=pivots[p][0];label^=pivots[p][1]
 return row,label

def refine(row):
 q=row['q'];n=row['n'];k=row['rank'];g=row['generator_polynomial'];target=row['target_distance']
 G=np.zeros((k,n),dtype=np.int64)
 for i in range(k):G[i,i:i+len(g)]=g
 binary=[sum(int(x%2)<<j for j,x in enumerate(v)) for v in G]
 words=[0]
 for b in binary:words += [v^b for v in words]
 weights=np.array([w.bit_count() for w in words]);weights[0]=n+1
 best_label=int(weights.argmin());minwt=int(weights[best_label])
 result=dict(n=n,q=q,rank=k,factor_subset=row['factor_subset'],binary_words_checked=len(words),binary_min_weight=minwt)
 def finish(coeff,method,details):
  W=(np.array(coeff,dtype=np.int64)@G%4).tolist();weight=wt(W)
  bad={}
  for s in range(4):
   a=(s*sum(W))%4
   if weight+wt([a])<target:bad[str(s)]=dict(extension_value=a,extended_weight=weight+wt([a]))
  return dict(**result,method=method,message=coeff,word=W,base_weight=weight,target_distance=target,
    excludes_every_extra_coordinate_value=weight+2<target,scalar_witnesses=bad,details=details)
 torsion=finish([2*(best_label>>i&1) for i in range(k)],'complete binary torsion enumeration',{})
 if torsion['excludes_every_extra_coordinate_value']:return torsion
 # Search residual classes in order of weight, at most64 classes. The even-coordinate
 # syndrome is exact; choose an error whose addition places the target inside the
 # projected binary residue code. Its preimage gives the 2*b correction to a.
 count=0
 for label in np.argsort(weights)[:8192 if n==31 else 64]:
  label=int(label);a=[label>>i&1 for i in range(k)];v=np.array(a)@G%4
  outside=[j for j in range(n) if not(v[j]&1)]
  projected=[sum(int(G[i,j]%2)<<t for t,j in enumerate(outside))for i in range(k)]
  pivots=basis(projected)
  baseline=sum(int(v[j]//2)<<t for t,j in enumerate(outside))
  maxerror=min(4,(target-1-int(weights[label]))//2)
  for ewt in range(maxerror+1):
   for support in combinations(range(len(outside)),ewt):
    error=sum(1<<t for t in support)
    residual,solution=reduce(baseline^error,pivots);count+=1
    if residual:continue
    coeff=[(a[i]+2*(solution>>i&1))%4 for i in range(k)]
    found=finish(coeff,'binary projected syndrome exact decoding',dict(residue_label=label,residue_weight=int(weights[label]),outside_error_weight=ewt,syndromes_tested=count))
    if found['excludes_every_extra_coordinate_value'] or len(found['scalar_witnesses'])==4:return found
    if found['base_weight']<torsion['base_weight']:torsion=found
 return {**torsion,'search_scope':'8192 smallest residue labels for length31, 64 otherwise; outside errors of weight<=4; no completeness claim beyond binary enumeration'}

def main():
 start=time.monotonic();data=json.loads((OUT/'composite_screen.json').read_text())
 rows=[r for r in data['examined']if r['status']=='unresolved_after_bounded_message_scan'and r['q']==4]
 results=[]
 for row in rows:
  r=refine(row);results.append(r)
  print(r['n'],r['rank'],r['factor_subset'],'base',r['base_weight'],'all extra',r['excludes_every_extra_coordinate_value'],'bad scalars',list(r['scalar_witnesses']),flush=True)
 receipt=dict(results=results,elapsed_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_sha256=hashlib.sha256((OUT/'composite_screen.json').read_bytes()).hexdigest())
 (OUT/'refine_residues.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('elapsed',receipt['elapsed_seconds'])
if __name__=='__main__':main()
