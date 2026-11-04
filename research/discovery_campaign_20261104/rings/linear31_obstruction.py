"""Search for finite inconsistent necessary conditions for arbitrary linear extensions."""
from pathlib import Path
from itertools import combinations
import json,time,hashlib
import numpy as np
from refine_residues import basis, reduce, wt
OUT=Path(__file__).resolve().parent

def inconsistent(constraints,k):
 M=[list(c['message'])+[2] for c in constraints]
 at=0
 for j in range(k):
  pivot=next((i for i in range(at,len(M))if M[i][j]&1),None)
  if pivot is None:continue
  M[at],M[pivot]=M[pivot],M[at]
  unit=M[at][j]
  M[at]=[x*unit%4 for x in M[at]]
  for i in range(len(M)):
   if i!=at:
    a=M[i][j];M[i]=[(x-a*y)%4 for x,y in zip(M[i],M[at])]
  at+=1
 # All remaining coefficients are even. Odd RHS is impossible; otherwise
 # reduce residual equations modulo2 after dividing by2.
 B=[]
 for row in M[at:]:
  if row[-1]&1:return True
  B.append([x//2 for x in row])
 at=0
 for j in range(k):
  pivot=next((i for i in range(at,len(B))if B[i][j]),None)
  if pivot is None:continue
  B[at],B[pivot]=B[pivot],B[at]
  for i in range(len(B)):
   if i!=at and B[i][j]:B[i]=[x^y for x,y in zip(B[i],B[at])]
  at+=1
 return any(not any(r[:-1])and r[-1]for r in B)

def main():
 start=time.monotonic();data=json.loads((OUT/'composite_screen.json').read_text());results=[]
 for row in data['examined']:
  if (row['n'],row['q'],row['rank'])!=(31,4,15)or row['status']!='unresolved_after_bounded_message_scan':continue
  k=15;n=31;g=row['generator_polynomial'];G=np.zeros((k,n),dtype=np.int64)
  for i in range(k):G[i,i:i+len(g)]=g
  binary=[sum(int(x%2)<<j for j,x in enumerate(v))for v in G];words=[0]
  for b in binary:words += [v^b for v in words]
  order=sorted(range(1,len(words)),key=lambda l:words[l].bit_count())
  constraints=[];found=False;tested=0
  for label in order[:8192]:
   rw=words[label].bit_count()
   if rw>14:break
   a=[label>>i&1 for i in range(k)];v=np.array(a)@G%4
   outside=[j for j in range(n)if not(v[j]&1)]
   pivots=basis([sum(int(G[i,j]%2)<<t for t,j in enumerate(outside))for i in range(k)])
   baseline=sum(int(v[j]//2)<<t for t,j in enumerate(outside))
   decoded=False
   for ewt in range((14-rw)//2+1):
    for support in combinations(range(len(outside)),ewt):
     tested+=1;error=sum(1<<t for t in support);residual,solution=reduce(baseline^error,pivots)
     if residual:continue
     c=[(a[i]+2*(solution>>i&1))%4 for i in range(k)];w=(np.array(c)@G%4).tolist()
     if wt(w)!=14:raise ValueError(('unexpected weight',wt(w)))
     constraints.append(dict(message=c,word=w,base_weight=14,required_extension_value=2))
     decoded=True;break
    if decoded:break
   if constraints and inconsistent(constraints,k):found=True;break
   if len(constraints)>=256:break
  multipliers=None
  if found:
   labels=np.arange(4**len(constraints),dtype=np.int64)
   coefficients=np.column_stack([(labels//4**i)%4 for i in range(len(constraints))])
   matrix=np.array([c['message']for c in constraints],dtype=np.int64)
   valid=np.flatnonzero(np.all(coefficients@matrix%4==0,axis=1)&(2*coefficients.sum(axis=1)%4==2))
   if len(valid):multipliers=coefficients[valid[0]].tolist()
  result=dict(contradiction_multipliers=multipliers,n=n,q=4,rank=k,factor_subset=row['factor_subset'],generator_polynomial=g,target_distance=16,
    constraints=constraints,inconsistent=found,residue_labels_screened=order.index(label)+1,syndromes_tested=tested)
  results.append(result);print('subset',row['factor_subset'],'inconsistent',found,'constraints',len(constraints),'syndromes',tested,flush=True)
 receipt=dict(results=results,elapsed_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (OUT/'linear31_obstruction.json').write_text(json.dumps(receipt,indent=2)+'\n');print('elapsed',receipt['elapsed_seconds'])
if __name__=='__main__':main()
