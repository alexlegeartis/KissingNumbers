"""New fixed-norm16 unequal-height alphabets; full blocker distributions first."""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import json,time
import numpy as np
import expanded as xp
P=Path(__file__).resolve().parent


def families():
 # h=+/-1, body norm15: one +/-2, eleven +/-1.
 A=[];B=[]
 for heavy in range(12):
  for signs in product((-1,1),repeat=12):
   for h in (-1,1):
    row=list(signs)+[h];row[heavy]*=2;A.append(row)
 yield 'h1_full_2_ones11',np.array(A,dtype=np.int16),None
 # h=+/-1, body norm15:3,2,1,1.
 A=[]
 for i in range(12):
  for j in range(12):
   if i==j:continue
   for pair in combinations([k for k in range(12) if k not in (i,j)],2):
    for signs in product((-1,1),repeat=4):
     for h in (-1,1):
      row=[0]*12+[h]
      for k,s,w in zip((i,j)+pair,signs,(3,2,1,1)):row[k]=s*w
      A.append(row)
 yield 'h1_sparse_3_2_1_1',np.array(A,dtype=np.int16),None
 # h=+/-sqrt3, integer body norm13:3,2.
 A=[];B=[]
 for i in range(12):
  for j in range(12):
   if i==j:continue
   for x,y,h in product((-1,1),repeat=3):
    row=[0]*13;rad=[0]*13;row[i]=3*x;row[j]=2*y;rad[-1]=h;A.append(row);B.append(rad)
 yield 'h_sqrt3_sparse_3_2',np.array(A,dtype=np.int16),np.array(B,dtype=np.int16)
 # h=+/-3, body norm7:2,1,1,1.
 A=[]
 for i in range(12):
  for supp in combinations([j for j in range(12) if j!=i],3):
   for signs in product((-1,1),repeat=4):
    for h in (-3,3):
     row=[0]*12+[h]
     for k,s,w in zip((i,)+supp,signs,(2,1,1,1)):row[k]=s*w
     A.append(row)
 yield 'h3_sparse_2_1_1_1',np.array(A,dtype=np.int16),None
 # h=+/-3, body norm7: seven +/-1.
 A=[]
 for supp in combinations(range(12),7):
  for signs in product((-1,1),repeat=7):
   for h in (-3,3):
    row=[0]*12+[h]
    for k,s in zip(supp,signs):row[k]=s
    A.append(row)
 yield 'h3_flat_ones7',np.array(A,dtype=np.int16),None


def profile_and_select(C,D,A,B,limit=300,rng=None):
 old=A+np.sqrt(3)*B;counts=np.empty(len(C),dtype=np.int16)
 for start in range(0,len(C),1024):
  x=C[start:start+1024].astype(float)
  if D is not None:x+=np.sqrt(3)*D[start:start+1024]
  counts[start:start+len(x)]=(x@old.T>8+1e-9).sum(axis=1)
 hist={int(k):int(v) for k,v in zip(*np.unique(counts,return_counts=True))}
 # Explicitly sample the cheapest300 from every family, including families
 # whose minimum blocker cost exceeds8 or20. No <=8 admission shortcut.
 shuffled=rng.permutation(len(C));order=shuffled[np.argsort(counts[shuffled],kind='stable')[:limit]]
 X=C[order];Y=np.zeros_like(X) if D is None else D[order]
 bad=(X+np.sqrt(3)*Y)@old.T>8+1e-9;blocks=[np.flatnonzero(row).tolist() for row in bad]
 return X,Y,blocks,{'total':len(C),'minimum_blockers':int(counts.min()),'maximum_blockers':int(counts.max()),'blocker_histogram':hist,'selected':len(order),'selected_blocker_min':min(map(len,blocks)),'selected_blocker_max':max(map(len,blocks)),'all_points_fixed_squared_norm':16,'admission_rule':'300 cheapest with deterministic shuffled tie breaks; no blocker cutoff'}


def main():
 source=xp.ex.parse_source();runs=[]
 for left,right in ((0,5),(2,3),(4,5)):
  start=time.time();A,B=xp.ex.seed(left,right,source);Xs=[];Ys=[];blocks=[];profiles={};rng=np.random.default_rng(20261104)
  for name,C,D in families():
   X,Y,b,profile=profile_and_select(C,D,A,B,rng=rng);Xs.append(X);Ys.append(Y);blocks.extend(b);profiles[name]=profile
   print(left,right,name,profile['total'],profile['minimum_blockers'],profile['selected_blocker_max'],flush=True)
  X=np.concatenate(Xs);Y=np.concatenate(Ys);norms=np.full(len(X),16,dtype=np.int64)
  result=xp.solve(X,Y,norms,blocks,A,B,maxc=len(X),seconds=20)
  result.update(left=left,right=right,families=profiles,full_alphabet_points=sum(p['total'] for p in profiles.values()),seconds=time.time()-start)
  runs.append(result);(P/'unequal_height_results.json').write_text(json.dumps({'outcome':'no SOTA without exact1155+ witness','runs':runs},indent=2));print('unequal result',left,right,result['gain'],result['solver_status'],result['seconds'],flush=True)
if __name__=='__main__':main()
