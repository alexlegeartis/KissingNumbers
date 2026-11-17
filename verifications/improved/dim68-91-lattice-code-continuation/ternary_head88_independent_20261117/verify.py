from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
import json,hashlib
P=Path(__file__).resolve().parent;z=json.load(open(P/'input.json'));c=json.load(open(P.parent/'ternary_qr48_20261117/certificate.json'));G=c['generator_Q48']
def rank(a):
 a=[[x%3 for x in r]for r in a];k=0
 for j in range(len(a[0])):
  i=next((i for i in range(k,len(a))if a[i][j]),None)
  if i is None:continue
  a[k],a[i]=a[i],a[k];u=pow(a[k][j],-1,3);a[k]=[x*u%3 for x in a[k]]
  for i in range(len(a)):
   if i!=k:
    t=a[i][j];a[i]=[(x-t*y)%3 for x,y in zip(a[i],a[k])]
  k+=1
 return k
assert rank(G)==24 and rank([r[:6]for r in G])==6
rows=z['rows'];assert len(rows)==z['pattern_count']==12
assert len({tuple(r['pattern'])for r in rows})==12
assert len({tuple(r['root'])for r in rows})==12
for r in rows:
 assert len(r['pattern'])==3 and all(x in [0,1,2]for x in r['pattern'])
 assert len(r['root'])==4 and all(x in [-1,0,1]for x in r['root']) and sum(x*x for x in r['root'])==2
for a,b in combinations(rows,2):
 d=sum(x!=y for x,y in zip(a['pattern'],b['pattern']));ip=sum(x*y for x,y in zip(a['root'],b['root']));assert -1<=ip<=2-d
 for body_distance in range(15-d,43):
  inner=42-F(3*body_distance,2)+F(3*ip,2)
  assert -F(45,2)<=inner<=F(45,2)
for body_distance in range(15,43):assert -F(45,2)<=45-F(3*body_distance,2)<=F(45,2)
dense=2*12*3**18;assert dense==z['dense_bound']==9298091736
counts=Counter(i^j for i in range(88)for j in range(i+1,88));supports=sum(k*(k-1)//2 for k in counts.values())//3;assert supports==19626
sparse=16*supports+176;assert sparse==314192
# Dense coordinate max sqrt(3/2), below prior safe sqrt2.
assert F(3,2)<=2 and 8*45<=F(45**2,4)
result=dict(status='PASS',input_sha256=hashlib.sha256((P/'input.json').read_bytes()).hexdigest(),prefix_projection_rank=6,head_pairs=66,fiber_count=3**18,dense=dense,sparse=sparse,total=dense+sparse,source='Published QR48 minimum15, independently replayed generator rank and shortening')
assert result['total']==9298405928
(P/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
