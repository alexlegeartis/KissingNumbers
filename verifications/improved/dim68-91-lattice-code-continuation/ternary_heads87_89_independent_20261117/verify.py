from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
import json,hashlib
P=Path(__file__).resolve().parent;c=json.load(open(P.parent/'ternary_qr48_20261117/certificate.json'));G=c['generator_Q48']
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
assert all(sum(a*b for a,b in zip(u,v))%3==0 for u in G for v in G)
results=[]
for n,m,pref,hashval,total in [(87,9,3,'56bc28f5641628611a0246b22f574c70ecd15dbe78c1763eb1273ec9e8c5720c',6973870208),(89,16,4,'7f4b94d67c5896ae9b39c4b6ebc6cfcc5aa12551547956bfa2f8d5019638c1b6',12397783154)]:
 raw=(P/f'input{n}.json').read_bytes();assert hashlib.sha256(raw).hexdigest()==hashval;z=json.loads(raw);rows=z['rows'];assert len(rows)==z['pattern_count']==m
 assert len({tuple(r['pattern'])for r in rows})==m and len({tuple(r['root'])for r in rows})==m
 assert n==84+len(rows[0]['root'])
 for r in rows:
  assert len(r['pattern'])==pref and all(x in [0,1,2]for x in r['pattern'])
  assert len(r['root'])==n-84 and all(x in [-1,0,1]for x in r['root']) and sum(x*x for x in r['root'])==2
 for i,a in enumerate(rows):
  for b in rows[:i+1]:
   d=sum(x!=y for x,y in zip(a['pattern'],b['pattern']));ip=sum(x*y for x,y in zip(a['root'],b['root']))
   assert ip<=2-d
   hmax=3*((42+d)//3)-d
   assert ip>=hmax-43
   for h in range(max(0,15-d),43):
    if (h+d)%3:continue
    inner=42-F(3*h,2)+F(3*ip,2);assert -F(45,2)<=inner<=F(45,2)
 dense=2*m*3**18;assert dense==z['dense_bound']
 pc=Counter(i^j for i in range(n)for j in range(i+1,n));partitions=sum(k*(k-1)//2 for k in pc.values());assert partitions%3==0;supports=partitions//3
 assert supports==sum(a^b^c^d==0 for a,b,c,d in combinations(range(n),4))
 sparse=16*supports+2*n;assert dense+sparse==total
 assert 8*45<=F(45**2,4)
 results.append(dict(n=n,input_sha256=hashval,patterns=m,all_pair_types=m*(m+1)//2,fiber_count=3**18,dense=dense,supports=supports,sparse=sparse,total=total))
(P/'report.json').write_text(json.dumps({'status':'PASS','rows':results},indent=2));print(json.dumps(results,indent=2))
