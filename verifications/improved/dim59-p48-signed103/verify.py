"""Exact integer verification of the 59D six-to-eight pure-tail exchange."""
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import json,sys
assert __debug__
P=Path(__file__).resolve().parent
base=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else P.parent/'dim59-p48-signed101'
def rd(p,sha=None):
 b=p.read_bytes()
 if sha:assert sha256(b).hexdigest()==sha
 return json.loads(b)
b=rd(base/'claims.json')
old=rd(base/'signed_winner.json','97d952397d5eebf8ffb53a35344131d6be4bbdf30cd54fbf4eda76b34c8baad8')
T=rd(base/'tet.json','cc4151354d9b6eaf6551275a2d9424eb5b8eba11e0a4c3bf5db0bec31e9aa053')['59']['supports']
Q=rd(base/'weighted.json','68f8687a602a51559bdbfecac75cacf90097a47e49b70bfdf711091c22527cf7')
d=rd(P/'witness.json'); c=rd(P/'claims.json')
def dot(x,y):return sum(a*b for a,b in zip(x,y))
D=set(d['removed']);A=old['poles'];N=d['added'];w=[50]*8+[49]*3
assert len(A)==101 and not old['removed_weighted']
assert len(D)==len(d['removed'])==6 and all(type(i) is int and 0<=i<len(A)for i in D)
assert len(N)==8 and all(len(r)==11 and all(type(i) is int for i in r)and any(r)for r in N)
R=[p for i,p in enumerate(A)if i not in D]
for r in N:
 n=dot(r,r)
 assert 3*dot(list(map(abs,r)),w)**2<=27203*n
 assert all(5*sum(abs(r[i])for i in t)**2<=6*n for t in T)
 for i,*z in Q['weighted_tails']:
  v=dot(r,z)
  assert v<=0 or 4*v*v<=Q['scale']**2*n
for a,b1 in list(combinations(N,2))+[(a,b1)for a in N for b1 in R]:
 z=dot(a,b1)
 assert z<=0 or 4*z*z<=dot(a,a)*dot(b1,b1)
previous=b['lower_bounds']['59']
assert c['previous']['59']==previous==56969164
assert c['lower_bounds']['59']==previous+len(N)-len(D)==56969166
print('PASS 59D: 56,969,166; all inherited mixed stock retained')
