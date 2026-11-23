"""Exact replay of orthogonal A1 and A2 sidecars in 81D and 82D."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations,product
import json,runpy,sys
import numpy as np
assert __debug__
here=Path(__file__).resolve().parent
up=Path(sys.argv[1]).resolve()if len(sys.argv)>1 else here.parent/'dim80-gamma72-partial-fibres'
base=up.parent/'dim77-79-gamma72-seven-charts'
pins={
  "verify.py": "c90daad443667a237429129bbf1a7e4cec67f7a224498409a1c47066240e5e26",
  "witness.json": "76271222491054d5596ed16c245371a7ac72053bb9370513ebe2d7c0b982b3d1"
}
assert all(sha256((up/n).read_bytes()).hexdigest()==h for n,h in pins.items())
sys.argv=[str(up/'verify.py'),str(base),str(base.parent/'dim73-82-gamma72-depth6-lattice')]
R=runpy.run_path(sys.argv[0])
S=R['S']; prior=R['bound']
assert prior==6228431712
U=[38860327224660,22886437556,22886438880]
assert [u for u in U if u not in R['p']]==U
assert all(S.get(u)==6 for u in U)
assert all(S[u^v] in (4,6) for u in U for v in R['p'])
assert [S[a^b]for a,b in combinations(U,2)]==[6,4,4]
chart=2; q=R['raw'][chart][U[-1]]
d=[R['raw'][chart][u^U[-1]] for u in U[:2]]
assert all(x[2]==0 and (x[0]or x[1]) and (not x[0]or not x[1]or x[0]==x[1])for x in d)
V=R['make_triple'](q);V2=R['two_occupied'](q)
assert 8*len(V)+4*len(V2)==136704
X=np.ones(len(V),dtype='bool');Y=np.ones(len(V2),dtype='bool')
for r in d:
 X=R['cutgood'](V,r,X);Y=R['cutgood'](V2,r,Y)
clip=8*int(X.sum())+4*int(Y.sum())
assert clip==104704
def inner(x,y):
 return 2*x[0]*y[0]+2*x[1]*y[1]-x[0]*y[1]-x[1]*y[0]
A=[(1,0),(0,1),(1,1)]
assert all(inner(a,a)==2 for a in A) and all(abs(inner(a,b))==1 for a,b in combinations(A,2))
expected={
 '81':{'prior':6228379802,'lower_bound':6228705122,'gain':325320,'pure_sidecar':2,'sidecar_counts_per_sign':[136704]},
 '82':{'prior':6228927150,'lower_bound':6229187942,'gain':260792,'pure_sidecar':6,'sidecar_counts_per_sign':[136704,136704,104704]}}
for a,r in expected.items():
 b=prior+r['pure_sidecar']+2*sum(r['sidecar_counts_per_sign'])
 assert b==r['lower_bound'] and b-r['prior']==r['gain']>0
 print('PASS',a,b,'gain',r['gain'])
assert expected==json.loads((here/'results.json').read_text())

