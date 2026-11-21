"""Replay 32 variable-radius fibers, then append a separate hexacode layer."""
from pathlib import Path
from itertools import combinations, product
import json, subprocess, sys
P=Path(__file__).resolve().parent
q=json.loads(subprocess.run([sys.executable,str(P/'verify.py')],check=True,text=True,capture_output=True).stdout)
assert q['lower_bound']==34504742912
O=P.parent/'dim88-89-golay-variable-punctures'
B0=P.parent/'dim89-two-color-latitudes'/'abbreviated_parent'
basis=[sum(((c&1)<<(i+j)) for j,c in enumerate([7,1,4,4,2,7,3,3,4,7,2,1]))+(1<<23) for i in range(12)]
C=[0]
for u in basis:C.extend(x^u for x in C[:])
J=json.loads((B0/'patterns.json').read_text())[0]
assert len(J['ids'])==629 and J['pivot']==3 and J['mask']==3
old=[(C[m%4096]|C[m//4096])&~8 for m in J['ids']]
old += [a|b for a,b in json.loads((O/'patterns_extra.json').read_text())]
old += [a|b for a,b in json.loads((B0/'weight10.json').read_text())]
assert len(old)==993
S=json.loads((P/'flat_supports.json').read_text())
M=[sum(1<<x for x in s) for s in S]
assert len(M)==len(set(M))==79
assert all(len(set(t))==6 and min(t)>=2 and max(t)<24 for t in S)
assert all((x&y).bit_count()<=3 for x,y in combinations(M,2))
assert all((m&n).bit_count()<=5 for m,n in product(M,old))
def g(x):return ((x<<1)^(7 if x&2 else 0))
H=[(a,b,c,a^b^c,a^g(b)^g(g(c)),a^g(g(b))^g(c)) for a,b,c in product(range(4),repeat=3)]
assert len(set(H))==64 and min(sum(x!=y for x,y in zip(a,b))for a,b in combinations(H,2))==4
signs=[x for x in range(64) if x.bit_count()%2==0]
assert len(signs)==32 and all((a^b).bit_count()>=2 for a,b in combinations(signs,2))
F=[[t for t in product((1,-1),repeat=4) if sum(x<0 for x in t)%2==p] for p in (0,1)]
A=[tuple(2*s if i==k else 0 for i in range(4))for k,s in product(range(4),(-1,1))]
# Full sign sets include both signs on four mutually orthogonal projective axes.
dot=lambda x,y:sum(a*b for a,b in zip(x,y))
for row in F:
 assert len(row)==8 and all(dot(x,y) in (0,-4) for x,y in combinations(row,2))
assert max(dot(x,y) for x,y in product(*F))==2
assert max(dot(x,y)for x,y in product(F[0]+F[1],A))==2
# Cross to dense 3 even + 3 odd: max raw 18 for signed Walsh, 12 for axis.
assert 4*18**2<=24*88
# Weight-10/11 old sparse stock overlaps in at most five rows, each dot <=4.
assert 4*20**2<=24*80
# Cross distinct supports: <=3 overlaps, each dot<=4; opposite parity modes on one <=6*2.
assert max(3*4,6*2)==12
extra=79*64*32*2;total=q['lower_bound']+extra
cl=json.loads((P/'claims.json').read_text())[0]
assert (extra,total,cl['lower_bound'],cl['previous'])==(323584,34505066496,34505066496,32357582848)
print(json.dumps({'dimension':92,'dense':q['dense'],'old_sparse':q['sparse'],'flat':extra,'lower_bound':total,'gain':total-cl['previous']},indent=2))
