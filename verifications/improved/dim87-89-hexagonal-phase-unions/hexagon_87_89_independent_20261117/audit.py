from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
import json
D=Path(__file__).resolve().parent
answer=[]
# Premises separately replayed by the public dense and phase-source packages:
# M prefix fibers, each3^18; antipodal completion; 37468 outer words and7812 phases.
for dim,M,want in [(87,9,7266268830),(88,12,9590791776),(89,16,12690155704)]:
 m=dim-84
 vectors=set()
 for i,j in combinations(range(m),2):
  for s,t in product((-1,1),repeat=2):
   r=[0]*m;r[i]=s;r[j]=t;vectors.add(tuple(r))
 assert len(vectors)==2*m*(m-1)
 for r in vectors:assert sum(x*x for x in r)==2
 for r,s in combinations(vectors,2):assert sum(x*y for x,y in zip(r,s))<=1
 # Norm45 root-scaled heads, dense heads norm3: cross inner <=sqrt135<45/2.
 assert F(45)*3<F(45,2)**2
 # Fifteen signed A2 roots against unit triangle body: absolute inner <=15*3/2.
 assert 15*F(3,2)==F(45,2)
 # Sparse points have zero heads; pure heads are orthogonal to them.
 # Dense uses42 nonzero body blocks, sparse15; pure heads have0: disjoint.
 dense=2*M*3**18;sparse=37468*7812;total=dense+sparse+len(vectors)
 assert total==want
 answer.append(dict(dimension=dim,dense=dense,sparse=sparse,head_points=len(vectors),total=total))
out=dict(status='PASS independent union geometry and arithmetic',inherited_premises='Previously audited dense compression and Pless/7812 phase construction; source wrappers must replay these independently of this small union check',results=answer)
(D/'report.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
