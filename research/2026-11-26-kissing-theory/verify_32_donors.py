"""Check the 32-donor signed tetrahedral module on the four upstream charts."""
import hashlib,itertools,json,sys
from pathlib import Path
raw=Path(sys.argv[1]).read_bytes()
assert hashlib.sha256(raw).hexdigest()=="5018425664146ee2ace418021537b8b13863fc57c7dc64376ab6f9e811aa5162"
d=json.loads(raw)
B=d["five_coordinate_supports"]
Q=d["four_coordinate_supports"]
O=[(0,3,1,2,9,8),(0,3,5,7,8,9),(4,6,1,2,8,9),(4,6,5,7,9,8)]
def sign(mask,i): return -1 if mask>>i&1 else 1
G=[]
for o in O:
    rows=[]
    for u,s in itertools.product((-7,7),(0,3,5,6)):
        x=[0]*10
        for i,v in zip(o,(u,5*sign(s,0),6*sign(s,1),6*sign(s,2),7,1)): x[i]=v
        assert sum(v*v for v in x)==196 and max(map(abs,x))==7
        assert sum(abs(x[i]) for i in o[:4])+2*abs(x[o[5]])==26
        rows.append(x)
    G.append(rows)
F=sum(G,[])
assert len({tuple(x) for x in F})==32
assert max(sum(a*b for a,b in zip(x,y)) for x,y in itertools.combinations(F,2))==98
paid=[(i,s) for i,q in enumerate(Q) for s in range(16)
      if any(sum(x[j]*sign(s,k) for k,j in enumerate(q))>14 for x in F)]
assert paid==[(5,s) for s in range(4)]
assert max(sum(abs(x[k]) for k in b) for i,m in enumerate(G) for x in m for j,b in enumerate(B) if i!=j)==13<14
assert max(len(set(a)&set(b)) for a,b in itertools.combinations(B,2))<=2
assert max(len(set(q)&set(b)) for q in Q for b in B)<=2
assert 200864+20+224-len(paid)+len(F)+4*19==201212
print("PASS: exact 32 site module; analytic Q4 chamber obstruction is in the note.")
