from pathlib import Path
import json,hashlib,itertools
P=Path(__file__).resolve().parent; S=P/'input'; source=S/'z8_generator.json'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73'
# Public source byte identity was verified before freezing; SHA above binds the input.
G=json.loads(source.read_text())['generator'];assert len(G)==12 and all(len(r)==24 for r in G)
assert all(G[i][i]%2 for i in range(12)) and all(G[i][j]==0 for i in range(12)for j in range(i))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
A=[tuple((-1)**((a//4)^((a%2)*x)^(((a//2)%2)*y))for x,y in [(0,0),(0,1),(1,0),(1,1)])for a in range(8)]
B=[(-v[0],)+v[1:]for v in A];C=[tuple((2 if a<4 else -2)if i==a%4 else 0 for i in range(4))for a in range(8)];bases=[A,B,C]
weight=[0,2,2,2,4,2,2,2]
for b in bases:
 assert len(set(b))==8
 for a,c in itertools.product(range(8),repeat=2):assert dot(b[a],b[c])==4-2*weight[(a-c)%8]
for b,c in itertools.combinations(bases,2):
 assert set(dot(x,y)for x in b for y in c)=={-2,2}
 assert set(b).isdisjoint(c)
# First-two-coordinate map is surjective with exactly one solution for first two coefficients.
assert all(G[i][0]==G[i][1]==0 for i in range(2,12))
assert len({((a*G[0][0]+b*G[1][0])%8,(a*G[0][1]+b*G[1][1])%8)for a,b in itertools.product(range(8),repeat=2)})==64
cert=json.loads((S/'certificate.json').read_text());out={96:3*8**12}
for rec in cert['compressed_results']:
 d=rec['head_dimension'];n=rec['dimension'];assert n==88+d
 labels=rec['labels'];heads=rec['heads'];assert len(labels)==len(heads)==2*d*(d-1)
 expected=set()
 for i,j in itertools.combinations(range(d),2):
  for a,b in itertools.product([-2,2],repeat=2):
   h=[0]*d;h[i]=a;h[j]=b;expected.add(tuple(h))
 assert set(map(tuple,heads))==expected
 assert len(set(map(tuple,labels)))==len(labels)
 for i,(color,a,b) in enumerate(labels):
  assert color in range(3) and a in range(8) and b in range(8)
  for j,(other,c,e) in enumerate(labels):
   ip=dot(heads[i],heads[j])
   if color==other:assert ip==dot(bases[color][a],bases[color][c])+dot(bases[color][b],bases[color][e])
   else:assert ip<=4 and 44+ip<=48
 assert rec['fiber_size']==8**10 and rec['count']==len(heads)*8**10
 out[n]=rec['count']
assert out=={96:206158430208,93:42949672960,94:64424509440,95:90194313216}
print(json.dumps({'status':'PASS elementary exact geometry/counts and public generator identity; inherited parent minimum24','counts':out},indent=2))
(P/'result.json').write_text(json.dumps({'status':'PASS conditional on published parent minimum24 proof','counts':out},indent=2))
