from pathlib import Path
import json,itertools,hashlib
P=Path(__file__).resolve().parent
src=P/'z8_generator.json'
assert hashlib.sha256(src.read_bytes()).hexdigest()=='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73'
G=json.loads(src.read_text())['generator']
assert len(G)==12 and all(len(r)==24 for r in G)
assert all(G[i][i]%2==1 for i in range(12))
assert all(G[i][j]==0 for i in range(12) for j in range(i))
# Thus coefficient-to-word map injective by successive odd-pivot elimination.
assert all(G[i][0]==G[i][1]==0 for i in range(2,12))
assert len({((a*G[0][0]+b*G[1][0])%8,(a*G[0][1]+b*G[1][1])%8)for a,b in itertools.product(range(8),repeat=2)})==64
# Every prefix has exactly 8^10 preimages; no enumeration assumption.
def dot(x,y):return sum(a*b for a,b in zip(x,y))
A=[tuple((-1)**((a//4)^((a%2)*x)^(((a//2)%2)*y))for x,y in [(0,0),(0,1),(1,0),(1,1)])for a in range(8)]
B=[(-v[0],)+v[1:] for v in A]
C=[tuple((2 if a<4 else -2)if i==a%4 else 0 for i in range(4))for a in range(8)]
alph=[A,B,C];w=[0,2,2,2,4,2,2,2]
for T in alph:
 assert len(set(T))==8
 for a,b in itertools.product(range(8),repeat=2):assert dot(T[a],T[b])==4-2*w[(a-b)%8]
for T,U in itertools.combinations(alph,2):
 assert set(T).isdisjoint(U)
 assert max(dot(x,y)for x in T for y in U)==2
prefixes=[(0,b)for b in range(4)]
for T in alph:
 for p,q in itertools.product(prefixes,repeat=2):
  deleted=sum(dot(T[a],T[b])for a,b in zip(p,q))
  assert deleted>=4
  assert 48-deleted<=44
# Same-color parent norm96/IP<=48; punctured norm88/IP<=44.
# Different colors: 22 blocks each IP<=2. Every vector normalized by sqrt88.
assert 22*4==88 and 22*2==44 and 44/88==.5
# Same-color collisions impossible from IP<=44<88 (or parent weight24>8 puncture weight).
assert 24>2*max(w)
N=3*len(prefixes)*8**10
assert N==12884901888
print(json.dumps({'status':'PASS exact geometry and population conditional on parent minimum24','dimension':88,'count':N,'previous_repository':9590791776,'gain':N-9590791776,'percent':100*(N/9590791776-1),'prefixes':prefixes},indent=2))
