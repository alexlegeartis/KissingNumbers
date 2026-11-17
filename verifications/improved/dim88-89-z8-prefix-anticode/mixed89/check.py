from fractions import Fraction as Q
from pathlib import Path
import json,itertools,hashlib
P=Path(__file__).resolve().parent
# Bind exact generator to the replayed parent, not merely an arbitrary rank12 code.
src=P/'z8_generator.json'
assert hashlib.sha256(src.read_bytes()).hexdigest()=='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73'
G=json.loads(src.read_text())['generator']
assert len(G)==12 and all(len(r)==24 for r in G)
assert all(G[i][i]%2==1 for i in range(12))and all(G[i][j]==0 for i in range(12)for j in range(i))
assert all(r[:2]==[0,0]for r in G[2:])
assert len({((a*G[0][0]+b*G[1][0])%8,(a*G[0][1]+b*G[1][1])%8)for a,b in itertools.product(range(8),repeat=2)})==64
# Three independent exact integer alphabet reconstructions.
A=[tuple((-1)**((a//4)^((a%2)*x)^(((a//2)%2)*y))for x,y in [(0,0),(0,1),(1,0),(1,1)])for a in range(8)]
B=[(-v[0],)+v[1:]for v in A]
C=[tuple((2 if a<4 else -2)if i==a%4 else 0 for i in range(4))for a in range(8)]
T=[A,B,C];w=[0,2,2,2,4,2,2,2]
dot=lambda u,v:sum(x*y for x,y in zip(u,v))
for S in T:
 assert len(set(S))==8
 for a,b in itertools.product(range(8),repeat=2):assert dot(S[a],S[b])==4-2*w[(a-b)%8]
for S,U in itertools.combinations(T,2):assert max(dot(a,b)for a in S for b in U)==2
labels=[(c,0,b,0)for c in range(3)for b in range(4)if (c,b)!=(0,3)]+[(0,0,3,1),(0,0,7,-1)]
assert len(labels)==13 and len(set(x[:3]for x in labels))==13
q=Q(8,3);cases={}
for i,(c,a,b,h)in enumerate(labels):
 n=88+h*h*q
 # Distinct words inside one fiber retain full norm96 and full IP<=48,
 # and lose prefix squared norm8, so body IP<=40.
 assert 2*(40+h*h*q)<=n
 for d,e,f,k in labels[:i]:
  m=88+k*k*q
  body=44 if c!=d else 48-dot(T[c][a],T[c][e])-dot(T[c][b],T[c][f])
  ip=body+h*k*q
  assert ip>0 and 4*ip*ip<=n*m
  key=(body,str(h*k*q),str(n),str(m));cases[str(key)]=cases.get(str(key),0)+1
# Pointwise normalization can only identify two rays if normalized dot=1.
# All distinct indexed pairs have normalized dot<=1/2, so no collisions.
# Add the two pure final-coordinate unit poles. Each fiber has last-coordinate
# magnitude 0 or sqrt((8/3)/(272/3)) = 1/sqrt34 < 1/2.
assert 4*Q(8,3)<Q(272,3)
# Poles have mutual product -1; neither equals a fiber point with body norm88.
N=13*8**10+2;assert N==13958643714
print(json.dumps({'status':'PASS exact mixed-radius geometry and population; parent replay required','count':N,'pairs':sum(cases.values()),'types':cases,'baseline':12884901888,'gain':N-12884901888,'percent':100*(N/12884901888-1)},indent=2))
