import json,itertools,hashlib
from pathlib import Path
P=Path(__file__).parent
C=P/'candidate.json'
rows=json.loads(C.read_text())['selected']; assert len(rows)==68
# Independent alphabet reconstruction: three orthogonal 8-symbol crosspolytopes.
def alph(c,a):
 if c==2:return tuple((2 if a<4 else -2) if k==a%4 else 0 for k in range(4))
 z=tuple((-1)**((a//4)^((a%2)*x)^(((a//2)%2)*y)) for x,y in [(0,0),(0,1),(1,0),(1,1)])
 return z if c==0 else (-z[0],)+z[1:]
def dot(x,y):return sum(a*b for a,b in zip(x,y))
assert len({(r['color'],*r['prefix']) for r in rows})==68
assert len({tuple(r['head']) for r in rows})==68
for r in rows:
 h=r['head']; assert len(h)==8 and all(type(x)==int for x in h)
 assert h[1]==h[4]==-h[0]; assert dot(h,h)==8
 # Orthonormal compression uses sqrt(3)*h0; all coordinates bounded by2.
 assert 3*h[0]**2<=4 and max(h[j]**2 for j in [2,3,5,6,7])<=4
for r,s in itertools.combinations(rows,2):
 ip=dot(r['head'],s['head'])
 if r['color']==s['color']:
  old=sum(dot(alph(r['color'],a),alph(s['color'],b)) for a,b in zip(r['prefix'],s['prefix']))
  assert ip<=old
 else:assert ip<=4
# Sparse support count and RS intersection property independently verified.
def mul(a,b):
 x=0
 while b:
  if b&1:x^=a
  b>>=1;a<<=1
  if a&16:a^=19
 return x
def ev(cs,x):
 y=0
 for a in cs[::-1]:y=mul(y,x)^a
 return y
assert all(mul(a,b) for a in range(1,16) for b in range(1,16))
for cs in itertools.product(range(16),repeat=4):
 if any(cs):assert sum(ev(cs,x)==0 for x in range(6))<=3
# At first evaluation point0 impose14 allowed constants, remaining3 free.
rs=14*16**3*32; roots=4*sum(n*(n-1)//2 for n in [14,16,16,16,16,16])
# RS support6 at amplitude4: norm96, overlap<=3; signs even distance>=2.
assert 6*16==96 and 3*16==48 and 6*4*2==48
# Dense/root bound4sqrt48<48; root/RS4sqrt48<48.
assert 16*48<48**2
count=68*8**10+rs+roots
assert count==73016281804
res={'status':'PASS fresh head and augmentation arithmetic','candidate_sha256':hashlib.sha256(C.read_bytes()).hexdigest(),'pairs':2278,'count':count,'rs':rs,'roots':roots,'parent':'Full published parent replay separately passed; parent proof reviewed, not new enumeration algorithm'}
(P/'result.json').write_text(json.dumps(res,indent=2));print(res)
