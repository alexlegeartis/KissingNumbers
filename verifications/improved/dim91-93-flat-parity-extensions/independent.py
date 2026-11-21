"""Independent exact audit against pinned explicit pattern and dense head inputs."""
from functools import reduce
from itertools import combinations, product
from operator import xor
from pathlib import Path
import json
root=Path(__file__).parent
def get(k):
 if k in ('patterns.json','weight10.json'): q=root.parent/'dim89-two-color-latitudes'/'abbreviated_parent'/k
 elif k=='patterns_extra.json': q=root.parent/'dim88-89-golay-variable-punctures'/k
 elif k in ('heads91.json','heads92.json','heads93.json'):
  d=int(k[5:7]); folder='dim91-z8-asymmetric-planes' if d==91 else f'dim{d}-z8-variable-radius-heads-continuation'
  q=root.parent/folder/'heads.json'
 else: q=root/k
 return json.loads(q.read_text())
w=int.bit_count
F=get('flat_supports.json')
add=[(s,p,k) for name,k in [('six',6),('seven',7)] for s,p in F[name]]
assert len(F['six'])==114 and len(F['seven'])==131 and len({s for s,p,k in add})==245
G=[0];coef=(1,1,0,0,0,1,1,1,0,1,0,1)
for i in range(12):
 x=sum(b<<(i+j) for j,b in enumerate(coef))|1<<23;G.extend(y^x for y in G[:])
assert len(set(G))==4096
p=get('patterns.json')[0]
assert p['mask']==p['pivot']==3 and len(p['ids'])==629
P=[(G[i&4095]|G[i>>12])&~8 for i in p['ids']]
P += [a|b for a,b in get('weight10.json')+get('patterns_extra.json')]
assert len(P)==993 and sum(w(p)==10 for p in P)==11 and sum(w(p)==11 for p in P)==982
# Check all four possible relative choices directly, not via min mismatch reduction.
for i,(s,p,n) in enumerate(add):
 assert w(s)==n and p&s==p and w(p)==n//2 and s<(1<<22)
 for q in P:
  assert 4*(4*w((s<<2)&q))**2 <= 4*n*8*w(q)
 for t,q,m in add[:i]:
  Z=s&t
  for p1,p2 in product((0,1),repeat=2):
   d=w(Z&((p^(s*p1))^(q^(t*p2))))
   cap=4*w(Z)-2*d
   assert 4*cap**2 <= 4*n*4*m
# Reproduce the two additive axis codes via direct coordinate recurrences.
def tr(x):return (x<<1)^((x>>1)*7)
A6=[(a,b,c,a^b^c,a^tr(b)^tr(tr(c)),a^tr(tr(b))^tr(c)) for a,b,c in product(range(4),repeat=3)]
A7=[]
for u in range(16):
 for v in range(8):
  sign=sum((w((u>>1)&i)%2^(u&1))<<(i-1) for i in range(1,8))
  sh=reduce(xor,(g for i,g in enumerate([15,0,8,1]) if u>>i&1),0)
  hi=sum((w(v&i)%2)<<j for j,i in enumerate([7,3,4,6,2,1,5]))^sh
  A7.append(tuple((sign>>i&1)+2*(hi>>i&1) for i in range(7)))
for n,A in [(6,A6),(7,A7)]:
 assert len(A)==len(set(A))==2**n
 assert min(sum(x!=y for x,y in zip(a,b))for a,b in combinations(A,2))==4
 assert 4*(n-4)<=2*n and 4*n-4*4<=2*n
 assert 4*(4*((n+1)//2)+2*(n//2))**2<=4*n*88
# Independent direct signed-squaring of dense head pairs from frozen inputs.
def bd(a,b): return 4 if a==b else -4 if (a^b)==4 else 0
k=[24,32,44]
for d,m in zip([91,92,93],k):
 x=get(f'heads{d}.json');D=2 if d<93 else 1;D*=D
 assert len(x)==len({(r['color'],*r['prefix'])for r in x})==m
 for u in x:
  assert len(u['head'])==d-88 and sum(t*t for t in u['head'])<=8*D
 for u,v in combinations(x,2):
  body=44 if u['color']!=v['color'] else 48-sum(bd(a,b)for a,b in zip(u['prefix'],v['prefix']))
  ip=D*body+sum(a*b for a,b in zip(u['head'],v['head']))
  NS=[88*D+sum(z*z for z in q['head'])for q in (u,v)]
  assert ip<=0 or 4*ip**2<=NS[0]*NS[1]
 rhs=m*8**10+982*1024*144+11*256*72+114*2*64*32+131*2*128*64
 print(d,rhs,rhs-(m*8**10+982*1024*144+11*256*72+79*4096))
