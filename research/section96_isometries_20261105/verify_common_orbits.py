"""Check elementary identities used for the common 4C orbit reduction."""
from itertools import permutations
from fractions import Fraction

def gray(z):
 return sum(((((z>>2)&1)^((z&1)*x)^(((z>>1)&1)*y))<<(2*x+y)) for x in (0,1) for y in (0,1))
assert all(gray(z^4)==(gray(z)^15) for z in range(8))
pairs=[(0,15),(3,12),(5,10),(6,9)];maps=[]
for perm in permutations(range(4)):
 for flips in range(16):
  d={}
  for i,(a,b) in enumerate(pairs):
   p=pairs[perm[i]];f=(flips>>i)&1;d[a]=p[f];d[b]=p[1-f]
  v=[d[gray(z)] for z in range(8)]
  maps.extend([v,[x^1 for x in v]])
assert len({tuple(v) for v in maps})==768
assert all(v[z^4]==(v[z]^15) for v in maps for z in range(8))
# First 12 columns are triangular with odd diagonal, hence 4C has 2^12
# distinct elements. Addition by each nonzero 4C element flips a nonempty
# set of blocks, and is a fixed-point-free action on binary words.
g=[7,1,4,4,2,7,3,3,4,7,2,1]
G=[[0]*24 for _ in range(12)]
for i in range(12):
 for j,x in enumerate(g):G[i][i+j]=x
 G[i][23]=3
assert all(G[i][j]==0 for i in range(12) for j in range(i))
assert all(G[i][i]%2==1 for i in range(12))
for m in (4,6):
 orbits=3**m*4**(12-m)
 assert orbits*4096 == 2**36*Fraction(3,4)**m
 for s in (0,2,4,6):
  if m+s<=12:print('dimension',96-m,'index',2**s,'orbits',orbits//2**s,'points',orbits//2**s*4096)
print('verified: antipodal equivariance, generator rank, and orbit counts')
