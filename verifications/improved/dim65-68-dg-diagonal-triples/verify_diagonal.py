from itertools import product,combinations
from fractions import Fraction as F
D=100000
P=[(-558229,105952,154047,40547,40547),
   (160809,-557420,100000,40547,40547),
   (103127,157232,-557870,-40547,-40547)]
N=[61*D*D+sum(a*a for a in p)for p in P]
V=[sum(abs(a-D)for a in p[:3])for p in P]
checks=[]
def gate(dot,n,m,label):
 good=dot<=0 or 4*dot*dot<=n*m
 checks.append((good,label))
 assert good,(F(dot,D*D),n,m,label)
 return good
for i,p in enumerate(P):
 n=N[i]; assert n<=96*D*D
 gate((32*D+V[i])*D,n,64*D*D,('ordinary',i))
 gate(sum(sorted([D]*61+list(map(abs,p)),reverse=True)[:4])*D,n,4*D*D,('fours',i))
 gate(max(D,*(abs(x)for x in p[-2:]))*D,n,D*D,('tail and ordinary poles',i))
 for val in (24*D+sum(max(0,x-D)for x in p[:3])+2*abs(p[3]),20*D+V[i]+2*abs(p[3])):
  assert 8*val*val<=72*n
 for j,q in enumerate(P):
  m=N[j]
  if i!=j:gate(61*D*D+sum(a*b for a,b in zip(p,q)),n,m,('same owner',i,j))
  for s in product((-1,1),repeat=3):
   gate(16*D*D+sum(a*(b*c-D*D) for a,b,c in zip(s,p[:3],q[:3]))+p[-2]*q[-2]+p[-1]*q[-1],n,m,('same triple',i,j))
  gate(32*D*D+(V[i]+V[j])*D+p[-2]*q[-2]+p[-1]*q[-1],n,m,('different triple',i,j))
 for mask in range(8):
  s=sum(abs(x)-D for j,x in enumerate(p[:3])if mask>>j&1)
  gate((16*D+s-2*D*(mask>>i&1))*D,n,16*D*D,('zero-tail',i,mask))
  gate((12*D+s+2*abs(p[3]))*D,n,16*D*D,('four-tail',i,mask))
  for overlap in (0,1):gate((14*D+s+overlap*abs(p[3]))*D,n,16*D*D,('short overlap',i,mask))
  unsafe=4*((14*D+s+2*abs(p[3]))*D)**2>n*16*D*D
  if unsafe:
   assert mask>>i&1
   gate((14*D+s)*D,n,16*D*D,('forced mismatch',i,mask))
   print('unsafe',i,mask)
print('PASS diagonal donor geometry for every support and profile type.')
