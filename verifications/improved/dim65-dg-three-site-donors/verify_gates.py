from itertools import product
D=200
P=[(-1129,250,250,0),(280,-1100,256,214),(280,256,-1100,-214)]
N=[61*D**2+sum(a*a for a in p)for p in P]
def chk(m,a,b): assert m<=0 or 4*m*m<=a*b
for i,p in enumerate(P):
 ni=N[i]
 chk((32*D+sum(abs(D-x)for x in p[:3]))*D, ni, 64*D**2)
 chk(abs(p[-1])*D,ni,D**2)
 # Safe for any support that uses unit entries, including 65th mark.
 chk(sum(sorted([D]*61+[abs(x)for x in p],reverse=True)[:4])*D,ni,4*D**2)
 for cap in (24*D+sum(max(0,x-D) for x in p[:3])+2*abs(p[-1]),
             20*D+sum(abs(x-D) for x in p[:3])+2*abs(p[-1])):
  assert 8*cap*cap<=72*ni
 # Fixed zero-tail old norm16 sign banks; guard enforces opposite parity if and only if
 # spike is inside. For its smallest-magnitude site, magnitude is at least one.
 for m in range(8):
  t=16*D+sum(abs(p[h])-D for h in range(3)if m>>h&1)
  if any(x<0 and m>>h&1 for h,x in enumerate(p[:3])): t-=2*D
  chk(t*D,ni,16*D**2)
 for j,q in enumerate(P):
  n_j=N[j]
  if i!=j: chk(61*D*D+sum(x*y for x,y in zip(p,q)),ni,n_j)
  # arbitrary sign agreement on distinguished coordinates, outside capped by coarse DG dot
  for s in product((-1,1),repeat=3):
   chk(16*D*D+sum(e*(u*v-D*D) for e,u,v in zip(s,p[:3],q[:3]))+p[-1]*q[-1],ni,n_j)
  # Tails only overlap between different marked triples; ordinary body baseline is ≤32.
  chk(32*D*D+(sum(abs(D-x)for x in p[:3])+sum(abs(D-y)for y in q[:3]))*D+p[-1]*q[-1],ni,n_j)
base=335544320+40414*2048+10741*16+2
assert base==418484050
assert base+180224+21*131072*2==424169298
print("independently checked all 65D profile pair/old-bank gates, numerators",N,"count 424169298")
