"""Independently generate the 24, 32, 44 Z8 head assignments in ranks 3, 4, 5."""
from itertools import combinations, product
def profiles(k):
 # Each pair is a specialized singleton to receive four extra small heads.
 special={3:[(0,1),(0,2),(1,2)],4:[(0,1),(2,3)],5:[(0,1)]}[k]
 all_edges=set(combinations(range(k),2))
 result=[]
 def v(i,j,a,b):
  s=[0]*k
  s[i]=a;s[j]=b
  return tuple(s)
 for c,(i,j) in enumerate(special):
  for t,(a,b) in enumerate(product((-1,1),repeat=2)):
   result.append((c,0,t,v(i,j,4*a,4*b)))
   result.append((c,0,t+4,v(i,j,-a,-b)))
 rest=all_edges-set(special)
 if k==4:
  classes=[((0,1),(2,3))]
 elif k==5:
  # Prefix code with first two axes sharing a 2-bit label. Take first differing bit.
  key=[0,0,1,2,3];classes=[]
  for bit in (1,0):
   left=tuple(i for i in range(k)if not(key[i]>>bit&1))
   right=tuple(i for i in range(k)if key[i]>>bit&1)
   classes.append((left,right))
 else:classes=[]
 for c,(L,R) in enumerate(classes,len(special)):
  assert max(len(L),len(R))<=4
  take={(i,j) for i in L for j in R if tuple(sorted((i,j))) in rest}
  rest.difference_update(tuple(sorted(t)) for t in take)
  for i,j in take:
   for a,b in product((-1,1),repeat=2):
    result.append((c,L.index(i)+4*(a<0),R.index(j)+4*(b<0),v(i,j,4*a,4*b)))
 assert not rest
 return result
def test(k):
 H=profiles(k); n={3:24,4:32,5:44}[k]
 assert len(H)==len({x[:3] for x in H})==n
 norm=lambda h:352+sum(x*x for x in h)
 ip=lambda x,y:0 if x%4!=y%4 else (4 if x==y else -4)
 for *_,h in H:
  assert sum(x*x for x in h) in (2,32)
  assert 2*(160+sum(x*x for x in h))<=norm(h)
 for (c,a,b,h),(e,f,g,t) in combinations(H,2):
  body=44 if c!=e else 48-ip(a,f)-ip(b,g)
  dd=4*body+sum(x*y for x,y in zip(h,t))
  assert dd<=0 or 4*dd*dd<=norm(h)*norm(t)
 print(k+88,len(H),len(H)*(len(H)-1)//2)
if __name__=='__main__':
 for k in (3,4,5):test(k)
