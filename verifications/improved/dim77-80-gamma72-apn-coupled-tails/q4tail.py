"""Exact E8 norm-four quotient classes and best singular-plane subsets."""
from itertools import combinations,product
from collections import defaultdict
import numpy as np
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def key(v):
 r=tuple(x%4 for x in v);s=tuple((x+2)%4 for x in v)
 return (min(r,s),sum(v)%8)
def coords(n):
 if n==2:
  for ij in combinations(range(8),2):
   for signs in product((-2,2),repeat=2):
    v=[0]*8
    for i,t in zip(ij,signs):v[i]=t
    yield tuple(v)
  for v in product((-1,1),repeat=8):
   if sum(v)%4==0:yield v
 elif n==4:
  for i in range(8):
   for t in (-4,4):yield tuple(t if j==i else 0 for j in range(8))
  for ix in combinations(range(8),4):
   for s in product((-2,2),repeat=4):
    v=[0]*8
    for i,t in zip(ix,s):v[i]=t
    yield tuple(v)
  for i in range(8):
   for s in product((-1,1),repeat=8):
    v=tuple(t*(3 if j==i else 1)for j,t in enumerate(s))
    if sum(v)%4==0:yield v

def groups():
 W=defaultdict(list);q={key((0,)*8):0};R={key((0,)*8):(0,)*8}
 for z in coords(2):
  s=key(z); R[s]=z;q[s]=1
 for z in coords(4):
  s=key(z); assert s not in q or q[s]==0
  q[s]=0; R[s]=z;W[s].append(z)
 assert len(R)==256 and len(W)==135 and all(len(z)==16 for z in W.values())
 L=sorted(R);I={x:i for i,x in enumerate(L)}; n=I[key((0,)*8)]
 A=np.array([R[s]for s in L],dtype='i8')
 assert not np.any((A@A.T)%4)
 B=((A@A.T)//4)%2
 V={frozenset([n])};cand=[i for i,s in enumerate(L) if s in W]
 for k in range(4):
  Z=set()
  for s in V:
   for i in cand:
    if i not in s and all(B[j,i]==0 for j in s):
     Z.add(s|frozenset(I[key(add(R[L[j]],R[L[i]]))]for j in s))
  V=Z
 assert len(V)==270 and all(len(s)==16 for s in V)
 for m in (5,6,7,8):
  stock={}
  for i,s in enumerate(L):
   stock[i]=W.get(s,())
   if m==5:stock[i]=[u for u in stock[i] if not any(u[5:])]
   elif m<8:stock[i]=[u for u in stock[i] if len(set(u[m-1:]))==1]
  _,ss=max((sum(len(stock[i])for i in sv),tuple(sorted(sv)))for sv in V)
  rows=[(u,j)for j,i in enumerate(ss,1)for u in stock[i]]
  assert all(sum(x*x for x in u)==16 for u,j in rows)
  assert all(sum(a*b for a,b in zip(x,y))<=(0 if j==i else 8)for (x,j),(y,i)in combinations(rows,2))
  yield m,rows
if __name__=='__main__':
 for m,u in groups():print(m,len(u),len({t for v,t in u}))
