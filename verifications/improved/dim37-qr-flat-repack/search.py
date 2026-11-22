from itertools import combinations,product
def eqsolve(eq,N):
 P={}
 for a,w in eq:
  while a:
   i=a.bit_length()-1
   if i in P:a^=P[i][0];w^=P[i][1]
   else:P[i]=(a,w);break
  else:
   if w:return None
 x=0
 for i,(a,w) in sorted(P.items()):x|=(w^((x&a).bit_count()&1))<<i
 return x
def direct(G,E):
 masks=[63^((1<<i)|(1<<j)) for i,j in E]
 eq=[]
 for i,j in combinations(range(len(E)),2):
  a=masks[i]&masks[j]
  if a.bit_count()!=3:continue
  h=next(h for h in range(1,64) if h&a==h and all((h&g).bit_count()%2==0 for g in G))
  eq.append(((h<<(6*i))|(h<<(6*j)),1))
 x=eqsolve(eq,len(E)*6)
 if x is None:return
 z=[(x>>(i*6))&masks[i] for i in range(len(E))]
 return masks,z
def phys(m):return sum(3<<(2*i) for i in range(6) if m>>i&1)
def sig(m):return sum(1<<(2*i) for i in range(6) if m>>i&1)
def run():
 for gg in product(range(3),repeat=6):
  if tuple(gg.count(i) for i in range(3))!=(2,2,2):continue
  G=[sum(1<<i for i in range(6) if gg[i]!=j) for j in (0,1)]
  E=list(product(range(3),range(3,6)))
  for k in [9]:
   for q in combinations(E,k):
    s=direct(G,q)
    if s:
     print(gg,G,s)
     with open('pairwin.txt','w') as f:
      f.write(' '.join(map(str,[3<<(2*i) for i in range(6)]+[sig(g) for g in G]))+'\n')
      for m,x in zip(*s):f.write(f'{phys(m)} {sig(x)} 64\n')
     return
run()
