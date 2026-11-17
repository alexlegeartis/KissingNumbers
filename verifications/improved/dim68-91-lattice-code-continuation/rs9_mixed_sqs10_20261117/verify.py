from itertools import combinations,product
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent
def add(x,y):return ((x%3+y%3)%3)+3*((x//3+y//3)%3)
def neg(x):return (-x%3)+3*((-(x//3))%3)
def mul(x,y):
 a,b=x%3,x//3;c,d=y%3,y//3
 return (a*c+2*b*d)%3+3*((a*d+b*c)%3)
def div(x,y):return mul(x,next(z for z in range(1,9)if mul(z,y)==1))
blocks=set()
for a,b,c,d in product(range(9),repeat=4):
 if add(mul(a,d),neg(mul(b,c)))==0:continue
 block=[]
 for x in [9,0,1,2]:
  if x==9:y=9 if c==0 else div(a,c)
  else:
   den=add(mul(c,x),d);y=9 if den==0 else div(add(mul(a,x),b),den)
  block.append(y)
 blocks.add(tuple(sorted(block)))
quads=sorted(blocks)
assert len(quads)==30 and all(len(set(a)&set(b))<=2 for a,b in combinations(quads,2))
assert sorted(t for q in quads for t in combinations(q,3))==list(combinations(range(10),3))
assert sum(9 in q for q in quads)==12
rounds=[sorted(tuple(sorted(((c+r)%9,(c-r)%9)))for r in range(1,5))for c in range(9)]
assert len(set(p for R in rounds for p in R))==36
assert all(len(set(a)&set(b))==0 for R in rounds for a,b in combinations(R,2))
words=[(a,b,c,a^b^c)for a,b,c in product(range(4),repeat=3)]
assert len(words)==64 and min(sum(a!=b for a,b in zip(x,y))for x,y in combinations(words,2))==2
signs=[x for x in range(256)if x.bit_count()%2==0]
assert len(signs)==128 and min((a^b).bit_count()for a,b in combinations(signs,2))==2
rows=[]
for n,M in [(87,9),(88,12),(89,16),(90,27),(91,27)]:
 k=min(9,n-81);supports=[]
 for Q in quads:
  for R in rounds:
   for word in words:
    sup=tuple(sorted(9*g+x for g,i in zip(Q,word)for x in R[i]))
    if max(sup)<81+k:supports.append(sup)
 assert len(supports)==len(set(supports))==18*576+12*16*(k*(k-1)//2)
 previous=2*M*3**18+72*k*9**5+2*n*(n-1);extra=128*len(supports)
 rows.append(dict(dimension=n,support8_count=len(supports),extra_points=extra,previous_total=previous,total=previous+extra,supports_sha256=hashlib.sha256(json.dumps(sorted(supports)).encode()).hexdigest()))
json.dump({'group_quads':quads,'near_perfect_matchings':rounds,'pair_label_code':words,'sign_masks':signs,'rows':rows},open(P/'certificate.json','w'),indent=2)
print(json.dumps(rows,indent=2))
