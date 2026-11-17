import itertools,json
from pathlib import Path
P=Path(__file__).resolve().parent
# F9=F3[u]/(u²+1). u²=2, and x²+1 has no F3 root.
assert all((x*x+1)%3 for x in range(3))
def add(x,y):return ((x%3+y%3)%3)+3*((x//3+y//3)%3)
def mul(x,y):
 a,b=x%3,x//3;c,d=y%3,y//3
 return (a*c+2*b*d)%3+3*((a*d+b*c)%3)
for a in range(1,9):assert any(mul(a,b)==1 for b in range(1,9))
# Shortened binary Hamming code. Ten distinct nonzero parity-check columns.
cols=[x for x in range(1,16) if x not in [1,2,4,8,15]]
signs=[]
for mask in range(1024):
 syndrome=0
 for i,h in enumerate(cols):
  if (mask>>i)&1:syndrome^=h
 if syndrome==0:signs.append(mask)
assert len(signs)==64 and 1023 in signs
assert min(m.bit_count() for m in signs if m)>2
assert all((m^1023)in signs for m in signs)
# The polynomial graph has evaluations at all9 field points plus coefficient a5.
# Distinct degree≤5 polynomials agree in≤5 finite positions; if a5 agrees,
# their difference has degree≤4, hence including infinity still≤5 agreements.
rows=[]
for n,M in [(87,9),(88,12),(89,16),(90,27),(91,27)]:
 k=min(9,n-81);supports=k*9**5;sparse=64*supports+2*n*(n-1);dense=2*M*3**18
 rows.append(dict(dimension=n,infinity_alphabet=k,supports=supports,signs_per_support=64,dense=dense,sparse=sparse,total=dense+sparse))
assert 3*9<2*25 # sqrt(3/2)<5/3: dense cross inequality
json.dump({'field':'F3[u]/(u²+1)','parity_columns':cols,'sign_masks':signs,'rows':rows},open(P/'certificate.json','w'),indent=2)
print(json.dumps(rows,indent=2))
