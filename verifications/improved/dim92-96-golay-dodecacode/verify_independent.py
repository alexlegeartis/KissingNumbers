import json
from itertools import product,combinations
from collections import Counter
from pathlib import Path
p=Path(__file__).parent
g=[7,1,4,4,2,7,3,3,4,7,2,1]
BI=[tuple((1 if j==23 else int(i<=j<i+12 and g[j-i]%2)) for j in range(24))for i in range(12)]
def allcode(gen):
 C=[0]
 for r in gen:C += [x^r for x in C]
 assert len(C)==len(set(C))
 return C
binary=allcode([sum(a<<j for j,a in enumerate(r))for r in BI])
assert Counter(x.bit_count()for x in binary)=={0:1,8:759,12:2576,16:759,24:1}
# Represent quaternary symbols as separated planes of binary bits.
seq=[2,1,0,1,0,0,1,0,0,1,0,1]
qua=allcode([sum((seq[(j-i)%12]&1)<<j | (seq[(j-i)%12]>>1)<<(12+j) for j in range(12)) for i in range(12)])
wt=lambda u:((u&4095)|(u>>12)).bit_count()
assert Counter(map(wt,qua))=={0:1,6:396,8:1485,10:1980,12:234}
sign=[s for s in range(4096)if not any(sum((s>>j&1)*((j+1)>>k&1)for j in range(12))%2 for k in range(4))]
assert len(sign)==256 and min(x.bit_count() for x in sign if x)==3
points={}
for entry in json.loads((p/'patterns.json').read_text()):
 ids=[int(x,36)for x in entry['encoding'].split(',')]
 assert min(ids)>=0 and max(ids)<2**24 and len(set(ids))==len(ids)
 words=[(binary[z%4096],binary[z//4096])for z in ids]
 assert all((x|y).bit_count()==12 and not ((x|y)&entry['zero'])for x,y in words)
 R=Counter(((x^a)|(y^b)).bit_count()for (x,y),(a,b) in combinations(words,2))
 assert min(R)>=12
 points[entry['zero']]=len(words)*len(qua)*len(sign)
 print('zero mask:',entry['zero'],'patterns:',len(words),'size:',points[entry['zero']],'distance histogram:',R)
assert points=={15:324009984,3:983564288,0:2708471808}
