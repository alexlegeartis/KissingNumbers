from itertools import combinations
from pathlib import Path
import numpy as np,json,sys
cells=range(2,32);blocks=[]
for a,b,c in combinations(cells,3):
 d=a^b^c
 if d>c and d<32:blocks.append((a,b,c,d))
assert len(blocks)==945
patterns=[[[0,1],[2,3]],[[0,2],[1,3]],[[0,3],[1,2]]]
supports=[]
for p in range(3):
 for block in blocks:
  rows=[x//2-1 for x in block]
  if p and len(set(rows))<4:continue
  m=0
  for x in block:
   row=x//2-1;side=x%2
   for c in patterns[p][side]:m|=1<<(4*row+c)
  assert m.bit_count()==8
  supports.append(m)
assert len(supports)==len(set(supports))==2625
assert max((a&b).bit_count() for i,a in enumerate(supports) for b in supports[:i])==4
assert all(sum(bool(m&(15<<(4*r))) for r in range(15))<=4 for m in supports)

Path(sys.argv[1] if len(sys.argv)>1 else 'bridge61_start.txt').write_text(''.join(str(x)+'\n' for x in supports))
print('PASS',len(supports),'new points',len(supports)*128,'new total',71310488+len(supports)*128)
