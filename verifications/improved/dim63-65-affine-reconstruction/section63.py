from pathlib import Path
import json,sys
from math import comb
from itertools import combinations
import numpy as np
H=Path(__file__).resolve().parent;P=H.parent/'dim61-67-compatible-sign-families';sys.path.insert(0,str(P));from subcode67 import kernel,span
S=P.parent/'dim61-68-95-96-code-layers/source'
f=json.loads((P/'data/global64.json').read_text());c=json.loads((H/'certificate.json').read_text())
rows={r['support']:r for r in f['old']+f['new']+c['modified_parent_classes']+c['new_classes']}
with np.load(S/'dim61_96/code64_supports.npz') as z:old=list(map(int,z['words']))
new=[r['support'] for r in f['new']+c['new_classes']]
ws={s:np.array([w^r['shift'] for w in span(kernel(r['checks']))],dtype=np.uint16) for s,r in rows.items()}
octads=json.loads((S/'shared_parity/dual_octads.json').read_text())
def count(s,h):
 m=(s&h).bit_count()
 if m%2:return 0
 if s not in rows:
  if m<8:return comb(m,m//2)*(1<<(11-m))
  assert m==8
  return 1120
 mask=sum(((h>>i)&1)<<j for j,i in enumerate(i for i in range(64) if s>>i&1))
 return int(np.count_nonzero(np.bitwise_count(ws[s]&mask)==m//2))
res=[]
for h in octads:
 terms=dict(dense=70*(1<<21),old_sparse=sum(count(s,h) for s in old),new_sparse=sum(count(s,h) for s in new),four_supports=102900,poles=112)
 res.append(dict(octad=h,terms=terms,total=sum(terms.values())))
best=max(res,key=lambda r:(r['total'],-r['octad']))
# This is a section of a valid relabeled variant, not necessarily the frozen
# 64D representative. Quadruple-support cross-layer bounds are label-independent.
h=best['octad']
perm=[i for i in range(64) if h>>i&1]+[i for i in range(64) if not(h>>i&1)]
blocks={tuple(sorted((a,b,c,a^b^c))) for a,b,c in combinations(range(64),3)}
assert len(blocks)==10416 and sorted(perm)==list(range(64))
four=0
for block in blocks:
    m=sum((h>>perm[i])&1 for i in block)
    if m%2==0:four+=comb(m,m//2)*(1<<(4-m))
assert four==102900==best['terms']['four_supports']
# For distinct quadruple supports intersection<=2; all sign vectors have norm64.
# Against weight16 amplitude2, dense amplitude1, and axis amplitude8,
# maximum scalar products are respectively32,16,32, regardless of labeling.
best['four_support_permutation']=perm
best['construction']='Relabel quadruple SQS independently; relabel unspecified full Hamming classes to put section coordinates first.'
out=dict(status='PASS',**best,gain=best['total']-184393988,octads=72,scope='Exact explicit affine-word count with published standard-class octad relabeling; requires full family geometry verification.')
(H/'section63.json').write_text(json.dumps(out,indent=2));(H/'all_sections.json').write_text(json.dumps(res,indent=2));print(out)
