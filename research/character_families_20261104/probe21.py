from pathlib import Path
import json,numpy as np
OUT=Path(__file__).resolve().parent
poly=sum(1<<i for i in (0,1,5,6,7,9,11))
G=[]
for m in range(4096):
 w=0
 for i in range(12):
  if m>>i&1:w^=poly<<i
 G.append(w | ((w.bit_count()&1)<<23))
D=[w&((1<<21)-1) for w in G]
line_ids=[i for i,w in enumerate(D) if w.bit_count()==5]
basis=[];piv={}
for m in line_ids:
 v=m
 while v:
  j=v.bit_length()-1
  if j in piv:v^=piv[j]
  else:piv[j]=v;basis.append(m);break
comp=[0]
for b in basis:comp += [m^b for m in comp]
pos={m:i for i,m in enumerate(comp)}
directions=[pos[m] for m in line_ids]
eta_rows=[(2,3),(1,10),(1,4,5,9),(0,2,7,13)]
def eta(w):return sum((sum(w>>j&1 for j in row)&1)<<i for i,row in enumerate(eta_rows))
labels={eta(D[m]):i for i,m in reversed(list(enumerate(line_ids)))}
seed=[labels[eta(D[m])] if eta(D[m]) else -1 for m in comp]
assert len(comp)==1024 and len(labels)==15 and 0 not in labels
for i,l in enumerate(seed):
 if l<0:continue
 for p,d in enumerate(directions):
  j=i^d;k=seed[j]
  assert k<0 or (l!=p and k!=p and l!=k)
obj={'code_polynomial_exponents':[0,1,5,6,7,9,11],'line_messages':line_ids,'component_basis':basis,'directions':directions,'seed':seed}
if __name__=='__main__':
 (OUT/'probe21.json').write_text(json.dumps(obj)+'\n')
 (OUT/'probe21.dat').write_text(f'{len(comp)} {len(line_ids)}\n'+' '.join(map(str,directions))+'\n'+' '.join(map(str,seed))+'\n')
 print({'vertices':len(comp),'labels':len(line_ids),'initial':sum(x>=0 for x in seed),'directions':directions})
