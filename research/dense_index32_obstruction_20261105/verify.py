"""Exact 32D subcode-coset obstruction; checks all six/five-anchor cases."""
import itertools,json,math,pathlib,numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[2];here=pathlib.Path(__file__).resolve().parent
source=json.loads((ROOT/'verifications/improved/dim32-35-signed-layers/source/classic/certificate32.json').read_text());parent=source['dense_generator'];H=json.loads((here/'seed.json').read_text())['basis'];assert len(H)==17
C=[0]
for g in parent:C += [x^g for x in C]
assert len(set(C))==1<<17 and min(x.bit_count() for x in C if x)>=8 and all(h in set(C) for h in H)
# Independently derive a full-rank parity check for the first twelve rows.
piv={}
for v0 in H:
 v=v0
 while v:
  i=v.bit_length()-1
  if i in piv:v^=piv[i]
  else:piv[i]=v;break
assert len(piv)==17
piv={}
for v0 in H[:12]:
 v=v0
 while v:
  i=v.bit_length()-1
  if i in piv:v^=piv[i]
  else:piv[i]=v;break
assert len(piv)==12
K=[]
for j in (j for j in range(32) if j not in piv):
 v=1<<j
 for i in sorted(piv):
  if (v&piv[i]).bit_count()%2:v|=1<<i
 K.append(v)
assert len(K)==20 and all((a&b).bit_count()%2==0 for a in H[:12] for b in K)
cols=np.array([sum(((K[i]>>j)&1)<<i for i in range(20)) for j in range(32)],np.uint32);D=np.full(1<<20,255,np.uint8);D[0]=0;front=np.array([0],np.uint32)
for depth in range(1,33):
 q=np.unique((front[:,None]^cols[None,:]).ravel());front=q[D[q]==255];D[front]=depth
 if not len(front):break
assert not np.any(D==255)
labels=[sum(((H[j]&K[i]).bit_count()%2)<<i for i in range(20)) for j in range(12,17)];O=[0]
for a in labels:O += [v^a for v in O]
O=np.array(O,np.uint32);assert len(set(map(int,O)))==32 and np.all(D[O[:,None]^O[None,:]][np.triu_indices(32,1)]>=8)
Q=np.flatnonzero(D>=8).astype(np.uint32);Q=Q[~np.isin(Q,O)];B=D[Q[:,None]^O[None,:]]<8;cost=B.sum(1);assert len(Q)==90057 and cost.min()==20
print('PASS quotient:',len(Q),'new cosets, each conflicts with at least 20 original cosets',flush=True)
# At least 33 selected cosets with k old ones have k<=12 and at least 21 new.
# If k>=6, normalize by one old coset and inspect a six-old anchor. If k=5,
# inspect the five-old anchor and require at least 28 new. Every other old
# coset remains available in the graph, making the core check conservative.
def check(anchor_size,threshold,required):
 counts={}
 for row in B:
  allowed=np.flatnonzero(~row).tolist()
  for anchor in itertools.combinations(allowed,anchor_size):counts[anchor]=counts.get(anchor,0)+1
 relevant=[(a,n) for a,n in counts.items() if n>=threshold];maxcore=0
 for ai,(anchor,num) in enumerate(relevant):
  elig=np.flatnonzero(~B[:,anchor].any(1));assert len(elig)==num
  V=np.r_[O[[j for j in range(32) if j not in anchor]],Q[elig]];good=D[V[:,None]^V[None,:]]>=8;np.fill_diagonal(good,False);alive=good.sum(1)>=required-1
  while alive.any():
   T=good[alive][:,alive];keep=T.sum(1)>=required-1
   if keep.all():break
   ids=np.flatnonzero(alive);alive[ids[~keep]]=False
  maxcore=max(maxcore,int(alive.sum()));assert not alive.any(),f'unresolved anchor {anchor}'
  if ai%5000==0:print('anchor size',anchor_size,'checked',ai,'of',len(relevant),flush=True)
 print('PASS anchor size',anchor_size,'relevant',len(relevant),'all cores empty',flush=True)
 return len(relevant)
assert check(6,21,27)==3954
assert check(5,28,28)==21051
print('PASS: any larger union containing an original coset must keep at most four original cosets; no global nonexistence claim.')
