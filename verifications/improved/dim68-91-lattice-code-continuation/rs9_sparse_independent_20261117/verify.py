from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import json,numpy as np,hashlib
P=Path(__file__).resolve().parent;c=json.load(open(P/'input.json'))
# Independent exhaustive extended linear evaluation code, every degree<=5 polynomial.
add=np.array([[(a%3+b%3)%3+3*((a//3+b//3)%3)for b in range(9)]for a in range(9)],dtype=np.uint8)
mul=np.array([[((a%3)*(b%3)+2*(a//3)*(b//3))%3+3*(((a%3)*(b//3)+(a//3)*(b%3))%3)for b in range(9)]for a in range(9)],dtype=np.uint8)
for a in range(9):
 for b in range(9):
  for d in range(9):
   assert add[a,add[b,d]]==add[add[a,b],d]
   assert mul[a,mul[b,d]]==mul[mul[a,b],d]
   assert mul[a,add[b,d]]==add[mul[a,b],mul[a,d]]
assert all(any(mul[a,b]==1 for b in range(1,9))for a in range(1,9))
ids=np.arange(9**6,dtype=np.int64);coeff=np.stack([(ids//9**i)%9 for i in range(6)],axis=1).astype(np.uint8);words=np.empty((len(ids),10),dtype=np.uint8)
for t in range(9):
 v=coeff[:,5].copy()
 for j in range(4,-1,-1):v=add[mul[v,t],coeff[:,j]]
 words[:,t]=v
words[:,9]=coeff[:,5]
weights=np.count_nonzero(words,axis=1);assert weights[0]==0 and min(weights[1:])==5
cols=c['parity_columns'];assert len(cols)==10 and len(set(cols))==10 and all(0<x<16 for x in cols)
masks=[]
for v in range(1024):
 syn=0
 for i,x in enumerate(cols):
  if v>>i&1:syn^=x
 if syn==0:masks.append(v)
assert masks==c['sign_masks']and len(masks)==64 and all((v^1023)in masks for v in masks)
assert all((a^b).bit_count()>=3 for a,b in combinations(masks,2))
# (9+sqrt(3/2))²*(45/10)<(45/2)² iff sqrt(3/2)<5/3.
assert F(3,2)<F(25,9)
assert F(45**2,5)<F(45**2,4)and 135<F(45**2,4)
for r in c['rows']:
 n=r['dimension'];k=min(9,n-81);supports=int(np.count_nonzero(coeff[:,5]<k));assert supports==k*9**5==r['supports']
 assert r['dense']+64*supports+2*n*(n-1)==r['total']
out={'status':'PASS','input_sha256':hashlib.sha256((P/'input.json').read_bytes()).hexdigest(),'all_extended_codewords_checked':len(ids),'extended_minimum_distance':int(min(weights[1:])),'extended_weight_histogram':np.bincount(weights,minlength=11).tolist(),'sign_codewords':64,'rows':c['rows']}
(P/'report.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
